/*
 * 并网变流器实现细节（单位：kV/kA/MW/MVAr）：
 * - abc↔dq 变换采用电工角度，q 轴项取负号保证功率表达一致；
 * - PLL 路径：过零测频（窄带带通+窗口鲁棒）→ SRF-PLL；外环 PQ 采用精确逆映射投影到 (id,iq)；
 * - VSG 路径：摆动方程（H、D）更新 ω/θ，Q-V 下垂决定电压幅值；
 * - SPWM：统一三角载波 LUT；调制/电流限幅，并带反 wind-up；
 * - FPGA_simu：1 μs 小步 LC 恒导纳网络，稠密一次求逆后每子步做 Z*I。
 *
 * */
#include "converter.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>

static void abc_to_dq(double theta, double A, double B, double C, double& d, double& q)
{
    double c0 = std::cos(theta), s0 = std::sin(theta);
    double c1 = std::cos(theta - 2.0 * M_PI / 3.0), s1 = std::sin(theta - 2.0 * M_PI / 3.0);
    double c2 = std::cos(theta + 2.0 * M_PI / 3.0), s2 = std::sin(theta + 2.0 * M_PI / 3.0);
    d = (2.0 / 3.0) * (A * c0 + B * c1 + C * c2);
    q = -(2.0 / 3.0) * (A * s0 + B * s1 + C * s2); // 这里一定要有负号
}

static void dq_to_abc(double theta, double d, double q, double& A, double& B, double& C)
{
    // 逆变换对应同一套约定：q 项前是 "−"
    A = std::cos(theta) * d - std::sin(theta) * q;
    B = std::cos(theta - 2.0 * M_PI / 3.0) * d - std::sin(theta - 2.0 * M_PI / 3.0) * q;
    C = std::cos(theta + 2.0 * M_PI / 3.0) * d - std::sin(theta + 2.0 * M_PI / 3.0) * q;
}

// === SPWM：三相共用同一三角载波，幅值 1.0（输出范围[-1,1]） ===
static void Cal_signal(double simu_time, double fre,
    const std::vector<double>& tri_LUT, double tri_step,
    const double sA, const double sB, const double sC,
    int& IGBT_a, int& IGBT_b, int& IGBT_c)
{
    const double Tsw = 1.0 / fre;
    // 周期内采样索引（避免昂贵的 trig/分段计算）
    double t_in = std::fmod(simu_time, Tsw);
    if (t_in < 0)
        t_in += Tsw;
    int idx = static_cast<int>(t_in / std::max(1e-9, tri_step));
    if (tri_LUT.empty()) { // 兜底（不应发生,停止仿真）
        // IGBT_a = (0.0 < sA) ? 1 : 0;
        // IGBT_b = (0.0 < sB) ? 1 : 0;
        // IGBT_c = (0.0 < sC) ? 1 : 0;
        exit(0);
    }
    if (idx >= static_cast<int>(tri_LUT.size()))
        idx = static_cast<int>(tri_LUT.size()) - 1;
    const double y0 = tri_LUT[idx];

    IGBT_a = (y0 < sA) ? 1 : 0;
    IGBT_b = (y0 < sB) ? 1 : 0;
    IGBT_c = (y0 < sC) ? 1 : 0;
}

// ================= FPGA_simu：内部 LC 恒定导纳模型 =================
void FPGA_simu::buid_model()
{
    // ===== 生成本 FPGA 独立的三角波 LUT =====
    // 约定：小步长固定为 1 us；一个周期内的采样点 N = round(Tsw / 1us)
    tri_step_ = 1e-6; // 与 cal_step 的小步一致
    const double Tsw = 1.0 / fre;
    tri_len_ = std::max(1, (int)std::llround(Tsw / tri_step_));
    tri_LUT_.resize(tri_len_);

    // 与旧 Cal_signal 完全对齐的波形：幅值 ±1，0~T 线性段 0→+1→-1→0
    for (int k = 0; k < tri_len_; ++k) {
        double t = k * tri_step_;
        double y0;
        if (t <= 0.25 * Tsw)
            y0 = 4.0 * t / Tsw;
        else if (t <= 0.75 * Tsw)
            y0 = 2.0 - 4.0 * t / Tsw;
        else
            y0 = -4.0 + 4.0 * t / Tsw;
        tri_LUT_[k] = y0; // ∈[-1,1]
    }
    // ===== 构建 LC 恒定导纳模型 =====
    V.resize(N_);
    I.resize(N_);
    for (int i = 0; i < N_; ++i) {
        V[i] = 0.0;
        I[i] = 0.0;
    }
    const double step = 1e-6; // 1 us
    value1_RL = (1 - step * R_ / 2.0 / L_) / (1 + step * R_ / 2.0 / L_);
    value2_RL = 1.0 / (R_ + (2.0 * L_) / step);
    C_value = step / (2.0 * C_uF_ * 0.001 * 0.001); // 以“常导纳+历史电流”形式

    RL_value = 1.0 / value2_RL;

    T_vec.resize(line_num_);

    // --- 元件类型 ---
    T_vec[0].type = T_type_R; // AC串并等效
    T_vec[1].type = T_type_R;
    T_vec[2].type = T_type_R;
    T_vec[3].type = T_type_RL; // A相滤波支路
    T_vec[4].type = T_type_RL; // B
    T_vec[5].type = T_type_RL; // C
    T_vec[6].type = T_type_I; // A上桥臂
    T_vec[7].type = T_type_I; // A下桥臂
    T_vec[8].type = T_type_I; // B上桥臂
    T_vec[9].type = T_type_I; // B下桥臂
    T_vec[10].type = T_type_I; // C上桥臂
    T_vec[11].type = T_type_I; // C下桥臂
    T_vec[12].type = T_type_R; // DC等效内阻
    T_vec[13].type = T_type_C; // DC母线电容
    T_vec[14].type = T_type_R; // DC+ 保持电阻到地（Norton 并联电阻）
    T_vec[15].type = T_type_R; // DC- 保持电阻到地（Norton 并联电阻）

    // --- 命名 ---
    T_vec[0].name = "Source0_A";
    T_vec[1].name = "Source0_B";
    T_vec[2].name = "Source0_C";
    T_vec[3].name = "PV_RL_A";
    T_vec[4].name = "PV_RL_B";
    T_vec[5].name = "PV_RL_C";
    T_vec[6].name = "PV_A_up";
    T_vec[7].name = "PV_A_down";
    T_vec[8].name = "PV_B_up";
    T_vec[9].name = "PV_B_down";
    T_vec[10].name = "PV_C_up";
    T_vec[11].name = "PV_C_down";
    T_vec[12].name = "PV_Rs";
    T_vec[13].name = "PV_line_C";
    T_vec[14].name = "DC_hold_plus";
    T_vec[15].name = "DC_hold_minus";

    // --- 拓扑节点 ---
    T_vec[0].node_I = 2;
    T_vec[0].node_J = -1;
    T_vec[1].node_I = 0;
    T_vec[1].node_J = -1;
    T_vec[2].node_I = 1;
    T_vec[2].node_J = -1;

    T_vec[3].node_I = 2;
    T_vec[3].node_J = 3; // AC A 与滤波节点
    T_vec[4].node_I = 0;
    T_vec[4].node_J = 4;
    T_vec[5].node_I = 1;
    T_vec[5].node_J = 5;

    T_vec[6].node_I = 3;
    T_vec[6].node_J = 7; // A上臂：滤波点-DC+
    T_vec[7].node_I = 3;
    T_vec[7].node_J = 6; // A下臂：滤波点-DC-
    T_vec[8].node_I = 4;
    T_vec[8].node_J = 7;
    T_vec[9].node_I = 4;
    T_vec[9].node_J = 6;
    T_vec[10].node_I = 5;
    T_vec[10].node_J = 7;
    T_vec[11].node_I = 5;
    T_vec[11].node_J = 6;

    T_vec[12].node_I = 7;
    T_vec[12].node_J = 6; // DC等效内阻
    T_vec[13].node_I = 7;
    T_vec[13].node_J = 6; // DC电容

    T_vec[14].node_I = 7;
    T_vec[14].node_J = -1; // DC+ 到地
    T_vec[15].node_I = 6;
    T_vec[15].node_J = -1; // DC- 到地

    // --- 互锁关系（用于历史项切换） ---
    T_vec[6].next_IGBT = 7; // A_up   <-> A_down
    T_vec[7].next_IGBT = 6;
    T_vec[8].next_IGBT = 9; // B_up   <-> B_down
    T_vec[9].next_IGBT = 8;
    T_vec[10].next_IGBT = 11; // C_up   <-> C_down
    T_vec[11].next_IGBT = 10;

    // --- 数值 ---
    T_vec[0].value = R_AC_;
    T_vec[1].value = R_AC_;
    T_vec[2].value = R_AC_;
    T_vec[3].value = RL_value;
    T_vec[4].value = RL_value;
    T_vec[5].value = RL_value;
    T_vec[6].value = 1.0;
    T_vec[7].value = 1.0;
    T_vec[8].value = 1.0;
    T_vec[9].value = 1.0;
    T_vec[10].value = 1.0;
    T_vec[11].value = 1.0;
    T_vec[12].value = 1e9;
    T_vec[13].value = C_value;
    T_vec[14].value = R_hold_;
    T_vec[15].value = R_hold_;

    // --- 构建导纳矩阵（稠密） ---
    Y_mat = Eigen::MatrixXd::Zero(N_, N_);
    for (int i = 0; i < line_num_; ++i) {
        int node_i = T_vec[i].node_I;
        int node_j = T_vec[i].node_J;
        double data = 1.0 / T_vec[i].value;
        if (node_j <= 0) {
            Y_mat(node_i, node_i) += data;
        } else {
            Y_mat(node_i, node_i) += data;
            Y_mat(node_j, node_j) += data;
            Y_mat(node_i, node_j) -= data;
            Y_mat(node_j, node_i) -= data;
        }
    }
    Y_mat_inv = Y_mat.inverse();
}

void FPGA_simu::cal_step(send_data& send_, receive_data& receive_)
{
    const int subSteps = 50;
    const double step = 1e-6;

    double Va = send_.Va;
    double Vb = send_.Vb;
    double Vc = send_.Vc;
    double Idc = send_.Idc;

    for (int iter = 0; iter < subSteps; ++iter) {
        I.setZero();
        simu_time_ += step;

        int Sa = 0, Sb = 0, Sc = 0;
        // === 查表版三角波比较 ===
        Cal_signal(simu_time_, fre, tri_LUT_, tri_step_,
            send_.Vdoor_va, send_.Vdoor_vb, send_.Vdoor_vc, Sa, Sb, Sc);
        std::array<int, 6> s_vec { Sa, 1 - Sa, Sb, 1 - Sb, Sc, 1 - Sc };

        // 交流端口等效注入（单位自洽：kV/Ω -> kA）
        I[2] = Va / R_AC_;
        I[0] = Vb / R_AC_;
        I[1] = Vc / R_AC_;
        // 直流侧电流源（Idc>0 表示 DC+ 向系统注入）
        I[6] = -Idc;
        I[7] = Idc;

        const double Vhalf = 0.5 * Vdc_target_kV_;
        I[7] += Vhalf / R_hold_; // 让 V7 ≈ +Vhalf
        I[6] += -Vhalf / R_hold_; // 让 V6 ≈ -Vhalf

        // 其余节点初始化
        I[3] = I[4] = I[5] = 0.0;

        int cc { 0 };
        for (int line_cc = 0; line_cc < line_num_; ++line_cc) {
            double value = 1.0 / T_vec[line_cc].value;

            // 支路端电压（kV）
            if (T_vec[line_cc].node_J < 0) {
                T_vec[line_cc].hisV = V[T_vec[line_cc].node_I];
            } else {
                T_vec[line_cc].hisV = V[T_vec[line_cc].node_I] - V[T_vec[line_cc].node_J];
            }

            if (T_vec[line_cc].type == T_type_R) {
                T_vec[line_cc].hisI = T_vec[line_cc].value * T_vec[line_cc].hisV;
            }
            if (T_vec[line_cc].type == T_type_C) {
                T_vec[line_cc].hisI = T_vec[line_cc].hist + value * T_vec[line_cc].hisV;
                T_vec[line_cc].hist = -T_vec[line_cc].hisI - value * T_vec[line_cc].hisV;
                I[T_vec[line_cc].node_I] -= T_vec[line_cc].hist;
                I[T_vec[line_cc].node_J] += T_vec[line_cc].hist;
            }
            if (T_vec[line_cc].type == T_type_RL) {
                T_vec[line_cc].hisI = T_vec[line_cc].hist + value * T_vec[line_cc].hisV;
                T_vec[line_cc].hist = T_vec[line_cc].hisI * value1_RL + T_vec[line_cc].hisV * value2_RL;
                I[T_vec[line_cc].node_I] -= T_vec[line_cc].hist;
                I[T_vec[line_cc].node_J] += T_vec[line_cc].hist;
            }
            if (T_vec[line_cc].type == T_type_I) {
                int c = s_vec[cc++];
                int next_IGBT = T_vec[line_cc].next_IGBT;
                T_vec[line_cc].hisI = T_vec[line_cc].hist + value * T_vec[line_cc].hisV;

                if (c == 1 && T_vec[line_cc].old_c == 1) {
                    T_vec[line_cc].hist = T_vec[line_cc].hisI + T_vec[line_cc].hisI * alpha_;
                }
                if (c == 0 && T_vec[line_cc].old_c == 0) {
                    T_vec[line_cc].hist = T_vec[line_cc].hisI * belta_ - T_vec[line_cc].hisV;
                }
                if (c == 1 && T_vec[line_cc].old_c == 0 && iter >= 1) {
                    T_vec[line_cc].hist = -T_vec[next_IGBT].hisI + T_vec[next_IGBT].hisV * alpha_;
                }
                if (c == 0 && T_vec[line_cc].old_c == 1 && iter >= 1) {
                    T_vec[line_cc].hist = -T_vec[next_IGBT].hisV * belta_ - T_vec[next_IGBT].hisV;
                }

                I[T_vec[line_cc].node_I] -= T_vec[line_cc].hist;
                I[T_vec[line_cc].node_J] += T_vec[line_cc].hist;
                T_vec[line_cc].old_c = c;

                const double I_HIST_LIM = 50; // kA 级别足够大，防数值飞天
                T_vec[line_cc].hist = std::clamp(T_vec[line_cc].hist, -I_HIST_LIM, I_HIST_LIM);
            }
        }
        V = Y_mat_inv * I;
    }

    // 输出给外部（单位：kA/kV）
    receive_.Ia = -T_vec[3].hisI;
    receive_.Ib = -T_vec[4].hisI;
    receive_.Ic = -T_vec[5].hisI;
    receive_.Vapp = V[3];
    receive_.Vbpp = V[4];
    receive_.Vcpp = V[5];
    receive_.Vdc1 = V[6];
    receive_.Vdc2 = V[7];
}

// ================= Converter：控制与并网接口 =================
Converter::Converter(int na, int nb, int nc, double P_ref_MW, double Q_ref_MVAr)
    : nA_(na)
    , nB_(nb)
    , nC_(nc)
    , Pref_MW_(P_ref_MW)
    , Qref_MVAr_(Q_ref_MVAr)
{
    // FPGA 模型初始化
    fpga_.buid_model();
    Vdc_kV_ = fpga_.Vdc_target_kV_;
}

static inline void abc_to_alphabeta(double A, double B, double C,
    double& alpha, double& beta)
{
    // 三相对称、无零序假设
    alpha = (2.0 / 3.0) * (A - 0.5 * B - 0.5 * C);
    beta = (1.0 / std::sqrt(3.0)) * (B - C); // 等价于 (2/3)*(√3/2)*(B - C)
}

void Converter::pll_updateHistory(Eigen::VectorXd& I, double t, double dt)
{
    // ---------------- 测频、dq 测量与一阶低通 ----------------
    freq_e2.step(t, vA_kV_, vB_kV_, vC_kV_);

    double ud, uq, id_meas, iq_meas;
    abc_to_dq(theta_, vA_kV_, vB_kV_, vC_kV_, ud, uq);
    abc_to_dq(theta_, iA_kA_, iB_kA_, iC_kA_, id_meas, iq_meas);

    const double a_f = dt / T_filter_;
    ud_f_ += a_f * (ud - ud_f_);
    uq_f_ += a_f * (uq - uq_f_);
    id_f_ += a_f * (id_meas - id_f_);
    iq_f_ += a_f * (iq_meas - iq_f_);

    // ----------------  SRF-PLL ----------------
    double vmag = std::hypot(ud_f_, uq_f_);
    if (vmag < 1e-6)
        vmag = 1e-6;
    const double e_pll = uq_f_ / vmag;

    constexpr double W_INT_MAX = 2.0 * PI * 10.0;
    pll_int_ = std::clamp(pll_int_ + Ki_pll_ * e_pll * dt, -W_INT_MAX, W_INT_MAX);
    double delta_w = Kp_pll_ * e_pll + pll_int_;
    constexpr double W_MIN = 2.0 * PI * 45.0, W_MAX = 2.0 * PI * 55.0;
    omega_ = std::clamp(omega_nom_ + delta_w, W_MIN, W_MAX);
    theta_ = std::remainder(theta_ + omega_ * dt, 2.0 * M_PI);

    // ----------------  功率计算与一阶低通 ----------------
    const double P_inst = 1.5 * (ud_f_ * id_f_ + uq_f_ * iq_f_);
    const double Q_inst = 1.5 * (ud_f_ * iq_f_ - uq_f_ * id_f_);
    P_inst_MW_ = P_inst;
    Q_inst_MVAr_ = Q_inst;
    P_f_ += (P_inst - P_f_) * (dt / Tpq_);
    Q_f_ += (Q_inst - Q_f_) * (dt / Tpq_);

    // ---------------- 目标 P/Q（含软启与下垂） ----------------
    const double U_mag = std::hypot(ud_f_, uq_f_);
    double kk = (t < 1.0) ? t : 1.0; // 0~1 s 软启
    double Pref_cmd = kk * Pref_MW_;
    double Qref_cmd = kk * Qref_MVAr_;

    if (droop_enabled_ && t > 1.0) {
        const double f_now = freq_e2.freq_hz(), f_N = 50.0;
        if (std::abs(f_now - f_N) > 0.02) {
            const double Kf_ = 50.0;
            Pref_cmd += -Kf_ * (f_now - f_N) / f_N * S_nom;
        }
        const double dV = U_mag * std::sqrt(3.0) / std::sqrt(2.0) - V_LL_kV_const_;
        if (std::abs(dV) > 0.06 * V_LL_kV_const_) {
            const double Kv_ = 0.6;
            Qref_cmd += -Kv_ * dV / V_LL_kV_const_ * S_nom;
        }
    }

    // ----------------  外环：用同一“逆矩阵”做前馈 + 误差投影 ----------------
    // 逆矩阵:  [id; iq] = 1/(1.5*(ud^2+uq^2)) * [[ ud, -uq],[ uq, ud]] * [P; Q]
    const double det = std::max(1e-6, ud_f_ * ud_f_ + uq_f_ * uq_f_);
    const double inv = 1.0 / (1.5 * det);

    // —— 前馈 ——（完全消除 θ 偏差带来的初值误差）
    const double id_ff = inv * (ud_f_ * Pref_cmd - uq_f_ * Qref_cmd);
    const double iq_ff = inv * (uq_f_ * Pref_cmd + ud_f_ * Qref_cmd);

    // —— 误差（基于低通功率）——
    const double eP = Pref_cmd - P_f_;
    const double eQ = Qref_cmd - Q_f_;

    // 外环 PI（先在功率域积分，再统一投影到电流域）
    eP_int_ += KiP_base_ * eP * dt;
    eQ_int_ += KiQ_base_ * eQ * dt;

    const double P_ctrl = KpP_base_ * eP + eP_int_;
    const double Q_ctrl = KpQ_base_ * eQ + eQ_int_;

    // 把 PI 输出来到 (id, iq)（与前馈同一逆矩阵）
    double id_star = id_ff + inv * (ud_f_ * P_ctrl - uq_f_ * Q_ctrl);
    double iq_star = iq_ff + inv * (uq_f_ * P_ctrl + ud_f_ * Q_ctrl);

    // ---------------- 6) 电流限幅 + 外环抗风-up ----------------
    const double I_rated_kA = 1.5 * S_nom / std::max(0.05, std::sqrt(3.0) * V_LL_kV_const_);
    const double I_max = std::max(0.01, i_limit_factor_ * I_rated_kA);
    const double i_mag = std::hypot(id_star, iq_star);

    saturated_ = (i_mag > I_max);
    if (saturated_) {
        const double s = I_max / i_mag;
        id_star *= s;
        iq_star *= s;
        // 回退本拍外环积分（功率域）以抗风-up
        eP_int_ -= KiP_base_ * eP * dt;
        eQ_int_ -= KiQ_base_ * eQ * dt;
        ++sat_hits_;
    }

    // ---------------- 7) 内环 PI（带解耦） ----------------
    const double ed = id_star - id_f_;
    const double eq = iq_star - iq_f_;

    if (!saturated_) { // 饱和时冻结内环积分
        id_int_ += Ki_i_ * ed * dt;
        iq_int_ += Ki_i_ * eq * dt;
    }

    const double wL = omega_ * fpga_.L_;
    const double vd_star = Kp_i_ * ed + id_int_ - wL * iq_f_ + ud_f_;
    const double vq_star = Kp_i_ * eq + iq_int_ + wL * id_f_ + uq_f_;

    // ---------------- abc、调制与饱和缩放 ----------------
    double va_ref, vb_ref, vc_ref;
    dq_to_abc(theta_, vd_star, vq_star, va_ref, vb_ref, vc_ref);

    const double vdc_half = std::max(0.05, 0.5 * Vdc_kV_);
    double mA = va_ref / vdc_half, mB = vb_ref / vdc_half, mC = vc_ref / vdc_half;
    const double m_peak = std::max({ std::fabs(mA), std::fabs(mB), std::fabs(mC) });
    constexpr double M_MAX = 0.95;
    if (m_peak > M_MAX) {
        const double s = M_MAX / m_peak;
        mA *= s;
        mB *= s;
        mC *= s;
        // 轻度反风-up：撤销本拍内环积分增长
        id_int_ -= Ki_i_ * ed * dt;
        iq_int_ -= Ki_i_ * eq * dt;
    }
    send_.Vdoor_va = std::clamp(mA, -1.0, 1.0);
    send_.Vdoor_vb = std::clamp(mB, -1.0, 1.0);
    send_.Vdoor_vc = std::clamp(mC, -1.0, 1.0);

    // ----------------  直流侧等效电流源 ----------------
    const double Idc_cmd_kA = Pref_cmd / std::max(0.1, Vdc_kV_);
    send_.Va = vA_kV_;
    send_.Vb = vB_kV_;
    send_.Vc = vC_kV_;
    send_.Idc = Idc_cmd_kA;

    // ----------------  调用“FPGA”小步仿真并回读 ----------------
    fpga_.cal_step(send_, receive_);
    iA_kA_ = receive_.Ia;
    iB_kA_ = receive_.Ib;
    iC_kA_ = receive_.Ic;

    // 慢速更新 Vdc 估计（用于下拍调制归一化/保护）
    double Vdc_meas = std::max(0.05, receive_.Vdc2 - receive_.Vdc1);
    Vdc_kV_ = 0.95 * Vdc_kV_ + 0.05 * Vdc_meas;

    // ---------------- 11) 注入到网络电流向量 ----------------
    if (nA_ > 0)
        I(nA_ - 1) += iA_kA_;
    if (nB_ > 0)
        I(nB_ - 1) += iB_kA_;
    if (nC_ > 0)
        I(nC_ - 1) += iC_kA_;
}

void Converter::vsg_updateHistory(Eigen::VectorXd& I, double t, double dt)
{
    freq_e2.step(t, vA_kV_, vB_kV_, vC_kV_);
    if (!vsg_theta_inited) {
        double alpha, beta;
        abc_to_alphabeta(vA_kV_, vB_kV_, vC_kV_, alpha, beta);
        if (std::abs(alpha) + std::abs(beta) > 1e-6) {
            theta_ = std::atan2(beta, alpha); // 让 θ 与并网电压对齐
            omega_ = omega_nom_; // 初速=额定
            // 初始化测量滤波器，避免第一拍突变
            double ud0, uq0;
            abc_to_dq(theta_, vA_kV_, vB_kV_, vC_kV_, ud0, uq0);
            ud_f_ = ud0;
            uq_f_ = uq0;
            id_f_ = iq_f_ = 0.0;
            vsg_theta_inited = true;
        }
    }
    //  dq 变换（使用自身 θ，不依赖 PLL）
    double ud, uq, id_meas, iq_meas;
    abc_to_dq(theta_, vA_kV_, vB_kV_, vC_kV_, ud, uq);
    abc_to_dq(theta_, iA_kA_, iB_kA_, iC_kA_, id_meas, iq_meas);

    //  与原代码一致的测量低通
    double alpha_filter = dt / T_filter_;
    ud_f_ += alpha_filter * (ud - ud_f_);
    uq_f_ += alpha_filter * (uq - uq_f_);
    id_f_ += alpha_filter * (id_meas - id_f_);
    iq_f_ += alpha_filter * (iq_meas - iq_f_);

    double Pref_cmd = vsg_.P0_MW; // 基准有功（和 PLL 的 Pref_MW_ 等价）
    double Qref_cmd = vsg_.Q0_MVAr; // 基准无功（和 PLL 的 Qref_MVAr_ 等价）

    // 启动软升：与 PLL 的 kk=t 保持一致
    if (t < 1.0) {
        Pref_cmd *= t;
        Qref_cmd *= t;
    }
    double U_mag = std::hypot(ud_f_, uq_f_);

    if (droop_enabled_ && t > 1.0) {
        // 频率用本机 ω（VSG 不依赖 PLL）
        double f_now = freq_e2.freq_hz();
        const double f_N = 50.0;

        // —— P-f 下垂：与 PLL 相同 ——
        if (std::fabs(f_now - f_N) > 0.02) {
            const double Kf_ = 50.0; // 与 PLL 相同
            const double dP = -Kf_ * (f_now - f_N) / f_N * S_nom; // MW
            Pref_cmd += dP;
        }

        // —— Q-V 下垂：与 PLL 相同 ——
        // PLL 里是：U_mag*(√3/√2) 与 V_LL_kV_const_ 做比较

        const double V_ll = U_mag * std::sqrt(3.0) / std::sqrt(2.0);
        const double dV = -V_ll + V_LL_kV_const_;
        if (std::fabs(dV) > 0.06 * V_LL_kV_const_) {

            double Kv_ = 0.6; // 调压系数
            double dQ = Kv_ * dV / V_LL_kV_const_ * S_nom; // MVAr
            Qref_cmd += dQ;
        }
    }

    // 瞬时功率 & 一阶低通（沿用 Tpq_）
    double P = 1.5 * (ud_f_ * id_f_ + uq_f_ * iq_f_);
    double Q = 1.5 * (ud_f_ * iq_f_ - uq_f_ * id_f_);
    P_inst_MW_ = P;
    Q_inst_MVAr_ = Q;
    P_f_ += (P - P_f_) * (dt / Tpq_);
    Q_f_ += (Q - Q_f_) * (dt / Tpq_);

    // 摆动方程：只积分一次 ω，再更新一次 θ
    const double w0 = 2.0 * M_PI * 50.0;
    double Pe_pu = P_f_ / std::max(1e-6, S_nom);
    double Pm_pu = Pref_cmd / std::max(1e-6, S_nom);
    double domega_pu = (Pm_pu - Pe_pu - vsg_.D_pu * (omega_ / w0 - 1.0))
        / (2.0 * std::max(1e-6, vsg_.H_s));
    omega_ += w0 * domega_pu * dt;
    omega_ = std::clamp(omega_, 2.0 * M_PI * vsg_.fmin_Hz, 2.0 * M_PI * vsg_.fmax_Hz);
    theta_ = std::remainder(theta_ + omega_ * dt, 2.0 * M_PI);

    //  无功-电压静态下垂 -> 电压幅值（pu）
    U_mag_f_ = 0.98 * U_mag_f_ + 0.02 * std::max(0.05, U_mag); // 平滑测量

    double k_ramp = (t < 1.0) ? t : 1.0; // 0~1 s 软升压
    double Ed_kV = U_mag_f_ * k_ramp; // d 轴内电势
    double Eq_kV = 0.0;

    // ——  P/Q -> i* 精确映射（不依赖角度完全对齐）——
    // [ P ]   = 1.5 * [  ud   uq ] [ id ]
    // [ Q ]             [ -uq   ud ] [ iq ]
    // 反解：
    // [ id ] = (1/(1.5*(ud^2+uq^2))) * [  ud  -uq ] [ P ]
    // [ iq ]                             [  uq   ud ] [ Q ]
    double det = ud_f_ * ud_f_ + uq_f_ * uq_f_;
    double inv = 1.0 / (1.5 * std::max(det, 1e-6));
    double id_star = inv * (ud_f_ * Pref_cmd - uq_f_ * Qref_cmd);
    double iq_star = inv * (uq_f_ * Pref_cmd + ud_f_ * Qref_cmd);

    // ----  电流限幅（同 PLL 路径），并做反风up ----
    double I_rated_kA = 1.5 * S_nom / std::max(0.05, std::sqrt(3.0) * V_LL_kV_const_);
    double I_max = std::max(0.01, i_limit_factor_ * I_rated_kA);
    double i_mag = std::hypot(id_star, iq_star);
    bool limited = (i_mag > I_max);
    if (limited) {
        double s = I_max / i_mag;
        id_star *= s;
        iq_star *= s;
    }

    // ----  电流 PI（带解耦），生成 vd*/vq* ----
    double ed = id_star - id_f_, eq = iq_star - iq_f_;
    if (!limited) {
        id_int_ += Ki_i_ * ed * dt;
        iq_int_ += Ki_i_ * eq * dt;
    } // 饱和冻结
    double wL = omega_ * fpga_.L_;
    double vd_star = Kp_i_ * ed + id_int_ - wL * iq_f_ + ud_f_;
    double vq_star = Kp_i_ * eq + iq_int_ + wL * id_f_ + uq_f_;

    // ----  abc、调制和抗饱和缩放 ----
    double va_ref, vb_ref, vc_ref;
    dq_to_abc(theta_, vd_star, vq_star, va_ref, vb_ref, vc_ref);

    double vdc_half = std::max(0.05, Vdc_kV_ * 0.5);
    double mA = va_ref / vdc_half, mB = vb_ref / vdc_half, mC = vc_ref / vdc_half;
    double m_peak = std::max({ std::fabs(mA), std::fabs(mB), std::fabs(mC) });
    const double M_MAX = 0.95; // 调制上限：留点裕量
    if (m_peak > M_MAX) {
        double s = M_MAX / m_peak;
        mA *= s;
        mB *= s;
        mC *= s;
        // 反风up：撤销本次积分增长
        id_int_ -= Ki_i_ * ed * dt;
        iq_int_ -= Ki_i_ * eq * dt;
    }
    send_.Vdoor_va = std::clamp(mA, -1.0, 1.0);
    send_.Vdoor_vb = std::clamp(mB, -1.0, 1.0);
    send_.Vdoor_vc = std::clamp(mC, -1.0, 1.0);

    // 直流侧功率：Idc = Pm / Vdc
    double Idc_cmd_kA = Pref_cmd / std::max(0.1, Vdc_kV_);
    send_.Va = vA_kV_;
    send_.Vb = vB_kV_;
    send_.Vc = vC_kV_;
    send_.Idc = Idc_cmd_kA;

    //  小步电路
    fpga_.cal_step(send_, receive_);
    iA_kA_ = receive_.Ia;
    iB_kA_ = receive_.Ib;
    iC_kA_ = receive_.Ic;

    // 更新 Vdc 估计
    double Vdc_meas = std::max(0.05, receive_.Vdc2 - receive_.Vdc1);
    Vdc_kV_ = 0.95 * Vdc_kV_ + 0.05 * Vdc_meas;

    //  将注入电流并到网络 I 向量
    if (nA_ > 0)
        I(nA_ - 1) += iA_kA_;
    if (nB_ > 0)
        I(nB_ - 1) += iB_kA_;
    if (nC_ > 0)
        I(nC_ - 1) += iC_kA_;
}

void Converter::updateHistory(Eigen::VectorXd& I, double t, double dt)
{
    if (mode_ == ControlMode::PLL) {
        pll_updateHistory(I, t, dt);
    } else {
        vsg_updateHistory(I, t, dt);
    }
}

void Converter::updateState(const Eigen::VectorXd& V, double dt)
{
    // 采样并网点电压（两种模式通用）
    vA_kV_ = (nA_ > 0) ? V(nA_ - 1) : 0.0;
    vB_kV_ = (nB_ > 0) ? V(nB_ - 1) : 0.0;
    vC_kV_ = (nC_ > 0) ? V(nC_ - 1) : 0.0;
}
