/* =============================================================================
 * 并网变流器（Converter）：PLL 并网型 / VSG 电网成形型
 * -----------------------------------------------------------------------------
 * - PLL 路径：SRF-PLL → 外环 PQ（精确前馈 + PI 投影）→ 内环电流 PI（解耦）→ SPWM
 * - VSG 路径：摆动方程（H, D） + P–f / Q–V 下垂；同样通过内环整形到可调制电压
 * - FPGA 子模块：1 μs 等效小步 LC+开关恒导纳网络，用于系统级耦合（非器件级开关仿真）
 * - 输出侧对外等效为注入型诺顿电流源（在 updateHistory() 向 I 向量注入 Ia/Ib/Ic）
 * ============================================================================= */

#pragma once

#include "Devices.h"
#include <Eigen/Dense>
#include <algorithm>
#include <array>
#include <bitset>
#include <deque>
#include <string>
#include <vector>

using std::bitset;
using std::string;
using std::vector;

static constexpr double PI = 3.14159265358979323846;

#pragma once
#include <algorithm>
#include <cfloat>
#include <cmath>
#include <deque>
#include <vector>

enum class ControlMode {
    PLL, // 基于PLL的控制
    VSG // 虚拟同步机控制
};
// 目前仅支持两种模式，默认为 PLL,还有VSG方式

// 过零频率测量（基于过零时刻的周期窗口法）
// 适合工频信号，抗谐波和噪声能力较强
// 需要信号有明显的过零点（即不适合直流偏置过大的情况）
struct ZeroCrossFreqWin {
    // ---------- 可调参数 ----------
    double f0_hz = 50.0; // 初始频率
    double f_min = 45.0, f_max = 55.0; // 合理物理范围
    double hyst_ratio = 0.02; // 自适应滞回 = 2%*Vpk
    int win_size = 10; // 周期窗口长度（单位：周期）
    double trim_ratio = 0.2; // 截尾比例
    double tau_ema = 0.05; // 输出 EMA 时间常数≈50ms

    // —— 预滤波（二阶带通，中心频率随估计频率跟踪）——
    bool prefilter_on = true; // 开关
    double bp_Q = 20.0; // 带通品质因数（15~30 较合适）
    double fc_track_alpha = 0.15; // 中心频率跟踪速度（0~1，小=更稳）

    // —— 新增：窗口清空策略（避免频繁清空）——
    double flush_thresh_hz = 0.02; // 单周期频率与 f_ema 差超阈值则清空
    int flush_hold_cycles = 10; // 清空后至少等待 N 个有效周期才再次清空

    // ---------- 内部状态 ----------
    // 自适应幅值估计（用于滞回阈值）
    double vrms2 = 0.0;
    double beta_env = 0.02;

    // 上一次采样
    double t_prev = 0.0, a_prev = 0.0;
    bool have_prev = false;

    // 过零武装/触发（上升沿）
    bool armed_rise = false;

    // 上次过零时刻
    double t_zc_prev = NAN;
    std::deque<double> Ts;

    // 频率输出
    double f_robust = 50.0;
    double f_ema = 50.0;

    // —— 带通滤波器状态（双二阶/双延迟线）——
    // y[n] = (b0*x + b1*x1 + b2*x2 - a1*y1 - a2*y2)/a0
    double b0 = 0, b1 = 0, b2 = 0, a0 = 1, a1 = 0, a2 = 0;
    double x1 = 0, x2 = 0, y1 = 0, y2 = 0;
    double fc_used = 50.0; // 正在使用的中心频率
    double last_dt = 0.0;
    bool biq_init = false;

    // —— 新增：清空冷却计数（按有效过零周期计数）——
    int flush_cooldown_left = 0;

    void reset(double t0 = 0.0, double f_init = 50.0)
    {
        f0_hz = f_init;
        f_robust = f_ema = f_init;
        t_prev = t0;
        a_prev = 0.0;
        have_prev = false;
        armed_rise = false;
        t_zc_prev = NAN;
        Ts.clear();
        vrms2 = 0.0;

        x1 = x2 = y1 = y2 = 0.0;
        fc_used = f_init;
        biq_init = false;
        last_dt = 0.0;

        flush_cooldown_left = 0;
    }

    // 每步调用；t 为绝对仿真时间（秒）
    void step(double t, double va, double vb, double vc)
    {
        double dt = have_prev ? std::max(1e-6, t - t_prev) : 1e-4;

        // 1) Clarke α（无零序）
        const double alpha_raw = (2.0 / 3.0) * (va - 0.5 * vb - 0.5 * vc);

        // 2) 可选：窄带带通（抑制谐波/载波），中心频率跟踪 f_ema
        double alpha = alpha_raw;
        if (prefilter_on) {
            // 跟踪中心频率（慢跟随，避免频繁换系数导致抖动）
            double fc_target = std::clamp(f_ema, f_min, f_max);
            fc_used += fc_track_alpha * (fc_target - fc_used);

            // 如首次/采样时间变化大/频率变动较多 => 重新计算系数
            if (!biq_init || std::fabs(dt - last_dt) > 1e-9 || std::fabs(fc_used - fc_target) > 0.2) {
                update_biquad(dt, fc_used, bp_Q);
                biq_init = true;
                last_dt = dt;
            }

            // 二阶带通
            double yn = (b0 * alpha_raw + b1 * x1 + b2 * x2 - a1 * y1 - a2 * y2) / a0;
            x2 = x1;
            x1 = alpha_raw;
            y2 = y1;
            y1 = yn;
            alpha = yn; // 用滤波后的 α 做后续逻辑
        }

        // 3) 自适应幅值估计（基于滤波后的 α）
        vrms2 = (1.0 - beta_env) * vrms2 + beta_env * (alpha * alpha);
        const double vpk = std::sqrt(std::max(1e-12, 2.0 * vrms2));
        const double H = hyst_ratio * vpk;

        // 4) 武装条件：信号充分低（< -H）才允许下一次“上升沿过零”
        if (alpha < -H)
            armed_rise = true;

        // 5) 上升沿过零 + 线性插值
        if (have_prev) {
            if (armed_rise && (a_prev <= 0.0) && (alpha > 0.0)) {
                const double s = a_prev / (a_prev - alpha + 1e-18);
                const double t_zc = t_prev + s * dt;

                if (!std::isnan(t_zc_prev)) {
                    double T = t_zc - t_zc_prev;
                    const double Tmin = 1.0 / f_max, Tmax = 1.0 / f_min;

                    if (T >= Tmin && T <= Tmax) {
                        // —— 计算最新单周期频率 ——
                        const double f_last = 1.0 / T;

                        // —— 窗口清空条件（带冷却）——
                        if (flush_cooldown_left == 0 && std::fabs(f_last - f_ema) > flush_thresh_hz) {
                            // 清空窗口并用 f_last 立即刷新鲁棒输出
                            Ts.clear();
                            f_robust = std::clamp(f_last, f_min, f_max);

                            // 进入冷却：未来 N 个有效周期不再清空
                            flush_cooldown_left = flush_hold_cycles;
                        } else {
                            // —— 正常路径：窗口更新 + 鲁棒估计 ——
                            bool ok = true;
                            if (!Ts.empty()) {
                                double med = median_period();
                                if (med > 0) {
                                    double rel = std::fabs(T - med) / med;
                                    double tol = (t < 0.3) ? 0.05 : 0.02; // 启动放宽5%，之后2%
                                    if (rel > tol)
                                        ok = false;
                                }
                            }
                            if (ok) {
                                if ((int)Ts.size() == win_size)
                                    Ts.pop_front();
                                Ts.push_back(T);
                            }
                            if (!Ts.empty())
                                f_robust = 1.0 / robust_T();
                        }

                        // 冷却计数递减：仅在捕获到“有效过零周期”时递减
                        if (flush_cooldown_left > 0)
                            --flush_cooldown_left;
                    }
                }
                t_zc_prev = t_zc;
                armed_rise = false;
            }
        }

        // 6) 输出 EMA 平滑
        double a = dt / std::max(1e-6, tau_ema);
        f_ema += a * (f_robust - f_ema);

        // 7) 更新上次样本
        a_prev = alpha;
        t_prev = t;
        have_prev = true;
    }

    // 读取频率（Hz）
    double freq_hz() const { return f_ema; }

private:
    // 计算带通系数（RBJ cookbook：bandpass constant skirt gain）
    void update_biquad(double dt, double fc, double Q)
    {
        const double Fs = 1.0 / std::max(1e-9, dt);
        const double w0 = 2.0 * M_PI * fc / Fs;
        const double cw0 = std::cos(w0);
        const double sw0 = std::sin(w0);
        const double alpha = sw0 / (2.0 * std::max(1e-3, Q));

        const double a0_ = 1.0 + alpha;
        b0 = alpha / a0_;
        b1 = 0.0 / a0_;
        b2 = -alpha / a0_;
        a0 = 1.0; // 归一化后 a0 = 1
        a1 = (-2.0 * cw0) / a0_;
        a2 = (1.0 - alpha) / a0_;
        // 不清零状态，保证系数渐变时平滑
    }

    // 中位数
    double median_period() const
    {
        std::vector<double> v(Ts.begin(), Ts.end());
        std::nth_element(v.begin(), v.begin() + v.size() / 2, v.end());
        return v[v.size() / 2];
    }

    // 截尾平均
    double robust_T() const
    {
        if (Ts.size() <= 3)
            return median_period();
        std::vector<double> v(Ts.begin(), Ts.end());
        std::sort(v.begin(), v.end());
        int n = (int)v.size();
        int k = (int)std::floor(n * trim_ratio);
        k = std::min(k, n / 3);
        double sum = 0.0;
        int cnt = 0;
        for (int i = k; i < n - k; ++i) {
            sum += v[i];
            ++cnt;
        }
        return (cnt > 0) ? (sum / cnt) : v[n / 2];
    }
};

/** FPGA 与主程序之间的数据结构 **/
struct receive_data {
    double Ia { 0.0 }; // 注入电网电流A (kA)
    double Ib { 0.0 }; // 注入电网电流B (kA)
    double Ic { 0.0 }; // 注入电网电流C (kA)
    double Vapp { 0.0 };
    double Vbpp { 0.0 };
    double Vcpp { 0.0 };
    double Vdc1 { 0.0 }; // 直流负极电压 (kV)
    double Vdc2 { 0.0 }; // 直流正极电压 (kV)
};

struct send_data {
    double Va { 0.0 }; // 并网点A相电压(kV)
    double Vb { 0.0 };
    double Vc { 0.0 };
    double Idc { 0.0 }; // 直流侧等效电流源(kA) = P_cmd/Vdc
    double Vdoor_va { 0 }; // 归一化调制信号(-1~1)
    double Vdoor_vb { 0 };
    double Vdoor_vc { 0 };
};

const int T_type_R = 1; // 电阻
const int T_type_RL = 2; // 梯形离散的串联RL（恒等效导纳+历史电流）
const int T_type_I = 3; // IGBT 支路（开关恒导纳+历史项）
const int T_type_C = 4; // 电容（恒等效导纳+历史电流）

struct T_para_ {
    int type { 0 };
    string name;
    int node_I { 0 };
    int node_J { 0 };
    double value { 0.0 };
    double hisI { 0.0 };
    double hisV { 0.0 };
    double hist { 0.0 };
    int old_c { 0 }; // IGBT前次通断状态
    int next_IGBT { 0 };
};

/**
 * @brief 采用 LC 恒定导纳模型并模拟 1us 小步长的“FPGA”计算
 *  - 内部电路不大，直接构建稠密导纳并求逆，每小步只做 Z*I
 *  - 三角波比较在此处完成
 */
class FPGA_simu {
public:
    void buid_model();
    void cal_step(send_data& send_, receive_data& receive_);
    friend class Converter;

private:
private:
    // ===== 三角波 LUT =====
    std::vector<double> tri_LUT_; // 一个开关周期内的 y0 样值（幅值[-1,1]）
    double tri_step_ { 1e-6 }; // 小步长 (s)，与 cal_step 内一致
    int tri_len_ { 0 }; // 一个周期内的采样点个数

    // ==== PWM/时间 ====
    double fre { 5000.0 }; // 开关频率 (Hz)
    double simu_time_ { 0 }; // 当前仿真时间 (s)

    // ==== LC 模型参数 ====
    const double L_ { 0.016 }; // H
    const double R_ { 0.01 }; // Ω
    const double C_uF_ { 5000 }; // μF
    const double R_AC_ { 0.05 }; // Ω  交流侧等效内阻
    const double Rs_DC_ { 0.001 }; // Ω  直流侧等效内阻
    const double alpha_ { 2.41421e-3 };
    const double belta_ { 0.414214e-3 };

    // ==== 拓扑 ====
    int N_ = 8;
    int line_num_ = 16;

    vector<T_para_> T_vec;
    Eigen::MatrixXd Y_mat;
    Eigen::MatrixXd Y_mat_inv;

    double value1_RL { 0.0 };
    double value2_RL { 0.0 };
    double C_value { 0.0 };
    double RL_value { 0.0 };

    Eigen::VectorXd V; // 电压向量（kV）
    Eigen::VectorXd I; // 电流向量（kA）

    const double Vdc_target_kV_ { 1.0 }; // 直流电压(我们假设直流侧有大电容，基本恒定)
    const double R_hold_ { 0.1 };
};

struct VSGParams {
    // 以 pu/Hz、pu/pu 的下垂系数；额定值由 Converter 内部 S_nom、V0/f0 决定
    double H_s = 0.5; // 惯量常数 H (s)
    double D_pu = 2.0; // 阻尼系数（pu）
    double Kf_pu_per_Hz = 0.2; // 有功-频率下垂（pu/Hz）
    double Kq_pu_per_pu = 0.1; // 无功-电压下垂（pu/pu）
    double P0_MW = 0.0; // 目标有功（MW）默认用 Pref_MW_
    double Q0_MVAr = 0.0; // 目标无功（MVAr）默认用 Qref_MVAr_
    double V0_LL_kV = 0.38; // 额定线电压(kV)
    double Rv = 0.01; // 虚拟电阻(Ω)
    double Lv = 0.0; // 虚拟电感(H)
    double Vmin_pu = 0.9, Vmax_pu = 1.1;
    double fmin_Hz = 45.0, fmax_Hz = 55.0;
}; // 虚拟同步机参数

/**
 * @brief 三相并网友侧变流器（PLL + 外环PQ + 内环电流PI + SPWM）
 *  - 对外体现为“注入型”诺顿电流源：在 updateHistory 中将 Ia/Ib/Ic 注入 I 向量
 *  - 单位体系：电压[kV]，电流[kA]，阻抗[Ω]，功率[MW/MVAr]
 */
class Converter : public Device {
public:
    // 三相并网：指定 A/B/C 三个节点；默认家庭光伏 ~10kW，Q=0
    Converter(int na, int nb, int nc, double P_ref_MW = 0.01, double Q_ref_MVAr = 0.0);

    // 运行前设置（初始化阶段选择，不在运行中切换）
    void set_control_mode(ControlMode m)
    {
        mode_ = m;
        if (mode_ == ControlMode::VSG) {
            if (vsg_.P0_MW == 0.0)
                vsg_.P0_MW = Pref_MW_;
            if (vsg_.Q0_MVAr == 0.0)
                vsg_.Q0_MVAr = Qref_MVAr_;
            // 额定相电峰值 = V_LL * sqrt(2/3)
            V_phase_base_kV_ = vsg_.V0_LL_kV * std::sqrt(2.0 / 3.0);
            // 初始化角速度为额定
            omega_ = omega_nom_;
            theta_ = 0.0;
        }
    }

    // 便于录波：角度与频率（Hz）
    double get_theta() const { return theta_; }
    double get_omega_Hz() const { return omega_ / (2.0 * M_PI); }

    // Device 接口
    void stamp(std::vector<Eigen::Triplet<double>>& triplets, double dt) override { /* 注流型，无导纳印章 */ }
    void updateHistory(Eigen::VectorXd& I, double t, double dt) override;
    void updateState(const Eigen::VectorXd& V, double dt) override; /// @brief 采样当前并网点三相电压（kV），供下一拍控制与注流计算使用
    void allocateNodes(Grid& grid) override { } // 无内部节点

    // 监测
    double get_Ia() const { return iA_kA_; }
    double get_Ib() const { return iB_kA_; }
    double get_Ic() const { return iC_kA_; }

    // 功率读数（MW / MVAr）
    double get_P_MW() const { return P_f_; }
    double get_Q_MVAr() const { return Q_f_; }

    // 下垂
    void enable_droop(bool en) { droop_enabled_ = en; }
    void set_droop_coeff(double Kf_MW_per_Hz, double Kv_MVAr_per_kV)
    {
        Kf_ = Kf_MW_per_Hz;
        Kv_ = Kv_MVAr_per_kV;
    }

private:
    // 额定参数
    double S_nom { 0.02 }; // 额定容量(MVA)
    double V_LL_kV_const_ { 0.38 }; // 额定线电压(kV)

    // === 模式选择 ===
    ControlMode mode_ { ControlMode::PLL };
    VSGParams vsg_;
    bool vsg_theta_inited = false;
    double V_phase_base_kV_ { V_LL_kV_const_ * std::sqrt(2.0 / 3.0) }; // 缺省按 0.38kV 初始化

    // === 将原 updateHistory/State 切成两套实现 ===
    void pll_updateHistory(Eigen::VectorXd& I, double t, double dt);
    void vsg_updateHistory(Eigen::VectorXd& I, double t, double dt);

    // ====== 连接 ======
    int nA_, nB_, nC_;

    // ====== 并网测量缓存（kV/kA） ======
    double vA_kV_ { 0 }, vB_kV_ { 0 }, vC_kV_ { 0 };
    double iA_kA_ { 0 }, iB_kA_ { 0 }, iC_kA_ { 0 };

    ZeroCrossFreqWin freq_e2; // 过零测频（窗口法）

    // ====== PLL（SRF-PLL） ======
    double theta_ { 0.0 };
    double omega_ { 2.0 * PI * 50.0 }; // rad/s
    const double omega_nom_ { 2.0 * PI * 50.0 };
    const double Kp_pll_ { 200 };
    const double Ki_pll_ { 40000 };
    double pll_int_ { 0.0 };

    // ====== 外环 PQ ======
    double Pref_MW_, Qref_MVAr_;
    double P_f_ { 0.0 }, Q_f_ { 0.0 };
    const double Tpq_ { 0.10 }; // s
    const double KpP_base_ { 0.6 };
    const double KiP_base_ { 94.2 }; // 2π·15
    const double KpQ_base_ { 0.6 };
    const double KiQ_base_ { 94.2 };
    double eP_int_ { 0.0 }, eQ_int_ { 0.0 };

    const double i_limit_factor_ { 1.5 };

    // ====== 内环电流 PI（带解耦 + 电压前馈） ======
    const double w_ci_ { 2.0 * M_PI * 300.0 }; // 目标带宽 300 Hz
    const double R_eq_for_design_ { 0.02 }; // 设计用等效电阻(Ω)，比原来0.1小很多
    const double Kp_i_ { w_ci_ * 0.016 }; // = 1885 * 0.016 ≈ 30.16（保持不变）
    const double Ki_i_ { w_ci_ * R_eq_for_design_ }; // = 1885 * 0.02 ≈ 37.7（原来是 188.5）

    double R_virtual_damp_ { 0.04 };

    double id_int_ { 0.0 }, iq_int_ { 0.0 };

    // ====== 直流侧 ======
    double Vdc_kV_; // 估计直流母线 (kV)

    // ====== 通信到“FPGA” ======
    send_data send_;
    receive_data receive_;
    FPGA_simu fpga_;

    // ====== 下垂控制 ======
    bool droop_enabled_ { true }; // 默认关闭以实现 Q≈0（需要时可 enable）
    // 注意两者通常为负值（f低/电压低时增加出力）
    double Kf_ { -0.004 }; // MW/Hz
    double Kv_ { -0.5 }; // MVAr/kV

    // 功率读数（即时量）
    double P_inst_MW_ { 0.0 }, Q_inst_MVAr_ { 0.0 };

    double U_mag_f_ { 0.0 }; // 每台机自己的 |ud| 平滑值

    // ====== 新增：用于滤波的成员变量 ======
    const double T_filter_ { 0.005 }; // 滤波器时间常数 (s)
    double ud_f_ { 0.0 }, uq_f_ { 0.0 }; // 滤波后的电压
    double id_f_ { 0.0 }, iq_f_ { 0.0 }; // 滤波后的电流

    bool saturated_ { false };
    int sat_hits_ { 0 };
};
