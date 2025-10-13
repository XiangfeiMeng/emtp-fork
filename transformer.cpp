#include "transformer.h"
#include <cmath>
#include <iostream>

std::atomic<bool> Transformer::s_needRefactor { false };

// ================== NLmag ==================
void NLmag::stamp(std::vector<Eigen::Triplet<double>>& triplets, double dt)
{
    // 使用上一轮判定的段（updateState 已更新 seg_）
    L_eff_ = (seg_ == 0) ? L_lin_ : L_air_;
    double G = dt / (2.0 * std::max(1e-9, L_eff_));
    if (n1_ > 0)
        triplets.emplace_back(n1_ - 1, n1_ - 1, G);
    if (n2_ > 0)
        triplets.emplace_back(n2_ - 1, n2_ - 1, G);
    if (n1_ > 0 && n2_ > 0) {
        triplets.emplace_back(n1_ - 1, n2_ - 1, -G);
        triplets.emplace_back(n2_ - 1, n1_ - 1, -G);
    }
    G_eq_ = G;
}

void NLmag::updateHistory(Eigen::VectorXd& I, double t, double dt)
{
    t_now_ = t;
    i_hist_ = i_hist_ + 2.0 * v_hist_ * G_eq_;
    if (n1_ > 0)
        I(n1_ - 1) -= i_hist_;
    if (n2_ > 0)
        I(n2_ - 1) += i_hist_;
}

void NLmag::updateState(const Eigen::VectorXd& V, double dt)
{
    double v1 = (n1_ > 0) ? V(n1_ - 1) : 0.0;
    double v2 = (n2_ > 0) ? V(n2_ - 1) : 0.0;
    double v = v1 - v2;

    // 当前电流（诺顿）
    double i_now = i_hist_ + G_eq_ * v;
    i_now_ = i_now;

    // 磁链积分
    // psi_ += v * dt;
    psi_ += 0.5 * dt * (v + v_hist_);

    // 判定段
    int seg_new = (std::abs(psi_) <= psi_knee_) ? 0 : 1;
    if (seg_new != seg_) {
        seg_ = seg_new;
        if (pNeedRefactor_)
            pNeedRefactor_->store(true); // 通知需要“数值重分解”
        // std::cout << " t =  " << t_now_ << " ; [NLmag] segment -> " << (seg_ ? "AIR" : "LIN")
        //<< "  psi=" << psi_ << "  knee=" << psi_knee_ << "\n";
    }

    v_hist_ = v;
}

// ================== Transformer ==================
Transformer::Transformer(std::vector<int> nodes_i, std::vector<int> nodes_j, Transformer_para& para)
    : nodes_i_(std::move(nodes_i))
    , nodes_j_(std::move(nodes_j))
    , para_(para)
{
    omega_ = 2.0 * 3.141592653589793 * std::max(1e-3, para_.freq);
}

void Transformer::computeParameters()
{
    auto phase_voltage = [](double Vll, int typeYorDelta) -> double {
        return (typeYorDelta == 0) ? (Vll / std::sqrt(3.0)) : Vll;
    };
    V1_ph_ = phase_voltage(para_.swinging1Voltage, para_.swinging1type);
    V2_ph_ = phase_voltage(para_.swinging2Voltage, para_.swinging2type);

    // 分接头装在#1侧
    double tap_scale = 1.0 + (para_.ratioTap * 0.01) * (para_.tap_now - para_.tapMid);
    a_eff_ = (V1_ph_ / std::max(1e-6, V2_ph_)) * tap_scale; // n1/n2

    // 额定线电流（#1侧）
    double I1_line = para_.load / std::max(1e-6, std::sqrt(3.0) * para_.swinging1Voltage);
    double P0_MW = para_.noLoadLosses / 1000.0;
    double Psc_MW = para_.short_circuit_loss / 1000.0;

    // 短路试验（折到#1侧）
    double Z_eq = (para_.swinging1Voltage * para_.swinging1Voltage)
        / std::max(1e-6, para_.load)
        * (para_.short_circuit_percentage * 0.01);
    double R_eq = Psc_MW / std::max(1e-9, 3.0 * I1_line * I1_line);
    R_eq_ = std::max(1e-6, R_eq);
    X_eq_ = std::sqrt(std::max(1e-12, Z_eq * Z_eq - R_eq_ * R_eq_));

    // 空载试验（#1侧相量）
    double I0 = para_.noLoadpercentage * 0.01 * I1_line;
    double Iw = P0_MW / std::max(1e-6, 3.0 * V1_ph_);
    double Im = std::sqrt(std::max(0.0, I0 * I0 - Iw * Iw));
    Rm_ = (Iw > 1e-9) ? (V1_ph_ / Iw) : 1e12; // 铁损电阻（并联）
    double Xm = (Im > 1e-9) ? (V1_ph_ / Im) : 1e12;
    Lm_lin_ = Xm / omega_;

    // 空心段
    // double Xm_air = std::max(1e-6, para_.airCoreReactance) * Xm;
    // Lm_air_ = Xm_air / omega_;
    double Z_base = (para_.swinging1Voltage * para_.swinging1Voltage)
        / std::max(1e-6, para_.load); // Ω
    double X_air = std::max(1e-6, para_.airCoreReactance) * Z_base; // pu(基值)→物理值
    Lm_air_ = X_air / omega_;

    double Vknee = std::max(0.1, para_.kneeVoltage) * V1_ph_;
    psi_knee_ = Vknee / omega_;
}

void Transformer::buildPhaseConnection(Grid& grid)
{
    // —— 中性点策略：若接地则直接用 0（地）；不接地才新建节点 ——
    if (para_.swinging1type == 0) {
        n1_neu_ = para_.swinging1_GROUNDING ? 0 : grid.addNode();
    }
    if (para_.swinging2type == 0) {
        n2_neu_ = para_.swinging2_GROUNDING ? 0 : grid.addNode();
    }

    auto pair_for_delta = [&](int phase, int sign) -> std::pair<int, int> {
        const int A = 0, B = 1, C = 2;
        if (sign >= 0) { // +30
            if (phase == A)
                return { nodes_j_[A], nodes_j_[B] };
            if (phase == B)
                return { nodes_j_[B], nodes_j_[C] };
            return { nodes_j_[C], nodes_j_[A] };
        } else { // -30
            if (phase == A)
                return { nodes_j_[B], nodes_j_[A] };
            if (phase == B)
                return { nodes_j_[C], nodes_j_[B] };
            return { nodes_j_[A], nodes_j_[C] };
        }
    };
    int sign = (para_.angle_delta >= 0) ? +1 : -1;

    for (int p = 0; p < 3; ++p) {
        auto& u = ph_[p];
        u.name = std::string("ph") + char('A' + p);
        u.a_eff = a_eff_;
        u.extra = grid.addNode();

        // #1侧端口
        if (para_.swinging1type == 0) { // Y
            u.k = nodes_i_[p];
            u.m = n1_neu_;
        } else { // Δ
            if (p == 0) {
                u.k = nodes_i_[0];
                u.m = nodes_i_[1];
            }
            if (p == 1) {
                u.k = nodes_i_[1];
                u.m = nodes_i_[2];
            }
            if (p == 2) {
                u.k = nodes_i_[2];
                u.m = nodes_i_[0];
            }
        }

        // #2侧端口
        if (para_.swinging2type == 0) {
            u.j = nodes_j_[p];
            u.l = n2_neu_;
        } else {
            auto pr = pair_for_delta(p, sign);
            u.j = pr.first;
            u.l = pr.second;
        }
    }
}

void Transformer::allocateNodes(Grid& grid)
{
    computeParameters();
    buildPhaseConnection(grid);

    double Rp_1 = 0.5 * R_eq_;
    double Xp_1 = 0.5 * X_eq_;
    double Lp_1 = Xp_1 / omega_;
    double Rs_2 = (Rp_1) / (a_eff_ * a_eff_);
    double Ls_2 = (Lp_1) / (a_eff_ * a_eff_);

    for (int p = 0; p < 3; ++p) {
        auto& u = ph_[p];
        u.Zp = std::make_unique<series_RL>(u.k, u.m, std::max(1e-6, Rp_1), std::max(1e-6, Lp_1));
        u.Zs = std::make_unique<series_RL>(u.j, u.l, std::max(1e-6, Rs_2), std::max(1e-6, Ls_2));
        u.mag = std::make_unique<NLmag>(u.k, u.m, Lm_lin_, Lm_air_, psi_knee_, &s_needRefactor);
    }
}

static inline void add_triplet(std::vector<Eigen::Triplet<double>>& T, int r, int c, double v)
{
    if (r > 0 && c > 0)
        T.emplace_back(r - 1, c - 1, v);
}

void Transformer::stampIdealConstraint(std::vector<Eigen::Triplet<double>>& triplets, const SP_Unit& u) const
{
    // 修正：理想关系为 (v_k - v_m) - a*(v_j - v_l) = 0
    add_triplet(triplets, u.k, u.extra, +1.0);
    add_triplet(triplets, u.m, u.extra, -1.0);
    add_triplet(triplets, u.j, u.extra, -u.a_eff);
    add_triplet(triplets, u.l, u.extra, +u.a_eff);

    add_triplet(triplets, u.extra, u.k, +1.0);
    add_triplet(triplets, u.extra, u.m, -1.0);
    add_triplet(triplets, u.extra, u.j, -u.a_eff);
    add_triplet(triplets, u.extra, u.l, +u.a_eff);

    // 为防止 0 主元，给 extra 加一个极小 gmin
    const double gmin_extra = 1e-6;
    add_triplet(triplets, u.extra, u.extra, gmin_extra);
}

void Transformer::stamp(std::vector<Eigen::Triplet<double>>& triplets, double dt)
{
    for (auto& u : ph_) {
        // 理想变比约束
        stampIdealConstraint(triplets, u);

        // 漏抗（两侧各一半）
        if (u.Zp)
            u.Zp->stamp(triplets, dt);
        if (u.Zs)
            u.Zs->stamp(triplets, dt);

        // 励磁非线性支路（并联在 k–m）
        if (u.mag)
            u.mag->stamp(triplets, dt);

        // —— 铁损电阻并联（不改头文件，直接以电导形式 stamp） ——
        if (Rm_ > 0.0) {
            double Gm = 1.0 / Rm_;
            add_triplet(triplets, u.k, u.k, +Gm);
            add_triplet(triplets, u.m, u.m, +Gm);
            add_triplet(triplets, u.k, u.m, -Gm);
            add_triplet(triplets, u.m, u.k, -Gm);
        }
    }

    // 若#1/#2侧是Y但不接地，则注入极小对地导纳避免悬浮
    if (para_.swinging1type == 0 && !para_.swinging1_GROUNDING && n1_neu_ > 0)
        triplets.emplace_back(n1_neu_ - 1, n1_neu_ - 1, 1.0 / 1e9);
    if (para_.swinging2type == 0 && !para_.swinging2_GROUNDING && n2_neu_ > 0)
        triplets.emplace_back(n2_neu_ - 1, n2_neu_ - 1, 1.0 / 1e9);
}

void Transformer::updateHistory(Eigen::VectorXd& I, double t, double dt)
{
    for (auto& u : ph_) {
        if (u.Zp)
            u.Zp->updateHistory(I, t, dt);
        if (u.Zs)
            u.Zs->updateHistory(I, t, dt);
        if (u.mag)
            u.mag->updateHistory(I, t, dt);
    }
}

void Transformer::updateState(const Eigen::VectorXd& V, double dt)
{
    for (auto& u : ph_) {
        if (u.Zp)
            u.Zp->updateState(V, dt);
        if (u.Zs)
            u.Zs->updateState(V, dt);
        if (u.mag)
            u.mag->updateState(V, dt); // 若段切换，将置 s_needRefactor=true
    }
}
