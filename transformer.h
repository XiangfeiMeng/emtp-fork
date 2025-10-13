/* =============================================================================
 * 三相双绕组变压器（Y/Δ，±30° 相移，可设中性点接地）
 * -----------------------------------------------------------------------------
 * - 漏抗：由短路试验参数折到 #1 侧，再按两侧串联各分一半（series_RL）
 * - 励磁：两段法（线性段/空心段），以磁链 ψ 为判据；段切换仅触发“数值重分解”
 * - 约束：理想变比用拉格朗日变量（额外结点）实现，不改变外部拓扑
 * - 单位：V[kV] / I[kA]；名牌参数按 MVA/kV/%/kW 输入（见 Transformer_para）
 * ============================================================================= */

#pragma once

#include "Devices.h"
#include "Grid.h"
#include <Eigen/Dense>
#include <algorithm>
#include <array>
#include <atomic>
#include <cmath>
#include <memory>
#include <string>
#include <vector>

using std::string;
using std::vector;

struct Transformer_para {
    // —— 名牌/额定参数 ——
    double load { 50.0 }; // 额定容量 MVA
    double freq { 50.0 }; // 额定频率 Hz
    double swinging1Voltage { 110.0 }; // #1侧额定线电压 kV
    double swinging2Voltage { 35.0 }; // #2侧额定线电压 kV

    // —— 饱和与励磁（两段法） ——
    int saturationEnabled { 1 }; // 0=线性；1=两段法
    double airCoreReactance { 0.2 }; // pu，相对线性励磁电抗
    double decayConstant { 1.0 }; // s（保留接口）
    double kneeVoltage { 1.25 }; // pu（相电压拐点）

    // —— 空载/短路试验近似（总值） ——
    double noLoadLosses { 50.0 }; // kW
    double noLoadpercentage { 0.5 }; // %（额定电流百分比）
    double short_circuit_loss { 300.0 }; // kW
    double short_circuit_percentage { 12.0 }; // %（短路电压）

    // —— 接线与相移 ——
    // 0=Y, 1=Δ
    int swinging1type { 0 }; // #1侧接线（默认Y）
    int swinging2type { 1 }; // #2侧接线（默认Δ）
    int swinging1_GROUNDING { 1 }; // #1侧中性点接地与否（Y侧有效）
    int swinging2_GROUNDING { 0 }; // #2侧中性点（Δ侧忽略）
    int angle_delta { +30 }; // {0,+30,-30}

    // —— 分接头（#1侧） ——
    int tapMid { 0 };
    double ratioTap { 1.25 }; // 每档电压变化率 %
    int tap_now { 0 };
    int tap_max { +9 };
    int tap_min { -9 };
};

// ================== 内部：两段线性励磁 ==================
// 非线性励磁支路（两段法）：TRAP 离散为诺顿等效；以 ψ=∫v dt 作段选
// 段切换通过原子标志通知 Simulation 做“仅数值重分解”，避免频繁结构重建。
class NLmag : public Device {
public:
    NLmag(int n1, int n2,
        double L_lin, double L_air, double psi_knee,
        std::atomic<bool>* refFlag)
        : n1_(n1)
        , n2_(n2)
        , L_lin_(std::max(1e-6, L_lin))
        , L_air_(std::max(1e-6, L_air))
        , psi_knee_(std::max(1e-6, psi_knee))
        , pNeedRefactor_(refFlag)
    {
    }

    void stamp(std::vector<Eigen::Triplet<double>>& triplets, double dt) override;
    void updateHistory(Eigen::VectorXd& I, double t, double dt) override;
    void updateState(const Eigen::VectorXd& V, double dt) override;

    int segment() const { return seg_; }

    // 便捷接口：设残磁（以拐点磁链为 1.0 pu）、观测
    void setResidualFluxPU(double pu)
    {
        psi_ = pu * psi_knee_;
        pickSegmentFromPsi();
    }
    bool isSaturated() const;
    double i_mag_now() const;

private:
    double t_now_;
    int n1_ { 0 }, n2_ { 0 };
    // 段0：线性段；段1：空心段
    int seg_ { 0 };
    double L_lin_, L_air_;
    double psi_knee_;
    double L_eff_ { 0.0 }, G_eq_ { 0.0 };

    // 历史量
    double i_hist_ { 0.0 };
    double v_hist_ { 0.0 };
    double psi_ { 0.0 };
    double i_now_ { 0.0 };

    std::atomic<bool>* pNeedRefactor_ { nullptr };

    inline void pickSegmentFromPsi()
    {
        seg_ = (std::abs(psi_) <= psi_knee_) ? 0 : 1;
        L_eff_ = (seg_ == 0) ? L_lin_ : L_air_;
    }
};

// =============== 内部：单相理想变压器相支路容器 ===============
class SP_Unit {
public:
    // int k { 0 }, m { 0 }, j { 0 }, l { 0 }, extra { 0 }; // #1端口(k,m)，#2端口(j,l)，extra=约束变量
    //  外部端子（与系统相连）
    int k { 0 }, m { 0 }, j { 0 }, l { 0 };
    // 绕组内部端子（实现串联泄漏）
    int kp { 0 }, mp { 0 }, js { 0 }, ls { 0 };
    // 理想变压器约束的拉格朗日变量（用作额外未知量）
    int extra { 0 };
    std::unique_ptr<class series_RL> Zp; // #1侧串联漏抗(一半)
    std::unique_ptr<class series_RL> Zs; // #2侧串联漏抗(一半，折算)
    std::unique_ptr<NLmag> mag; // 励磁（接在#1侧）
    double a_eff { 1.0 }; // 本相有效变比
    std::string name;
};

// ================== 双绕组三相变压器 ==================
class Transformer : public Device {
public:
    Transformer(std::vector<int> nodes_i, std::vector<int> nodes_j, Transformer_para& para_);

    void allocateNodes(Grid& grid) override; /// @brief 申请中性点/额外约束节点，构建相支路并创建漏抗与励磁元件
    void stamp(std::vector<Eigen::Triplet<double>>& triplets, double dt) override;
    void updateHistory(Eigen::VectorXd& I, double t, double dt) override;
    void updateState(const Eigen::VectorXd& V, double dt) override;

    // —— 非线性段切换全局标志（Simulation 轮询） ——
    static bool consumeRefactorFlag() { return s_needRefactor.exchange(false); }
    void setResidualFluxPU(double a, double b, double c)
    {
        if (ph_[0].mag)
            ph_[0].mag->setResidualFluxPU(a);
        if (ph_[1].mag)
            ph_[1].mag->setResidualFluxPU(b);
        if (ph_[2].mag)
            ph_[2].mag->setResidualFluxPU(c);
    }

private:
    std::vector<int> nodes_i_; // #1侧 A,B,C
    std::vector<int> nodes_j_; // #2侧 A,B,C
    Transformer_para para_;

    double omega_ { 2.0 * 3.141592653589793 * 50.0 };

    // #1/#2 相电压（按接线）
    double V1_ph_ { 0.0 }, V2_ph_ { 0.0 };
    double a_eff_ { 1.0 };

    // 中性点（若 Y）
    int n1_neu_ { 0 };
    int n2_neu_ { 0 };

    std::array<SP_Unit, 3> ph_;

    // 由名牌得到（折到#1侧）
    double R_eq_ { 0.0 }, X_eq_ { 0.0 };
    double Rm_ { 0.0 }, Lm_lin_ { 0.0 };
    double Lm_air_ { 0.0 };
    double psi_knee_ { 0.0 };

    void computeParameters();
    void buildPhaseConnection(Grid& grid);
    void stampIdealConstraint(std::vector<Eigen::Triplet<double>>& triplets, const SP_Unit& u) const;

    // 非线性切段 -> 通知 Simulation 需要“仅数值重分解”
    static std::atomic<bool> s_needRefactor;
};
