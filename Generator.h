/* =============================================================================
 * Generator.h - 同步发电机模型
 * -----------------------------------------------------------------------------
 * 支持PV节点和平衡节点的同步发电机实现，用于IEEE 9节点系统等电力系统电磁暂态仿真
 * 
 * PV节点：已知有功功率P和电压幅值|V|，计算无功功率Q和电压相角θ
 * 平衡节点：已知电压幅值|V|和相角θ，计算有功功率P和无功功率Q
 * ============================================================================= */
 
#pragma once

#include "Devices.h"
#include <Eigen/Dense>
#include <cmath>
#include <string>
#include <vector>

using std::string;
using std::vector;
/**
 * @brief 同步发电机模型，支持PV节点和平衡节点。
 * 
 * PV节点：已知有功功率P和电压幅值|V|，计算无功功率Q和电压相角θ
 * 平衡节点：已知电压幅值|V|和相角θ，计算有功功率P和无功功率Q
 */
class Generator : public VoltageSource {
public:
    /// @brief 发电机节点类型
    enum class NodeType {
        PV,         ///< PV节点：已知P和|V|
        SLACK       ///< 平衡节点：已知|V|和θ
    };
 
    /**
     * @brief 构造一个同步发电机。
     * @param node 连接节点编号。
     * @param type 节点类型（PV或SLACK）。
     * @param P_MW 有功功率 (MW)，PV节点为设定值，SLACK节点为初值。
     * @param V_kV 电压幅值 (kV)，两种节点都是设定值。
     * @param theta_deg 电压相角 (度)，SLACK节点为设定值，PV节点为初值。
     * @param Xd 直轴同步电抗 (pu)。
     * @param Xq 交轴同步电抗 (pu)。
     * @param Ra 电枢电阻 (pu)。
     * @param V_base_kV 基准电压 (kV)，用于标幺化。
     * @param S_base_MVA 基准功率 (MVA)，用于标幺化。
     */
    Generator(int node, NodeType type,
        double P_MW, double V_kV, double theta_deg,
        double Xd, double Xq, double Ra,
        double V_base_kV = 110.0, double S_base_MVA = 100.0);
 
    void stamp(std::vector<Eigen::Triplet<double>>& triplets, double dt) override;
    void updateHistory(Eigen::VectorXd& I, double t, double dt) override;
    void updateState(const Eigen::VectorXd& V, double dt) override;
 
    // --- 状态查询 ---
    double getP_MW() const { return P_calc_MW_; }
    double getQ_MVAr() const { return Q_calc_MVAr_; }
    double getVoltageAngle_deg() const { return theta_calc_deg_; }
    double getInternalVoltage_pu() const { return E_internal_pu_; }
    double getPowerAngle_deg() const { return delta_deg_; }
 
    // --- 约束设置（用于运行中调整） ---
    void setPVConstraints(double P_MW, double V_kV);
    void setSlackConstraints(double V_kV, double theta_deg = 0.0);
 
    double get_I() { return I_inject_kA_ * 1000.0; } ///< 返回注入电流 (A)
 
private:
    // 基本参数
    int node_;                  ///< 连接节点
    NodeType type_;             ///< 节点类型
    double V_base_kV_;          ///< 基准电压 (kV)
    double S_base_MVA_;         ///< 基准功率 (MVA)
    double Z_base_ohm_;         ///< 基准阻抗 (Ω)
 
    // 发电机参数
    double Xd_pu_, Xq_pu_, Ra_pu_;  ///< 同步机参数 (pu)
 
    // 约束条件
    double P_set_MW_;           ///< 设定有功功率 (MW)
    double V_set_kV_;           ///< 设定电压幅值 (kV)
    double theta_set_deg_;      ///< 设定电压相角 (度)，仅SLACK节点使用
 
    // 计算结果
    double P_calc_MW_;          ///< 计算的有功功率 (MW)
    double Q_calc_MVAr_;        ///< 计算的无功功率 (MVAr)
    double theta_calc_deg_;     ///< 计算的电压相角 (度)
    double E_internal_pu_;      ///< 内部电势幅值 (pu)
    double delta_deg_;          ///< 功角 (度)
 
    // 电流注入
    double I_inject_kA_;        ///< 注入电流幅值 (kA)
 
    // 收敛控制
    static constexpr double CONV_TOL_ = 1e-6;   ///< 收敛容差
    static constexpr int MAX_ITER_ = 20;        ///< 最大迭代次数
 
    // 内部计算方法
    void solvePVNode(double V_terminal_kV, double dt);
    void solveSlackNode(double dt);
    void calculateInternalVoltage();
    void calculateInjectionCurrent();
};