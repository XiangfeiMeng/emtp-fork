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
#include "Grid.h"
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
class Generator : public Device {
public:
    enum class NodeType {
        PV,         // 定P、定|V|，计算Q、θ
        SLACK       // 定|V|、定θ，计算P、Q
    };

    Generator(int n_pos, int n_neu, NodeType type,
        double P_MW, double V_kV, double phase_deg,
        double V_base_kV, double S_base_MVA);

    void stamp(std::vector<Eigen::Triplet<double>>& triplets, double dt) override;
    void updateHistory(Eigen::VectorXd& I, double t, double dt) override;
    void updateState(const Eigen::VectorXd& V, double dt) override;

    // 状态查询（简化）
    double getP_MW() const { return P_calc_MW_; }
    double getQ_MVAr() const { return Q_calc_MVAr_; }
    double getVoltageAngle_deg() const { return phase_calc_deg_; }
    double get_I() const { return I_inject_kA_ * 1000.0; } // 返回电流(A)

    // 约束调整接口
    void setPVConstraints(double P_MW, double V_kV) {
        if (type_ == NodeType::PV) {
            P_set_MW_ = P_MW;
            V_set_kV_ = V_kV;
        }
    }
    void setSlackConstraints(double V_kV, double phase_deg) {
        if (type_ == NodeType::SLACK) {
            V_set_kV_ = V_kV;
            phase_set_deg_ = phase_deg;
        }
    }

private:
    int n_pos_, n_neu_;                  // 连接节点
    NodeType type_;             // 节点类型

    // 基准参数
    double V_base_kV_;          // 基准电压(kV)
    double S_base_MVA_;         // 基准功率(MVA)

    // 设定值
    double P_set_MW_;           // 有功设定值(MW)
    double V_set_kV_;           // 电压幅值设定值(kV)
    double phase_set_deg_;      // 相角设定值(度)

    // 计算值
    double P_calc_MW_;          // 计算有功(MW)
    double Q_calc_MVAr_;        // 计算无功(MVAr)
    double phase_calc_deg_;     // 计算相角(度)
    double I_inject_kA_;        // 注入电流(kA)

     // 新增：对Grid的引用（用于平衡节点功率计算）
    const Grid* grid_ptr_ = nullptr;  
    
    // 修改后的函数声明
    void solveSlackNode(const Eigen::VectorXd& V);  // 需要传入电压向量

    // 内部计算
    void solvePVNode(double V_terminal_kV);   // 简化PV节点求解
    std::complex<double> calculateInjectionCurrent(const Eigen::VectorXd& V) ;
    void calculateSlackPowerFromSystemState(const Eigen::VectorXd& V) ;
};