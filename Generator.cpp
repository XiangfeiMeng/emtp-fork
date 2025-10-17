/* =============================================================================
 * Generator.cpp - 同步发电机模型实现
 * -----------------------------------------------------------------------------
 * PV节点和平衡节点发电机的电磁暂态仿真实现
 * ============================================================================= */
 
#include "Generator.h"
#include <cmath>
#include <iostream>
// --- Generator ---
Generator::Generator(int node, NodeType type,
    double P_MW, double V_kV, double theta_deg,
    double Xd, double Xq, double Ra,
    double V_base_kV, double S_base_MVA)
    : node_(node)
    , type_(type)
    , V_base_kV_(V_base_kV)
    , S_base_MVA_(S_base_MVA)
    , Xd_pu_(Xd)
    , Xq_pu_(Xq)
    , Ra_pu_(Ra)
    , P_set_MW_(P_MW)    , V_set_kV_(V_kV)
    , theta_set_deg_(theta_deg)
    , P_calc_MW_(P_MW)
    , Q_calc_MVAr_(0.0)
    , theta_calc_deg_(theta_deg)
    , E_internal_pu_(1.0)
    , delta_deg_(0.0)
    , I_inject_kA_(0.0)
{
    Z_base_ohm_ = V_base_kV_ * V_base_kV_ / S_base_MVA_;
    
    // 初始化内部电势
    calculateInternalVoltage();
}
 
void Generator::stamp(std::vector<Eigen::Triplet<double>>& triplets, double dt)
{
    // 发电机作为电流源，不贡献导纳矩阵
    // 仅在需要阻尼时添加小电导（类似VoltageSource的处理）
    if (node_ > 0) {
        const double g_damp = 1e-12; // 数值阻尼
        triplets.emplace_back(node_ - 1, node_ - 1, g_damp);
    }
}
 
void Generator::updateHistory(Eigen::VectorXd& I, double t, double dt)
{
    if (node_ <= 0) return;
 
    // 注入计算得到的电流
    I(node_ - 1) += I_inject_kA_;
}
 
void Generator::updateState(const Eigen::VectorXd& V, double dt)
{
    if (node_ <= 0) return;
 
    double V_terminal_kV = V(node_ - 1);
 
    switch (type_) {
    case NodeType::PV:
        solvePVNode(V_terminal_kV, dt);
        break;
    case NodeType::SLACK:
        solveSlackNode(dt);
        break;
    }
 
    calculateInjectionCurrent();
}
 
void Generator::solvePVNode(double V_terminal_kV, double dt)
{
    // PV节点：已知P和|V|，求解Q和θ
    // 使用牛顿-拉夫逊法迭代求解
    
    double V_pu = V_terminal_kV / V_base_kV_;
    double P_pu = P_set_MW_ / S_base_MVA_;
    
    // 初始猜值
    double theta_rad = theta_calc_deg_ * M_PI / 180.0;
    double delta_rad = delta_deg_ * M_PI / 180.0;
    
    for (int iter = 0; iter < MAX_ITER_; ++iter) {
        // 计算功率方程残差
        double cos_delta = std::cos(delta_rad);
        double sin_delta = std::sin(delta_rad);
        
        // P = (E*V/Xd)*sin(δ) - (V²/Xd - V²/Xq)*sin(2θ)/2
        double P_calc = (E_internal_pu_ * V_pu / Xd_pu_) * sin_delta;
        double P_error = P_calc - P_pu;
        
        // 电压幅值约束：|V| = V_set
        double V_error = V_pu - (V_set_kV_ / V_base_kV_);
        
        // 检查收敛
        if (std::abs(P_error) < CONV_TOL_ && std::abs(V_error) < CONV_TOL_) {
            break;
        }
        
        // 雅可比矩阵求解（简化处理）
        double dP_ddelta = (E_internal_pu_ * V_pu / Xd_pu_) * cos_delta;
        
        if (std::abs(dP_ddelta) > 1e-12) {
            delta_rad -= P_error / dP_ddelta;
        }
        
        // 更新功角
        delta_deg_ = delta_rad * 180.0 / M_PI;
        if (delta_deg_ > 180.0) delta_deg_ -= 360.0;
        if (delta_deg_ < -180.0) delta_deg_ += 360.0;
    }
    
    // 计算无功功率
    double cos_delta = std::cos(delta_rad);
    Q_calc_MVAr_ = ((E_internal_pu_ * V_pu / Xd_pu_) * cos_delta - V_pu * V_pu / Xd_pu_) * S_base_MVA_;
    
    // 更新相角（这里简化处理，实际需要考虑网络约束）
    theta_calc_deg_ = 0.0; // 参考节点相角
    P_calc_MW_ = P_set_MW_;
}
 
void Generator::solveSlackNode(double dt)
{
    // 平衡节点：已知|V|和θ，计算P和Q
    // 直接设置电压约束，功率由网络决定
    
    double V_pu = V_set_kV_ / V_base_kV_;
    theta_calc_deg_ = theta_set_deg_;
    
    // 这里简化处理，实际需要通过网络方程求解功率
    // 在电磁暂态仿真中，平衡节点通常作为无穷大母线处理
    double theta_rad = theta_calc_deg_ * M_PI / 180.0;
    double delta_rad = delta_deg_ * M_PI / 180.0;
    
    // 计算功率输出
    double cos_delta = std::cos(delta_rad);
    double sin_delta = std::sin(delta_rad);
    
    P_calc_MW_ = ((E_internal_pu_ * V_pu / Xd_pu_) * sin_delta) * S_base_MVA_;
    Q_calc_MVAr_ = ((E_internal_pu_ * V_pu / Xd_pu_) * cos_delta - V_pu * V_pu / Xd_pu_) * S_base_MVA_;
}
 
void Generator::calculateInternalVoltage()
{
    // 根据端电压和电流计算内部电势
    // 这里使用简化的稳态模型
    double V_pu = V_set_kV_ / V_base_kV_;
    double P_pu = P_calc_MW_ / S_base_MVA_;
    double Q_pu = Q_calc_MVAr_ / S_base_MVA_;
    
    if (V_pu > 1e-6) {
        double I_pu = std::sqrt(P_pu * P_pu + Q_pu * Q_pu) / V_pu;
        double phi_rad = std::atan2(Q_pu, P_pu);
        
        // E = V + jX*I（简化为直轴分量）
        E_internal_pu_ = std::sqrt(V_pu * V_pu + (Xd_pu_ * I_pu) * (Xd_pu_ * I_pu) 
                                  + 2 * V_pu * Xd_pu_ * I_pu * std::sin(phi_rad));
    } else {
        E_internal_pu_ = 1.0;
    }
}
 
void Generator::calculateInjectionCurrent()
{
    // 计算注入电流
    double V_pu = (node_ > 0) ? V_set_kV_ / V_base_kV_ : 1.0;
    double S_pu = std::sqrt(P_calc_MW_ * P_calc_MW_ + Q_calc_MVAr_ * Q_calc_MVAr_) / S_base_MVA_;
    
    if (V_pu > 1e-6) {
        double I_pu = S_pu / V_pu;
        I_inject_kA_ = I_pu * S_base_MVA_ / (std::sqrt(3.0) * V_base_kV_); // 转换为kA
    } else {
        I_inject_kA_ = 0.0;
    }
}
 
void Generator::setPVConstraints(double P_MW, double V_kV)
{
    if (type_ == NodeType::PV) {
        P_set_MW_ = P_MW;
        V_set_kV_ = V_kV;
    }
}
 
void Generator::setSlackConstraints(double V_kV, double theta_deg)
{
    if (type_ == NodeType::SLACK) {
        V_set_kV_ = V_kV;
        theta_set_deg_ = theta_deg;
    }
}