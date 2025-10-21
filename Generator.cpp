/* =============================================================================
 * Generator.cpp - 同步发电机模型实现
 * -----------------------------------------------------------------------------
 * PV节点和平衡节点发电机的电磁暂态仿真实现
 * ============================================================================= */
 
#include "Generator.h"
#include <cmath>
#include <iostream>
// --- Generator ---
Generator::Generator(int n_pos, int n_neu, NodeType type,
        double P_MW, double V_kV, double phase_deg,
        double V_base_kV, double S_base_MVA)
        : n_pos_(n_pos)
        , n_neu_(n_neu)
        , type_(type)
        , P_set_MW_(P_MW)
        , V_set_kV_(V_kV)
        , phase_set_deg_(phase_deg)
        , V_base_kV_(V_base_kV)
        , S_base_MVA_(S_base_MVA) 
{}

void Generator::stamp(std::vector<Eigen::Triplet<double>>& triplets, double dt) {
    // 仅添加数值阻尼，避免矩阵奇异（类似电流源处理）
    if (n_pos_ > 0) {
        const double g_damp = 1e-12;  // 极小电导，不影响潮流
        triplets.emplace_back(n_pos_ - 1, n_pos_ - 1, g_damp);
    }
}

void Generator::updateHistory(Eigen::VectorXd& I, double t, double dt) {
    if (n_pos_ <= 0) return;
    // 注入计算好的电流（电流源特性）
    I(n_pos_ - 1) += I_inject_kA_;
}

void Generator::updateState(const Eigen::VectorXd& V, double dt) {
    if (type_ == NodeType::PV) {
        // PV节点：从电压向量中获取终端电压
        double v_pos = (n_pos_ > 0) ? V(n_pos_ - 1) : 0.0;
        double v_neu = (n_neu_ > 0) ? V(n_neu_ - 1) : 0.0;
        double V_terminal_kV = std::abs(v_pos - v_neu);
        
        solvePVNode(V_terminal_kV);
    } 
    else if (type_ == NodeType::SLACK) {
        // 平衡节点：基于系统求解后的电压计算功率
        solveSlackNode(V);
    }
    
    calculateInjectionCurrent();
}
 
void Generator::solveSlackNode(const Eigen::VectorXd& V) {
    // 平衡节点的电压幅值和相角固定为设定值
    phase_calc_deg_ = phase_set_deg_;
    
    // 从电压向量获取实际节点电压
    double v_pos = (n_pos_ > 0) ? V(n_pos_ - 1) : 0.0;
    double v_neu = (n_neu_ > 0) ? V(n_neu_ - 1) : 0.0;
    double v_terminal = v_pos - v_neu;  // 终端电压(kV)
    
    // 计算注入电流（这是关键的修正）
    // 注入电流 = 终端电压 / 内阻 （这里简化为理想电压源）
    // 实际情况下应该根据发电机的等效电路计算
    
    // 方法1：基于功率平衡的简化计算
    // 对于理想平衡节点，假设发电机能提供任意所需功率
    // 通过KCL计算注入电流
    
    // 获取与此节点相连的所有设备电流（需要Grid提供接口）
    calculateSlackPowerFromSystemState(V);
}
 
void Generator::calculateSlackPowerFromSystemState(const Eigen::VectorXd& V) {
    // 这是修正的核心函数
    // 平衡节点的功率 = 系统总负荷 - 其他发电机功率 + 网损
    
    // 方法：通过节点电压和系统导纳矩阵计算注入电流
    double v_pos = (n_pos_ > 0) ? V(n_pos_ - 1) : 0.0;
    double v_neu = (n_neu_ > 0) ? V(n_neu_ - 1) : 0.0;
    
    // 将电压转换为复数形式
    double V_mag = std::abs(v_pos - v_neu);
    double V_angle_rad = phase_set_deg_ * M_PI / 180.0;
    
    std::complex<double> V_complex(V_mag * cos(V_angle_rad), 
                                  V_mag * sin(V_angle_rad));
    
    // 简化版本：假设发电机内阻很小，注入电流主要由外部负载决定
    // 这需要知道与该节点相连的所有支路电流
    
    // 临时简化计算（应该由Grid类提供更精确的接口）
    std::complex<double> I_injection = calculateInjectionCurrent(V);
    
    // 计算复功率 S = V * I*
    std::complex<double> S_complex = V_complex * std::conj(I_injection);
    
    // 转换为MW和MVAr
    double S_base_pu = S_base_MVA_ * 1000.0;  // 转换为kVA
    P_calc_MW_ = S_complex.real() * S_base_pu / 1000.0;      // MW
    Q_calc_MVAr_ = S_complex.imag() * S_base_pu / 1000.0;    // MVAr
}
 
std::complex<double> Generator::calculateInjectionCurrent(const Eigen::VectorXd& V) {
    // 这个函数计算节点的总注入电流
    // 注意：这需要访问系统的导纳矩阵，理想情况下应该由Grid类提供
    
    // 简化版本：假设我们能够通过某种方式获取连接到该节点的支路电流
    // 实际实现需要Grid类提供相应的接口
    
    // 对于现在的简化实现，我们基于功率设定值估算
    if (type_ == NodeType::SLACK) {
        // 对于平衡节点，这里应该通过 I = Y * V 计算
        // 但由于接口限制，暂时返回基于电压的估算值
        
        double v_pos = (n_pos_ > 0) ? V(n_pos_ - 1) : 0.0;
        double v_neu = (n_neu_ > 0) ? V(n_neu_ - 1) : 0.0;
        double v_terminal = v_pos - v_neu;
        
        // 这里需要更精确的实现，应该访问Grid的导纳矩阵
        // 暂时使用简化的估算
        std::complex<double> I_estimated(0.1, 0.05);  // 临时值，需要实际计算
        return I_estimated;
    }
    
    return std::complex<double>(0.0, 0.0);
}