/* =============================================================================
 * Simulation：主时间步驱动器
 * -----------------------------------------------------------------------------
 * 主循环步骤：
 *   1) 事件触发（故障/源幅值/频偏）
 *   2) 断路器定时动作（可能触发拓扑变更 → 重建 G 并数值分解）
 *   3) updateHistory() 汇总 I
 *   4) 解 G·V=I（KLU）
 *   5) updateState()
 *   6) 断路器过零开断（如触发 → 仅重建/重分解一次）
 *   7) 变压器励磁段切换（仅数值重分解 + 单次对齐）
 *   8) 采样
 * ============================================================================= */

#pragma once
#include "Control.h"
#include "Curve.h"
#include "Devices.h"
#include "Grid.h"
#include "KluSolverEigen.h"
#include <Eigen/Sparse>
#include <vector>

/**
 * @class Simulation
 * @brief 仿真器主类，负责驱动整个电磁暂态仿真过程。
 *
 * 管理主时间步循环，处理事件，调用 Grid 和求解器，并记录数据。
 */
class Simulation {
public:
    /**
     * @brief 构造仿真器。
     * @param ctrl 仿真控制参数对象。
     * @param grid 电路网络对象。
     * @param curve 曲线记录对象。
     * @param fault_switch 指向用于施加故障的开关元件的指针。
     * @param fault_node 故障施加的节点编号。
     */
    Simulation(Control& ctrl, Grid& grid, Curve& curve, TwoValueResistor* fault_switch, int fault_node);
    Simulation(Control& ctrl, Grid& grid, Curve& curve);

    /// @brief 启动并运行整个仿真。
    void run();
    void run(const Eigen::VectorXd& initial_predictor_voltage);

    /**
     * @brief 注册一个断路器到仿真器中，以便在循环中处理其时序逻辑。
     * @param brk 指向 CircuitBreakerPhase 对象的指针。
     */
    void addBreaker(CircuitBreakerPhase* brk) { breakers.push_back(brk); }

    void addVoltageSource(VoltageSource* src) { sources.push_back(src); } /// 注册电压源指针，用于在事件触发时修改幅值/频偏

private:
    Control& ctrl;
    Grid& grid;
    Curve& curve;
    TwoValueResistor* fault_switch { nullptr }; ///< 指向故障开关
    int fault_node { -1 }; ///< 故障节点

    /// 稀疏矩阵线性方程组求解器 (LU分解)
    KluSolverEigen solver;

    /// @brief 存储所有需要时序控制的断路器
    std::vector<CircuitBreakerPhase*> breakers;

    std::vector<VoltageSource*> sources; // 新增

};
