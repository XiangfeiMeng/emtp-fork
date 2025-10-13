#include "Simulation.h"
#include "transformer.h" // 为了轮询 Transformer::consumeRefactorFlag()
#include <chrono>
#include <iostream>
#include <sstream>
#include <variant>

Simulation::Simulation(Control& ctrl, Grid& grid, Curve& curve, TwoValueResistor* fault_switch, int fault_node)
    : ctrl(ctrl)
    , grid(grid)
    , curve(curve)
    , fault_switch(fault_switch)
    , fault_node(fault_node)
{
}

Simulation::Simulation(Control& ctrl, Grid& grid, Curve& curve)
    : ctrl(ctrl)
    , grid(grid)
    , curve(curve)
{
}

// 主循环关键步骤摘要（见 README“仿真流程”）：
// 1) 事件 → 2) 断路器定时动作 → 3) 历史注入 → 4) 解 G·V=I
// 5) 元件状态更新 → 6) 断路器过零开断 → 7) 变压器励磁段切换（仅数值重分解）→ 8) 采样
// 注：拓扑改变需“重建矩阵+factorize”；仅数值改变（如励磁段切换）可直接 refactor

void Simulation::run()
{
    std::cout << "Starting simulation..." << std::endl;
    auto start_time = std::chrono::high_resolution_clock::now();

    // 初始化
    grid.allocateInternalNodes();
    grid.buildMatrix(ctrl.dt);
    if (false) {
        grid.dumpMatrix("导纳矩阵"); // 调试用
    }

    solver.analyzePattern(grid.getG());
    solver.factorize(grid.getG());

    int num_steps = static_cast<int>(ctrl.t_end / ctrl.dt);
    Eigen::VectorXd V(grid.getNumNodes());
    V.setZero();
    Eigen::VectorXd I(grid.getNumNodes());
    I.setZero();

    std::cout << "----begin----simulation " << std::endl;
    for (int i = 0; i <= num_steps; ++i) {
        I.setZero();
        V.setZero();
        double t = i * ctrl.dt;
        bool matrix_changed = false; // 仅用于拓扑/开关改变

        // 处理预定事件（安全访问 variant，避免 bad_variant_access）
        for (const auto& ev : ctrl.events) {
            if (std::holds_alternative<FaultEvent>(ev)) {
                const auto& event = std::get<FaultEvent>(ev);
                if (event.node == this->fault_node && std::abs(t - event.time) < ctrl.dt / 2.0) {
                    if (fault_switch) {
                        std::cout << "Time: " << t << "s - Fault event on node " << event.node
                                  << ": " << (event.apply ? "ON" : "OFF") << std::endl;
                        fault_switch->setState(event.apply);
                        matrix_changed = true;
                    }
                }
            } else if (std::holds_alternative<VoltageRampEvent>(ev)) {
                const auto& e = std::get<VoltageRampEvent>(ev);
                if (std::abs(t - e.time) < ctrl.dt / 2.0) {
                    if (e.idx < 0) {
                        for (auto* s : sources)
                            if (s)
                                s->set_voltage_scale(e.scale);
                    } else if (e.idx >= 0 && e.idx < (int)sources.size() && sources[e.idx]) {
                        sources[e.idx]->set_voltage_scale(e.scale);
                    }
                    // 仅源等效电流改变，不改导纳；无需重建矩阵
                }
            } else if (std::holds_alternative<FrequencyOffsetEvent>(ev)) {
                const auto& e = std::get<FrequencyOffsetEvent>(ev);
                if (std::abs(t - e.time) < ctrl.dt / 2.0) {
                    if (e.idx < 0) {
                        for (auto* s : sources)
                            if (s)
                                s->set_freq_offset(e.df_Hz);
                    } else if (e.idx >= 0 && e.idx < (int)sources.size() && sources[e.idx]) {
                        sources[e.idx]->set_freq_offset(e.df_Hz);
                    }
                    // 同上，无需重建矩阵
                }
            }
        }

        // 断路器动作（只在动作时重建）
        for (auto* brk : breakers) {
            if (brk && brk->applyScheduledAt(t)) {
                matrix_changed = true;
            }
        }

        // 如需（拓扑或开断）则重建并数值分解
        if (matrix_changed) {
            grid.buildMatrix(ctrl.dt);
            solver.factorize(grid.getG());
        }

        // 更新历史注入 I
        grid.updateHistoryVector(I, t, ctrl.dt);

        // 求解 G*V=I
        V = solver.solve(I);

        //  更新设备状态（这里可能触发“励磁段切换”标志）
        grid.updateDeviceStates(V, ctrl.dt);

        //  断路器过零开断（若改变，重建一次）
        bool changed_post = false;
        for (auto* brk : breakers) {
            if (brk && brk->checkZeroCrossAndOpen()) {
                changed_post = true;
            }
        }
        if (changed_post) {
            grid.buildMatrix(ctrl.dt);
            solver.factorize(grid.getG());
        }

        // 非线性励磁段切换：仅数值重分解（单次内迭代以对齐当前步）
        if (Transformer::consumeRefactorFlag()) {
            // 用新等效导纳重建数值并再解一次当前步（单次内迭代）
            grid.buildMatrix(ctrl.dt); // 结构不变，仅数值变
            solver.factorize(grid.getG()); // 仅数值分解
            grid.updateHistoryVector(I, t, ctrl.dt);
            V = solver.solve(I);
            grid.updateDeviceStates(V, ctrl.dt);
        }

        // 采样
        curve.sample(t, V);
    }

    auto end_time = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end_time - start_time;
    std::cout << "Simulation finished in " << elapsed.count() << " seconds." << std::endl;
}
