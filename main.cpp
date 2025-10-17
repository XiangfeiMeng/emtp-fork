// ============================================================================
// 本工程的全局单位体系（务必统一）：电压[kV]，电流[kA]，阻抗[Ω]，电感[H]，电容[F]，功率[MW/MVAr]
// Node 编号：1-based（0 为地）；仿真步长建议 50 μs（工频场景），更高速场景可酌情减小。
// 三个算例：
//   1) test_simple() —— 110 kV 简单系统 + 单相永久故障/重合闸；验证断路器过零开断逻辑
//   2) test_converter() —— 0.38 kV 三段馈线 + 3 台变流器（PLL/VSG）；附源幅值/频率阶跃事件
//   3) test_two_transformers_inrush() —— 两台 110/35 kV Yd±30° 变压器；观察励磁涌流与相位差影响
// 输出：curve_V.dat（kV），curve_I.dat（kA），curve_S.dat（MW/MVAr，每 0.01 s）
// ============================================================================

#include "Control.h"
#include "Curve.h"
#include "Devices.h"
#include "Grid.h"
#include "Simulation.h"
#include "converter.h"
#include "Generator.h"
#include "transformer.h"
#include <iostream>
#include <memory>
#include <vector>

/* ------------------- 算例说明：110 kV 简单系统 -------------------
 * 拓扑：三相正弦源(1/2/3) -> 断路器(7/8/9) -> PI 线路 -> 负荷母线(4/5/6)
 * 事件：
 *   t=2.00 s：在节点 4 施加永久故障（TwoValueResistor 闭合到地）
 *   t=2.10 s：A 相断路器分闸（等待过零开断）
 *   t=3.10 s：A 相重合闸
 *   t=3.20 s：三相断路器分闸（合于故障失败加速解列）
 * 记录：母线(4/5/6)电压、三相电源电流
 * 现象：A 相短路引起不对称电压/电流；分闸时刻以“电流过零”为准，而非指令时刻。
 */
void test_simple()
{
    std::cout << "--- EMTP Simulation for Simple Case ---\n";
    Control control;
    control.dt = 50e-6;
    control.t_end = 3.3;
    control.plot_t_start = 1.9;
    control.plot_t_end = 3.3;

    // 记录节点 4/5/6 电压
    control.traces = { { "voltage_A", 4 }, { "voltage_B", 5 }, { "voltage_C", 6 } };
    // 2.0 s 施加永久故障（节点4）
    control.events.push_back(FaultEvent { 2.0, 4, true });

    // —— 以下内容与原始算例一致（略做变量名保持） ——
    const double freq = 50.0;
    const double v_line = 110.0;
    const double line_length_km = 5.0;
    const double R1_per_km = 0.2, R0_per_km = 0.6;
    const double X1_per_km = 0.26, X0_per_km = 1.13;
    const double C1_per_km = 0.012e-6, C0_per_km = 0.006e-6;

    const double R1 = R1_per_km * line_length_km;
    const double R0 = R0_per_km * line_length_km;
    const double L1 = (X1_per_km * line_length_km) / (2 * 3.14159 * freq);
    const double L0 = (X0_per_km * line_length_km) / (2 * 3.14159 * freq);
    const double C1 = C1_per_km * line_length_km;
    const double C0 = C0_per_km * line_length_km;

    Grid grid(9);

    double v_peak = v_line * sqrt(2.0 / 3.0);
    double rs = 0.001, isc_peak = 48;
    auto srcA = std::make_unique<VoltageSource>(1, 0, v_peak, freq, 0.0, rs, isc_peak);
    auto srcB = std::make_unique<VoltageSource>(2, 0, v_peak, freq, -120.0, rs, isc_peak);
    auto srcC = std::make_unique<VoltageSource>(3, 0, v_peak, freq, +120.0, rs, isc_peak);
    auto srcA_ptr = srcA.get(), srcB_ptr = srcB.get(), srcC_ptr = srcC.get();
    grid.addDevice(std::move(srcA));
    grid.addDevice(std::move(srcB));
    grid.addDevice(std::move(srcC));

    const double R_closed = 1e-5, R_open = 1e9;
    auto brkA = std::make_unique<CircuitBreakerPhase>(1, 7, R_closed, R_open, true, freq);
    auto brkB = std::make_unique<CircuitBreakerPhase>(2, 8, R_closed, R_open, true, freq);
    auto brkC = std::make_unique<CircuitBreakerPhase>(3, 9, R_closed, R_open, true, freq);
    auto brkA_ptr = brkA.get(), brkB_ptr = brkB.get(), brkC_ptr = brkC.get();
    grid.addDevice(std::move(brkA));
    grid.addDevice(std::move(brkB));
    grid.addDevice(std::move(brkC));

    auto line = std::make_unique<PI_line>(
        std::vector<int> { 7, 8, 9 }, std::vector<int> { 4, 5, 6 },
        R1, L1, C1, R0, L0, C0, freq, control.dt);
    grid.addDevice(std::move(line));

    auto load = std::make_unique<Load>(std::vector<int> { 4, 5, 6 }, 60, 0, v_line, freq, 0);
    grid.addDevice(std::move(load));

    const int fault_node = 4;
    auto fault_switch = std::make_unique<TwoValueResistor>(fault_node, 0, 1e-3, 1e9);
    TwoValueResistor* fault_switch_ptr = fault_switch.get();
    grid.addDevice(std::move(fault_switch));

    Curve curve(control, fault_switch_ptr);
    curve.addCurrentTrace("I_src_A", [srcA_ptr]() { return srcA_ptr->get_I(); });
    curve.addCurrentTrace("I_src_B", [srcB_ptr]() { return srcB_ptr->get_I(); });
    curve.addCurrentTrace("I_src_C", [srcC_ptr]() { return srcC_ptr->get_I(); });

    Simulation simulation(control, grid, curve, fault_switch_ptr, fault_node);
    simulation.addBreaker(brkA_ptr);
    simulation.addBreaker(brkB_ptr);
    simulation.addBreaker(brkC_ptr);
    simulation.addVoltageSource(srcA_ptr);
    simulation.addVoltageSource(srcB_ptr);
    simulation.addVoltageSource(srcC_ptr);

    const double t_trip_A = 2.10;
    const double t_recloseA = 3.10;
    const double t_trip_3ph = 3.20;
    brkA_ptr->schedule(t_trip_A, CircuitBreakerPhase::Change::CloseToOpen);
    brkA_ptr->schedule(t_recloseA, CircuitBreakerPhase::Change::OpenToClose);
    brkA_ptr->schedule(t_trip_3ph, CircuitBreakerPhase::Change::CloseToOpen);
    brkB_ptr->schedule(t_trip_3ph, CircuitBreakerPhase::Change::CloseToOpen);
    brkC_ptr->schedule(t_trip_3ph, CircuitBreakerPhase::Change::CloseToOpen);

    simulation.run();
}

/* -------------- 算例说明：0.38 kV 馈线 + 三台并网变流器 --------------
 * 拓扑：源(1/2/3) -> 3 段 0.5 km PI 线路 -> 母线段(4-6),(7-9),(10-12)
 *      每段节点处挂 R//L 负荷；在 4/5/6、7/8/9、10/11/12 分别接入三台 Converter。
 * 模式：默认 PLL；如需 VSG，将 if_VSG 置 true。
 * 事件：
 *   2~5 s：源线电压幅值缩放至 0.88 p.u.（电压跌落）
 *   5~8 s：系统频率偏移 -0. Hz
 * 记录：三台变流器 P/Q（MW/MVAr）与三相电流；各段末端电压。
 * 现象：PLL 与 VSG 在低电压/频偏下的有功/无功能力与动态差异；内置电流限幅。
 */
void test_converter()
{
    std::cout << "--- EMTP Simulation for test_converter ---\n";

    Control control;
    control.dt = 50e-6;
    control.t_end = 10.0; // 仿真时间
    control.plot_t_start = 1.5; // 录波启动时间
    control.plot_t_end = 10.0; // 录波结束时间

    const double freq = 50.0;
    const double v_line = 0.38;
    const double line_length_km = 0.5;
    const double R1_per_km = 0.02, R0_per_km = 0.03;
    const double X1_per_km = 0.26, X0_per_km = 1.13;
    const double C1_per_km = 0.012e-6, C0_per_km = 0.006e-6;

    const double R1 = R1_per_km * line_length_km;
    const double R0 = R0_per_km * line_length_km;
    const double L1 = (X1_per_km * line_length_km) / (2 * 3.14159 * freq);
    const double L0 = (X0_per_km * line_length_km) / (2 * 3.14159 * freq);
    const double C1 = C1_per_km * line_length_km;
    const double C0 = C0_per_km * line_length_km;

    Grid grid(12);

    // 电源
    double v_peak = 1.05 * v_line * std::sqrt(2.0 / 3.0);
    double rs = 0.0001, isc_peak = 12;
    auto srcA = std::make_unique<VoltageSource>(1, 0, v_peak, freq, 0.0, rs, isc_peak);
    auto srcB = std::make_unique<VoltageSource>(2, 0, v_peak, freq, -120., rs, isc_peak);
    auto srcC = std::make_unique<VoltageSource>(3, 0, v_peak, freq, +120., rs, isc_peak);
    auto srcA_ptr = srcA.get(), srcB_ptr = srcB.get(), srcC_ptr = srcC.get();
    grid.addDevice(std::move(srcA));
    grid.addDevice(std::move(srcB));
    grid.addDevice(std::move(srcC));

    // 线路
    grid.addDevice(std::make_unique<PI_line>(
        std::vector<int> { 1, 2, 3 }, std::vector<int> { 4, 5, 6 },
        R1, L1, C1, R0, L0, C0, freq, control.dt));
    grid.addDevice(std::make_unique<PI_line>(
        std::vector<int> { 4, 5, 6 }, std::vector<int> { 7, 8, 9 },
        R1, L1, C1, R0, L0, C0, freq, control.dt));
    grid.addDevice(std::make_unique<PI_line>(
        std::vector<int> { 7, 8, 9 }, std::vector<int> { 10, 11, 12 },
        R1, L1, C1, R0, L0, C0, freq, control.dt));

    // 负荷
    grid.addDevice(std::make_unique<Load>(std::vector<int> { 4, 5, 6 }, 0.015, 0.0, v_line, freq, 0));
    grid.addDevice(std::make_unique<Load>(std::vector<int> { 7, 8, 9 }, 0.015, 0.0, v_line, freq, 0));
    grid.addDevice(std::make_unique<Load>(std::vector<int> { 10, 11, 12 }, 0.015, 0.0, v_line, freq, 0));

    // 并网变流器
    bool if_VSG { false }; // 是否虚拟同步机模型
    auto conv1 = std::make_unique<Converter>(4, 5, 6, 0.01, 0.0);
    auto conv2 = std::make_unique<Converter>(7, 8, 9, 0.015, 0.0);
    auto conv3 = std::make_unique<Converter>(10, 11, 12, 0.005, 0.0);
    auto conv_ptr1 = conv1.get(), conv_ptr2 = conv2.get(), conv_ptr3 = conv3.get();
    grid.addDevice(std::move(conv1));
    grid.addDevice(std::move(conv2));
    grid.addDevice(std::move(conv3));

    // 曲线：功率采样（每 0.01s）
    Curve curve(control);
    curve.addPowerTrace("P1", [conv_ptr1]() { return conv_ptr1->get_P_MW(); }, "Q1", [conv_ptr1]() { return conv_ptr1->get_Q_MVAr(); });
    curve.addPowerTrace("P2", [conv_ptr2]() { return conv_ptr2->get_P_MW(); }, "Q2", [conv_ptr2]() { return conv_ptr2->get_Q_MVAr(); });
    curve.addPowerTrace("P3", [conv_ptr3]() { return conv_ptr3->get_P_MW(); }, "Q3", [conv_ptr3]() { return conv_ptr3->get_Q_MVAr(); });

    if (if_VSG) {
        conv_ptr1->set_control_mode(ControlMode::VSG);
        conv_ptr2->set_control_mode(ControlMode::VSG);
        conv_ptr3->set_control_mode(ControlMode::VSG);
    }

    curve.addCurrentTrace("Ia_1", [conv_ptr1]() { return conv_ptr1->get_Ia(); });
    curve.addCurrentTrace("Ib_1", [conv_ptr1]() { return conv_ptr1->get_Ib(); });
    curve.addCurrentTrace("Ic_1", [conv_ptr1]() { return conv_ptr1->get_Ic(); });

    curve.addCurrentTrace("Ia_2", [conv_ptr2]() { return conv_ptr2->get_Ia(); });
    curve.addCurrentTrace("Ib_2", [conv_ptr2]() { return conv_ptr2->get_Ib(); });
    curve.addCurrentTrace("Ic_2", [conv_ptr2]() { return conv_ptr2->get_Ic(); });

    curve.addCurrentTrace("Ia_3", [conv_ptr3]() { return conv_ptr3->get_Ia(); });
    curve.addCurrentTrace("Ib_3", [conv_ptr3]() { return conv_ptr3->get_Ib(); });
    curve.addCurrentTrace("Ic_3", [conv_ptr3]() { return conv_ptr3->get_Ic(); });

    control.addVoltageTrace("10_A", 10);
    control.addVoltageTrace("11_B", 11);
    control.addVoltageTrace("12_C", 12);

    control.addVoltageTrace("7_A", 7);
    control.addVoltageTrace("8_B", 8);
    control.addVoltageTrace("9_C", 9);

    control.addVoltageTrace("4_A", 4);
    control.addVoltageTrace("5_B", 5);
    control.addVoltageTrace("6_C", 6);

    // 仿真器
    Simulation simulation(control, grid, curve);
    simulation.addVoltageSource(srcA_ptr);
    simulation.addVoltageSource(srcB_ptr);
    simulation.addVoltageSource(srcC_ptr);

    // 下发电压/频率时程（-1 表示对全部三相源生效）
    control.events.push_back(VoltageRampEvent { 2.0, -1, 0.88 }); // 2s→5s：电压-10%
    control.events.push_back(VoltageRampEvent { 5.0, -1, 1.00 });

    control.events.push_back(FrequencyOffsetEvent { 5.0, -1, -0.1 }); // 5s→8s：系统频率-0.1 Hz
    control.events.push_back(FrequencyOffsetEvent { 8.0, -1, 0.00 });

    simulation.run();
}

/* ------- 算例说明：两台 110/35 kV Yd±30° 变压器励磁涌流 -------
 * T1：Yd +30°，高压侧通过相间断路器空投（2.005/2.007/2.009 s 分相合闸），低压侧开路；
 * T2：Yd -30°，长期接在 4/5/6（低压侧 7/8/9 带负荷 ~20 MW）。
 * 模型要点：
 *   - 漏抗由短路试验参数折算，分摊至一、二侧串联；
 *   - 励磁采用两段法（线性段/空心段）；在段切换时仅做数值重分解，保持稀疏结构稳定；
 * 记录：T1 高压侧三相断路器电流（观察合闸瞬时与衰减过程中的涌流峰值与波形）。
 */
void test_two_transformers_inrush()
{
    std::cout << "--- EMTP Simulation: Two Transformers Inrush ---\n";

    Control control;
    control.dt = 50e-6;
    control.t_end = 3.5;
    control.plot_t_start = 1.0;
    control.plot_t_end = 3.5;

    // 电源侧（110 kV）三相
    const double freq = 50.0;
    const double V_HV = 110.0; // kV
    Grid grid(15);

    // 源（接到 1/2/3）
    double v_peak = V_HV * std::sqrt(2.0 / 3.0);
    double rs = 0.001, isc_peak = 48;
    auto srcA = std::make_unique<VoltageSource>(1, 0, v_peak, freq, 0.0, rs, isc_peak);
    auto srcB = std::make_unique<VoltageSource>(2, 0, v_peak, freq, -120.0, rs, isc_peak);
    auto srcC = std::make_unique<VoltageSource>(3, 0, v_peak, freq, +120.0, rs, isc_peak);
    auto srcA_ptr = srcA.get(), srcB_ptr = srcB.get(), srcC_ptr = srcC.get();
    grid.addDevice(std::move(srcA));
    grid.addDevice(std::move(srcB));
    grid.addDevice(std::move(srcC));

    // 源到母线：简短线路（等效 2km）
    const double line_km = 2.0;
    const double R1_per_km = 0.2, R0_per_km = 0.6;
    const double X1_per_km = 0.26, X0_per_km = 1.13;
    const double C1_per_km = 0.012e-6, C0_per_km = 0.006e-6;
    auto add_line = [&](std::vector<int> ni, std::vector<int> nj) {
        double R1 = R1_per_km * line_km;
        double R0 = R0_per_km * line_km;
        double L1 = (X1_per_km * line_km) / (2 * 3.14159 * freq);
        double L0 = (X0_per_km * line_km) / (2 * 3.14159 * freq);
        double C1 = C1_per_km * line_km;
        double C0 = C0_per_km * line_km;
        grid.addDevice(std::make_unique<PI_line>(ni, nj, R1, L1, C1, R0, L0, C0, freq, control.dt));
    };
    add_line({ 1, 2, 3 }, { 4, 5, 6 });

    // 断路器用于 T1 高压侧空投
    const double R_closed = 5e-4, R_open = 1e7;
    auto brkA = std::make_unique<CircuitBreakerPhase>(4, 13, R_closed, R_open, /*closed*/ false, freq);
    auto brkB = std::make_unique<CircuitBreakerPhase>(5, 14, R_closed, R_open, /*closed*/ false, freq);
    auto brkC = std::make_unique<CircuitBreakerPhase>(6, 15, R_closed, R_open, /*closed*/ false, freq);
    auto brkA_ptr = brkA.get(), brkB_ptr = brkB.get(), brkC_ptr = brkC.get();
    // 把 13/14/15 当作 T1 高压侧接点（通过断路器与4/5/6相连）
    grid.addDevice(std::move(brkA));
    grid.addDevice(std::move(brkB));
    grid.addDevice(std::move(brkC));

    // —— 两台变压器参数（110/35 kV，50 MVA，Yd ±30°）——
    Transformer_para tp;
    tp.load = 50.0;
    tp.freq = 50.0;
    tp.swinging1Voltage = 110.0;
    tp.swinging2Voltage = 35.0;
    tp.saturationEnabled = 1; // 1-两段法;0-线性
    tp.airCoreReactance = 0.05;
    tp.kneeVoltage = 0.95;
    tp.noLoadLosses = 50.0;
    tp.noLoadpercentage = 1.0;
    tp.short_circuit_loss = 300.0;
    tp.short_circuit_percentage = 12.0;
    tp.swinging1type = 0; // Y
    tp.swinging1_GROUNDING = 1;
    tp.swinging2type = 1; // Δ
    tp.swinging2_GROUNDING = 0;
    tp.angle_delta = +30;
    tp.tapMid = 0;
    tp.tap_now = 0;
    tp.tap_max = +9;
    tp.tap_min = -9;
    tp.ratioTap = 1.25;

    // T1：高压侧通过断路器接 13/14/15；低压侧 10/11/12（开路，无负荷）
    auto T1 = std::make_unique<Transformer>(std::vector<int> { 13, 14, 15 }, std::vector<int> { 10, 11, 12 }, tp);
    grid.addDevice(std::move(T1));

    // T2：高压侧直接接 4/5/6；低压侧 7/8/9（带负荷）
    tp.angle_delta = -30; // 另一台取 -30°，便于观察合应涌流
    auto T2 = std::make_unique<Transformer>(std::vector<int> { 4, 5, 6 }, std::vector<int> { 7, 8, 9 }, tp);
    grid.addDevice(std::move(T2));

    // 低压侧负荷（35 kV 侧，三相合计 ~20 MW）
    grid.addDevice(std::make_unique<Load>(std::vector<int> { 7, 8, 9 }, 20.0, 2.0, 35.0, freq, 0));

    // 曲线

    control.addVoltageTrace("4_A", 4);
    control.addVoltageTrace("5_B", 5);
    control.addVoltageTrace("6_C", 6);

    control.addVoltageTrace("7_A", 7);
    control.addVoltageTrace("8_B", 8);
    control.addVoltageTrace("9_C", 9);

    Curve curve(control);

    curve.addCurrentTrace("I_T1_HV_A", [brkA_ptr]() { return brkA_ptr->get_I(); });
    curve.addCurrentTrace("I_T1_HV_B", [brkB_ptr]() { return brkB_ptr->get_I(); });
    curve.addCurrentTrace("I_T1_HV_C", [brkC_ptr]() { return brkC_ptr->get_I(); });

    // 仿真器
    Simulation sim(control, grid, curve);
    sim.addVoltageSource(srcA_ptr);
    sim.addVoltageSource(srcB_ptr);
    sim.addVoltageSource(srcC_ptr);
    sim.addBreaker(brkA_ptr);
    sim.addBreaker(brkB_ptr);
    sim.addBreaker(brkC_ptr);

    // 空投时刻：t=2.0 s 合上 T1 三相断路器
    brkA_ptr->schedule(2.005, CircuitBreakerPhase::Change::OpenToClose);
    brkB_ptr->schedule(2.007, CircuitBreakerPhase::Change::OpenToClose);
    brkC_ptr->schedule(2.009, CircuitBreakerPhase::Change::OpenToClose);

    sim.run();
}
void test_ieee9bus() {
    std::cout << "--- EMTP Simulation for IEEE 9-Bus System ---\n";

    // 1. 仿真控制参数设置
    Control control;
    control.dt = 50e-6;          // 50μs时间步长（适配短线路模型）
    control.t_end = 6.0;         // 总仿真时间6秒
    control.plot_t_start = 0.0;  // 从0秒开始录波
    control.plot_t_end = 4.0;    // 录波至4秒
    const double freq = 50.0;    // 系统频率50Hz
    const double base_kV = 230.0; // 高压侧基准线电压230kV
    const double base_MVA = 100.0; // 基准功率100MVA

    // 2. 建立电网模型（9个节点+接地节点0）
    Grid grid(27);  // 9节点系统含27个相节点（A/B/C各9个）

    // 3. 发电机模型
    // IEEE 9节点系统标准参数：
    // 节点1：平衡节点，V=1.04pu ∠0°
    // 节点2：PV节点，P=163MW, V=1.025pu
    // 节点3：PV节点，P=85MW, V=1.025pu
 
    double V1_kV = 1.04 * base_kV / sqrt(3.0);   // 节点1相电压（平衡节点）
    double V2_kV = 1.025 * base_kV / sqrt(3.0);  // 节点2相电压（PV节点）
    double V3_kV = 1.025 * base_kV / sqrt(3.0);  // 节点3相电压（PV节点）
 
    // 创建三相发电机组
    // 节点1 - 平衡节点发电机组（A/B/C三相）
    auto gen1A = std::make_unique<Generator>(1, Generator::NodeType::SLACK,
        0.0, V1_kV, 0.0,           // P初值=0MW, V=相电压, θ=0°
        0.0576, 0.0576, 0.0,       // Xd=Xq=0.0576pu, Ra=0 (简化)
        base_kV / sqrt(3.0), base_MVA);
 
    auto gen1B = std::make_unique<Generator>(2, Generator::NodeType::SLACK,
        0.0, V1_kV, -120.0,        // B相滞后120°
        0.0576, 0.0576, 0.0,
        base_kV / sqrt(3.0), base_MVA);
 
    auto gen1C = std::make_unique<Generator>(3, Generator::NodeType::SLACK,
        0.0, V1_kV, 120.0,         // C相超前120°
        0.0576, 0.0576, 0.0,
        base_kV / sqrt(3.0), base_MVA);
 
    // 节点2 - PV节点发电机组（A/B/C三相）
    double P2_MW = 163.0 / 3.0;  // 每相功率（三相均分）
    auto gen2A = std::make_unique<Generator>(4, Generator::NodeType::PV,
        P2_MW, V2_kV, 0.0,         // P=54.33MW/相, V=相电压
        0.1292, 0.1292, 0.0,       // Xd=Xq=0.1292pu
        base_kV / sqrt(3.0), base_MVA);
 
    auto gen2B = std::make_unique<Generator>(5, Generator::NodeType::PV,
        P2_MW, V2_kV, -120.0,
        0.1292, 0.1292, 0.0,
        base_kV / sqrt(3.0), base_MVA);
 
    auto gen2C = std::make_unique<Generator>(6, Generator::NodeType::PV,
        P2_MW, V2_kV, 120.0,
        0.1292, 0.1292, 0.0,
        base_kV / sqrt(3.0), base_MVA);
 
    // 节点3 - PV节点发电机组（A/B/C三相）
    double P3_MW = 85.0 / 3.0;   // 每相功率（三相均分）
    auto gen3A = std::make_unique<Generator>(7, Generator::NodeType::PV,
        P3_MW, V3_kV, 0.0,         // P=28.33MW/相, V=相电压
        0.1813, 0.1813, 0.0,       // Xd=Xq=0.1813pu
        base_kV / sqrt(3.0), base_MVA);
 
    auto gen3B = std::make_unique<Generator>(8, Generator::NodeType::PV,
        P3_MW, V3_kV, -120.0,
        0.1813, 0.1813, 0.0,
        base_kV / sqrt(3.0), base_MVA);
 
    auto gen3C = std::make_unique<Generator>(9, Generator::NodeType::PV,
        P3_MW, V3_kV, 120.0,
        0.1813, 0.1813, 0.0,
        base_kV / sqrt(3.0), base_MVA);

    // 保存发电机指针用于录波
    auto gen1A_ptr = gen1A.get();
    auto gen1B_ptr = gen1B.get();
    auto gen1C_ptr = gen1C.get();
    auto gen2A_ptr = gen2A.get();
    auto gen2B_ptr = gen2B.get();
    auto gen2C_ptr = gen2C.get();
    auto gen3A_ptr = gen3A.get();
    auto gen3B_ptr = gen3B.get();
    auto gen3C_ptr = gen3C.get();

    // 添加发电机到电网
    grid.addDevice(std::move(gen1A));
    grid.addDevice(std::move(gen1B));
    grid.addDevice(std::move(gen1C));
    grid.addDevice(std::move(gen2A));
    grid.addDevice(std::move(gen2B));
    grid.addDevice(std::move(gen2C));
    grid.addDevice(std::move(gen3A));
    grid.addDevice(std::move(gen3B));
    grid.addDevice(std::move(gen3C));

    // 4. 变压器模型 (3台双绕组变压器)
    Transformer_para t_para;
    t_para.load = 100.0;          // 100MVA容量
    t_para.swinging1Voltage = 13.8; // 发电机侧13.8kV
    t_para.swinging2Voltage = 230.0; // 系统侧230kV
    t_para.swinging1type = 0;     // Y接
    t_para.swinging2type = 1;     // Δ接
    t_para.short_circuit_percentage = 10.0;
    t_para.freq = freq;

    // 变压器1: 节点1(A/B/C) -> 节点4(A/B/C)
    grid.addDevice(std::make_unique<Transformer>(
        std::vector<int>{1, 2, 3}, std::vector<int>{10, 11, 12}, t_para));
    // 变压器2: 节点2(A/B/C) -> 节点5(A/B/C)
    grid.addDevice(std::make_unique<Transformer>(
        std::vector<int>{4, 5, 6}, std::vector<int>{13, 14, 15}, t_para));
    // 变压器3: 节点3(A/B/C) -> 节点6(A/B/C)
    grid.addDevice(std::make_unique<Transformer>(
        std::vector<int>{7, 8, 9}, std::vector<int>{16, 17, 18}, t_para));

    // 5. 输电线路模型 (6条分布参数线路)
    // 线路参数 (每公里值)
    const double R1_per_km = 0.026;    // 正序电阻(Ω/km)
    const double X1_per_km = 0.26;     // 正序电抗(Ω/km)
    const double C1_per_km = 0.012e-6; // 正序电容(F/km)
    const double R0_per_km = 0.08;     // 零序电阻(Ω/km)
    const double X0_per_km = 0.8;      // 零序电抗(Ω/km)
    const double C0_per_km = 0.006e-6; // 零序电容(F/km)

    // 线路拓扑: 起点(A/B/C)、终点(A/B/C)、长度(km)
    std::vector<std::tuple<std::vector<int>, std::vector<int>, double>> lines = {
        {{10, 11, 12}, {19, 20, 21}, 40.0},  // 线路4-7
        {{10, 11, 12}, {25, 26, 27}, 20.0},  // 线路4-9
        {{13, 14, 15}, {19, 20, 21}, 20.0},  // 线路5-7
        {{13, 14, 15}, {22, 23, 24}, 20.0},  // 线路5-8
        {{16, 17, 18}, {22, 23, 24}, 20.0},  // 线路6-8
        {{16, 17, 18}, {25, 26, 27}, 30.0}   // 线路6-9
    };

    for (const auto& line : lines) {
        const auto& from_nodes = std::get<0>(line);
        const auto& to_nodes = std::get<1>(line);
        double length = std::get<2>(line);
        
        // 计算线路总参数
        double R1 = R1_per_km * length;
        double L1 = (X1_per_km * length) / (2 * M_PI * freq);
        double C1 = C1_per_km * length;
        double R0 = R0_per_km * length;
        double L0 = (X0_per_km * length) / (2 * M_PI * freq);
        double C0 = C0_per_km * length;

        // 使用混合Bergeron-PI-ZIM模型
        grid.addDevice(std::make_unique<PI_line>(
            from_nodes, to_nodes, R1, L1, C1, R0, L0, C0, freq, control.dt));
    }

    // 6. 负荷模型 (3个三相负荷)
    // 负荷1: 节点7(A/B/C)
    grid.addDevice(std::make_unique<Load>(
        std::vector<int>{19, 20, 21}, 100.0, 30.0, base_kV, freq, 0));
    // 负荷2: 节点8(A/B/C)
    grid.addDevice(std::make_unique<Load>(
        std::vector<int>{22, 23, 24}, 100.0, 30.0, base_kV, freq, 0));
    // 负荷3: 节点9(A/B/C)
    grid.addDevice(std::make_unique<Load>(
        std::vector<int>{25, 26, 27}, 100.0, 30.0, base_kV, freq, 0));

    // 7. 故障设置 (2.0秒7A相接地，3.0秒切除)
    const int fault_node = 19; // 节点7A相
    auto fault_switch = std::make_unique<TwoValueResistor>(fault_node, 0, 1e9, 0.001); // 节点7A相
    TwoValueResistor* fault_switch_ptr = fault_switch.get();
    grid.addDevice(std::move(fault_switch));
    control.events.push_back(FaultEvent{2.0, fault_node, true});   // 施加故障
    control.events.push_back(FaultEvent{3.0, fault_node, false});  // 移除故障
    // 8. 录波设置
    Curve curve(control,fault_switch_ptr);
    // 记录关键节点电压
    control.addVoltageTrace("Bus7A", 19);
    // control.addVoltageTrace("Bus7B", 20);
    // control.addVoltageTrace("Bus7C", 21);
    // control.addVoltageTrace("Bus8A", 22);
    // 记录发电机电流
    curve.addCurrentTrace("Gen1A", [gen1A_ptr]() { return gen1A_ptr->get_I(); });
    
    // curve.addCurrentTrace("Gen2A", [gen2A_ptr]() { return gen2A_ptr->get_I(); });
    // curve.addCurrentTrace("Gen3A", [gen3A_ptr]() { return gen3A_ptr->get_I(); });


    // 9. 启动仿真
    Simulation simulation(control, grid, curve,fault_switch_ptr,fault_node);
    // 添加电压源到仿真器
    // simulation.addVoltageSource(gen1A_ptr);
    // simulation.addVoltageSource(gen1B_ptr);
    // simulation.addVoltageSource(gen1C_ptr);
    // simulation.addVoltageSource(gen2A_ptr);
    // simulation.addVoltageSource(gen2B_ptr);
    // simulation.addVoltageSource(gen2C_ptr);
    // simulation.addVoltageSource(gen3A_ptr);
    // simulation.addVoltageSource(gen3B_ptr); 
    // simulation.addVoltageSource(gen3C_ptr);
    simulation.run();
}


int main()
{
    // test_converter();
    // test_two_transformers_inrush();
    // test_simple();
    test_ieee9bus();
    return 0;
}