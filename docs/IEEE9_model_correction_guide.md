# IEEE9 EMTP 模型修正与迁移指南

## 1. 目标

本文说明如何把 IEEE9 项目中已经验证的物理候选模型迁移到 `EMTP_1`，使 C++ 节点电压法模型尽量复现 Simulink/Powergui 的行为，同时保留现有旧模型用于回退和对照。

本次迁移不是调参，也不是读取 MATLAB 输出波形进行拟合。所有参数和状态均来自 Simulink 元件配置、Powergui 初始状态或 IEEE9 工程已经导出的物理模型数据。

目标模型的关键结果如下：

- 时间步长：`dt = 1e-6 s`。
- 样本数：`4,000,001`，对应 `0–4 s`。
- 节点数：`63`。
- 最终历史状态数：`99`。
- 故障投入拓扑起点：`step_idx = 2,000,001`。
- fault command 清除：`step_idx = 3,000,000`。
- 断路器实际开断：`step_idx = 3,003,332`。
- 故障闭合导纳：`1000 S`。

已验证的 double 模型对 MATLAB 误差：

| 信号 | RMSE | 最大绝对误差 |
|---|---:|---:|
| Bus7A | 约 `0.006078 kV` | 约 `0.07051 kV` |
| Gen1A | 约 `0.000908 kV` | 约 `0.002939 kV` |
| fault_I | 约 `0.000327 kA` | 约 `0.003975 kA` |

## 2. 当前 C++ 模型的主要差异

### 2.1 初始状态不一致

当前 `test_ieee9bus()` 主要依靠负时间预热获得初始状态。Simulink 模型使用 Powergui 在 `t=-1 us` 的精确初始节点电压和器件状态。

预热只能得到近似稳态，无法保证变压器、线路电容、负载和电源的所有历史电流与 Powergui 一致。

### 2.2 变压器状态定义不同

当前通用 `Transformer` 根据空载试验和短路试验参数重新计算绕组及励磁参数，并包含非线性励磁切段。新候选模型使用 Simulink 的准确参数和每相两个动态状态：

1. 二次绕组历史电流；
2. 励磁支路历史电流。

`Rm` 与 `Lm` 是并联励磁支路，不应改成串联支路。

### 2.3 L1/L3 并联电容拓扑不同

旧 `PI_line` 将 `Cabc` 等效矩阵直接接地。L1/L3 的物理模型应为每端：

```text
A ---- Cp_A ----+
B ---- Cp_B ----+---- neutral ---- Cgnd ---- ground
C ---- Cp_C ----+
```

内部 neutral 可以通过 KCL 消去，不必增加全局节点。

### 2.4 断路器开断过早

当前故障开关通常在 `3.0 s` command 清除时立即打开。Simulink 中 command 在第 `3,000,000` 步清除，但 A 相实际 status 到第 `3,003,332` 步才打开，并且该点与闭合开关电流过零对齐。

因此必须分开保存：

- `fault_command_active`；
- `fault_topology_closed`。

## 3. 推荐的代码组织

不要直接删除旧模型。建议增加并行入口：

```text
test_ieee9bus()           原模型
test_ieee9bus_physical()  新物理模型
```

通用 `Transformer` 和原 `PI_line` 构造函数继续保留。新模型通过重载构造函数或 IEEE9 专用类实现。

## 4. PI-line 修正

### 4.1 保留原 aggregate 构造函数

L2/L4/L5/L6 暂时继续使用原来的正序/零序参数构造函数和 aggregate `Cabc`。

L1/L3 使用新的相域矩阵构造函数：

```cpp
PI_line(std::vector<int> nodes_i, std::vector<int> nodes_j,
    const Eigen::Matrix3d& Rabc,
    const Eigen::Matrix3d& Labc,
    const Eigen::Vector3d& Cphase,
    double Cground,
    double dt);
```

类中增加：

```cpp
bool physical_neutral_caps = false;
Eigen::Vector3d Cphase = Eigen::Vector3d::Zero();
double Cground = 0.0;
Eigen::Vector3d Gphase = Eigen::Vector3d::Zero();
double Gground = 0.0;
```

### 4.2 精确串联矩阵

L1/L3 使用：

```cpp
Eigen::Matrix3d Rabc;
Rabc <<
    13.5293831769454, 12.2635478388026, 12.2635478388026,
    12.2635478388026, 13.5293831769454, 12.2635478388026,
    12.2635478388026, 12.2635478388026, 13.5293831769454;

Eigen::Matrix3d Labc;
Labc <<
    0.198642487743882, 0.105535006650437, 0.105535006650437,
    0.105535006650437, 0.198642487743882, 0.105535006650437,
    0.105535006650437, 0.105535006650437, 0.198642487743882;
```

串联 companion 保持：

```text
Gseries = inverse(Rabc + 2*Labc/dt)
Alpha   = (I - dt/2*inverse(Labc)*Rabc)
          * inverse(I + dt/2*inverse(Labc)*Rabc)
Beta    = (I + Alpha)*Gseries
h_new   = Alpha*h_old + Beta*v_old
```

### 4.3 物理中性点电容

每端参数：

```cpp
const Eigen::Vector3d Cphase =
    Eigen::Vector3d::Constant(0.637898940812743e-6);
const double Cground = 2.99136421769046e-6;
```

离散导纳：

```cpp
Gphase = (2.0 / dt) * Cphase;
Gground = (2.0 / dt) * Cground;
const double denom = Gphase.sum() + Gground;

G_sh_half = Gphase.asDiagonal()
    - (Gphase * Gphase.transpose()) / denom;
```

每一步先根据上一拍相电压计算内部 neutral：

```cpp
const double neutral = Gphase.dot(v_prev) / denom;
const Eigen::Vector3d branch_voltage =
    v_prev.array() - neutral;

h_new = -h_old
    - 2.0 * Gphase.cwiseProduct(branch_voltage);
```

节点注入仍采用 `I_phase -= h_new`。每端只保存三个相历史状态，不再单独保存 `Cgnd` 历史状态。

### 4.4 L1/L3 创建方式

```cpp
grid.addDevice(std::make_unique<PI_line>(
    std::vector<int>{10, 11, 12},
    std::vector<int>{19, 20, 21},
    Rabc, Labc, Cphase, Cground, control.dt)); // L1

grid.addDevice(std::make_unique<PI_line>(
    std::vector<int>{13, 14, 15},
    std::vector<int>{19, 20, 21},
    Rabc, Labc, Cphase, Cground, control.dt)); // L3
```

必须从旧六线路循环中删除 L1/L3，避免重复 stamp。

## 5. 变压器修正

### 5.1 参数

| 变压器 | V1 | R1 | L1 | V2 | R2 | L2 | Rm | Lm |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| T1 | 16500 | 2.7225e-6 | 0 | 230000 | 0.000529 | 0.080825 | 1361.2 | 3.6108 |
| T2 | 18000 | 3.24e-6 | 0 | 230000 | 0.000529 | 0.087701 | 1620.0 | 4.2972 |
| T3 | 13800 | 1.9044e-6 | 0 | 230000 | 0.000529 | 0.082228 | 952.2 | 2.5258 |

### 5.2 动态状态

每相采用两个状态。推荐直接实现已经验证的 injection-state 形式：

```text
w_new = A*w_old + B*v_old
Ihist = w_new
```

其中：

```text
A = H*P*inverse(H)
B = H*(P*Q + Q)
```

不要只替换 `Rm/Lm` 数值而保留旧状态递推；这样仍不能与 Simulink 对齐。第一版应关闭当前模型中的非线性饱和切段，先验证线性 Powergui 等效模型。

变压器的精确分组和系数位于：

```text
/mnt/d/Desktop/code/IEEE9/Code/reports/sim/
physical_stage_q40_injection_20260713/fixture/
transformer_predictor_groups.csv
```

## 6. 精确初始状态

### 6.1 删除预热依赖

最终模型不应依赖：

```cpp
control.t_start = -warmup_time;
```

应从 `step_idx=0` 开始，并在仿真前设置：

- 63 个节点的 `t=-1 us` 电压；
- 99 个器件历史状态；
- 9 个电源的 Norton 初始历史；
- 电源正弦/余弦振荡器状态。

### 6.2 器件初始化接口

建议为下列器件增加 `setInitialState()`：

- `VoltageSource`；
- `Inductor` / `series_RL`；
- `Load`；
- `PI_line`；
- IEEE9 精确变压器。

`Simulation::run()` 应接收初始节点电压，不能只执行 `V.setZero()` 后依赖预热。

注意：中间诊断模型曾使用 103 个历史状态；最终 KCL 消元模型使用 99 个历史状态。迁移时以最终 fixture 为准，不能混用两套下标。

## 7. 断路器和故障时序

### 7.1 使用整数步事件

不要使用浮点条件：

```cpp
std::abs(t - event.time) < dt / 2.0
```

使用 `step_idx`：

```cpp
constexpr std::int64_t FAULT_START_STEP = 2'000'001;
constexpr std::int64_t COMMAND_CLEAR_STEP = 3'000'000;
constexpr std::int64_t BREAKER_OPEN_STEP = 3'003'332;
```

第一版确定性规则：

```cpp
fault_command_active =
    step_idx >= FAULT_START_STEP &&
    step_idx < COMMAND_CLEAR_STEP;

fault_topology_closed =
    step_idx >= FAULT_START_STEP &&
    step_idx < BREAKER_OPEN_STEP;
```

### 7.2 故障 stamp

```cpp
if (fault_topology_closed) {
    G(fault_node, fault_node) += 1000.0;
}
```

打开后应完全删除该支路，不建议使用 `R_open=1e9` 留下残余导纳。

完成固定 status replay 后，再验证基于电流过零的 FSM 是否也能稳定得到 `3,003,332`。不要在第一版同时修改阈值和开断规则。

## 8. fault_I 输出

同时输出：

```text
fault_I_physical
fault_I_exported
fault_command_active
fault_topology_closed
```

物理电流：

```cpp
fault_I_physical = fault_topology_closed
    ? 1000.0 * Bus7A
    : 0.0;
```

`fault_I_exported` 可继续遵守现有外部协议。比较 MATLAB 时必须明确使用的是哪一种电流，不能混用。

## 9. 实施顺序

1. 保留旧入口，新增 `test_ieee9bus_physical()`。
2. 将所有事件改为整数 `step_idx`。
3. 实现固定 `3,003,332` 开断和两种 fault_I 输出。
4. 修改 L1/L3 的精确串联矩阵和物理中性点电容。
5. 增加精确初始状态接口，删除新入口中的预热。
6. 实现变压器双状态 companion。
7. 先跑短窗和关键步检查。
8. 短窗通过后再跑完整 `0–4 s`。

不要在同一次测试中同时修改多个尚未验证的物理模块，否则无法判断误差改善来自哪一项。

## 10. 验收步骤

### 10.1 编译

```bash
cd ~/coding/EMTP_1
cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Release
cmake --build build --target emtp_sim -j
```

### 10.2 基础检查

- 输出必须正好有 `4,000,001` 行。
- 时间范围必须为 `0–4 s`。
- 故障拓扑必须在 `2,000,001` 投入。
- command 必须在 `3,000,000` 清除。
- breaker status 必须在 `3,003,332` 打开。
- 不应出现 NaN、Inf 或求解失败。

### 10.3 关键窗口

至少比较：

- 初始窗口：`step 0` 附近；
- 故障投入：`step 2,000,001`；
- 故障稳定：`step 2,500,000`；
- command 清除：`step 3,000,000`；
- 实际开断：`step 3,003,332`；
- 恢复：`step 3,500,000`。

### 10.4 完整误差目标

以 Simulink 为参考：

- Bus7A RMSE 应接近 `0.0061 kV`；
- Gen1A RMSE 应接近 `0.00091 kV`；
- fault_I RMSE 应接近 `0.00033 kA`；
- 最大误差不应重新集中在 `3.0–3.005 s` 并达到旧模型数量级。

## 11. 禁止的处理方式

- 不使用 MATLAB golden 波形反推自由参数。
- 不添加 POD、ARX 或输出修正。
- 不继续调 breaker 电阻或 snubber。
- 不把 `Rm` 和 `Lm` 错接为串联支路。
- 不把中间 103 状态 fixture 与最终 99 状态 fixture 混用。
- 不在确认 C++ double 模型前修改 RTL production ROM。

## 12. 参考数据

模型合同：

```text
/mnt/d/Desktop/code/IEEE9/Code/reports/sim/
full_physical_initial_state_20260713/
candidate_cpp_handoff_contract.json
```

最终 Q40/RTL 交接：

```text
/mnt/d/Desktop/code/IEEE9/Code/reports/sim/
physical_candidate_rtl_handoff_20260713/diagnosis.md
```

最终 99 状态 fixture：

```text
/mnt/d/Desktop/code/IEEE9/Code/reports/sim/
physical_stage_q40_injection_20260713/fixture/
```

已经验证的独立 C++ double 参考程序：

```text
/mnt/d/Desktop/code/IEEE9/Code/scripts/full_physical_cpp_replay.cpp
```

该参考程序用于确认新物理模型本身是否正确；最终仍应使用 `EMTP_1` 的手工器件类逐项复现相同结果。
