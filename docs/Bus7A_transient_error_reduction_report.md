# IEEE9 C++ 暂态误差由约 8 kV 降至约 80 V 的改进报告

## 1. 报告目的

本文说明 IEEE9 电磁暂态 C++ 模型与 Simulink/Powergui 模型对比时，如何将 BUS7A 电压的暂态最大绝对误差从 kV 级降低到几十 V 级，并解释各项修改降低误差的机理。

报告中的“约 8 kV 降至约 80 V”是工程化表述。现有 `summary.json` 给出的精确结果为：

| 指标 | 改进前候选模型 | 最终过零模型 | 改善 |
|---|---:|---:|---:|
| BUS7A 最大绝对误差 | `7.202744975 kV` | `0.075644573 kV` | 降低 `95.22` 倍 |
| BUS7A 最大绝对误差（V） | `7202.745 V` | `75.645 V` | 降低 `98.95%` |
| BUS7A RMSE | `0.887272162 kV` | `0.005753788 kV` | 降低 `154.21` 倍 |
| BUS7A RMSE（V） | `887.272 V` | `5.754 V` | 降低 `99.35%` |
| 最大误差步号 | `3,006,769` | `3,005,938` | 均位于故障切除后的恢复暂态 |

最终计算使用 `dt = 1 us`、`0~4 s`、`4,000,001` 个样本。以 230 kV 系统的相电压基准 `230/sqrt(3) = 132.79 kV` 计算，最终最大误差约为 `0.057%`，RMSE 约为 `0.0043%`。

改进前证据：

```text
reports/emtp_final_candidate_20260715/comparison/summary.json
```

最终证据：

```text
reports/emtp_zero_cross_final_20260715/comparison/summary.json
```

## 2. 问题本质

原始 kV 级误差不是单一参数错误造成的，而是以下误差叠加：

1. C++ 与 Powergui 的离散时间步、事件步号和故障拓扑切换时刻不一致；
2. PI-line 的相域串联矩阵、并联电容拓扑和历史源递推与 SPS 不一致；
3. 通用变压器由多个独立器件拼接，其动态状态定义与 SPS 变压器不同；
4. C++ 从零状态或预热状态启动，而 Simulink 使用精确的 `t=-1 us` 节点电压和器件历史状态；
5. 故障命令清除时刻与交流断路器实际电流过零开断时刻被混为同一时刻；
6. 比较输出中物理故障电流、外部协议电流和拓扑状态没有明确分离。

故障期间 BUS7A 被 `0.001 ohm` 故障支路钳位，电压接近零，因此两侧误差本来就较小。最大误差集中在故障切除后的恢复阶段：一旦故障支路释放，线路电容、线路串联电感、变压器和电源内部保存的历史状态同时参与恢复。如果这些状态或拓扑切换时刻略有差异，误差会在开断后几毫秒内被放大到 kV 级。

因此，本次工作的重点不是对 BUS7A 输出做补偿，而是让 C++ 的离散网络方程、历史状态和事件顺序与 Simulink/Powergui 对齐。

## 3. 改进总览

| 类别 | 原模型 | 改进后 |
|---|---|---|
| 仿真步长 | 旧入口曾使用 `50 us` | 与 Powergui 一致的 `1 us` |
| 事件定位 | 浮点时间比较 | 以整数步号定义关键协议 |
| PI-line 串联支路 | 通用序参数换算和旧历史递推 | 精确相域矩阵及直接梯形伴随递推 |
| L1/L3 并联电容 | `Cabc` 直接接地 | 相电容—浮动中性点—接地电容结构 |
| 变压器 | 漏抗、理想变比和励磁支路拼接 | 每相二状态 SPS 等效伴随模型 |
| 初始状态 | 零状态或预热近似 | 63 节点、99 历史量的 Q40 精确预测状态 |
| 故障开关 | `TwoValueResistor + FaultEvent` | `CircuitBreakerPhase` 自然过零开断 |
| 故障输出 | 命令和物理状态混合 | 命令、拓扑、物理电流分别输出 |
| 验证方式 | 只看整段曲线 | 样本数、关键步号、RMSE 和峰值同时验收 |

## 4. 时间步和事件协议对齐

### 4.1 统一为 1 us 固定步长

最终模型使用：

```cpp
control.dt = 1e-6;
control.t_end = 4.0;
```

这样 C++ 与 Powergui 在每一个微秒上都有一一对应的离散步。最终输出应为：

```text
step 0          -> 0.000000 s
step 2,000,001  -> 2.000001 s
step 3,000,000  -> 3.000000 s
step 4,000,000  -> 4.000000 s
```

如果 C++ 使用 `50 us` 而 Simulink 使用 `1 us`，断路器过零点最多可能相差几十个微秒；在故障恢复阶段，这种偏差足以产生明显的瞬时电压差。

### 4.2 关键事件使用整数步号

最终协议为：

```cpp
constexpr std::int64_t fault_start_step = 2'000'001;
constexpr std::int64_t command_clear_step = 3'000'000;
```

关键动作不再依赖 `abs(t-event.time)<dt/2` 这种可能受浮点表示影响的判断。输出比较也按照行号重建 C++ 时间轴，避免文本时间列有效数字不足造成重复时间戳或 `10 us` 跳变。

该项修改消除了“模型本身正确，但比较点错位”的伪误差。

## 5. PI-line 离散模型修正

### 5.1 精确串联相域矩阵

BUS7 通过 L1、L3 分别连接 BUS4、BUS5，因此这两条线路的误差直接进入 BUS7 节点方程。最终模型优先将 L1、L3 改为显式 `Rabc/Labc`：

```text
L1: BUS4 -- BUS7
L3: BUS5 -- BUS7
```

精确相域矩阵保留了 SPS 提取参数的完整有效数字和相间耦合，避免序参数舍入后再转换回 ABC 相域产生差异。

串联 RL 方程为：

```text
v(t) = R*i(t) + L*di(t)/dt
```

使用梯形积分后写成诺顿形式：

```text
i(k) = Gseries*v(k) + h(k)
Gseries = inverse(Rabc + 2*Labc/dt)
```

历史源直接递推为：

```text
A = Gseries*(2*Labc/dt - Rabc)
B = (I + A)*Gseries
h(k) = A*h(k-1) + B*v(k-1)
```

对应代码位于：

```text
Devices.cpp
PI_line::computeMatrices()
PI_line::updateHistory()
```

直接递推历史源避免了“先算物理电流、再反推等效历史源”的多次变换，使运算顺序更接近 Powergui 的伴随模型。

### 5.2 L1/L3 浮动中性点电容

旧模型将 `Cabc` 等效矩阵直接接地。L1/L3 的目标物理结构为：

```text
A ---- Cp_A ----+
B ---- Cp_B ----+---- neutral ---- Cgnd ---- ground
C ---- Cp_C ----+
```

内部中性点通过 KCL 消元：

```text
Gphase = 2*Cphase/dt
Gground = 2*Cground/dt
denom = sum(Gphase) + Gground
Gsh = diag(Gphase) - Gphase*Gphase^T/denom
```

每一步还要根据上一拍三相电压恢复内部中性点电压，再更新电容历史电流。这样保留了零序暂态下的中性点偏移及三相耦合。

L2/L4/L5/L6 暂时保留 aggregate 序参数构造。原因不是它们物理上不需要精确模型，而是本轮目标为 BUS7A 且要求尽量小改动；L1/L3 直接连接 BUS7，优先精确化的收益最大。

## 6. 变压器改为每相二状态 SPS 等效模型

旧通用变压器由以下部分组合：

```text
原边漏抗 + 理想变比约束 + 副边漏抗 + 励磁支路
```

这种连接在物理含义上合理，但内部未知量、状态数量、励磁历史量和离散运算顺序与 Simulink/SPS 不相同。只替换 `Rm/Lm` 数值而保留旧递推，仍不能消除暂态误差。

最终物理分支为每台变压器使用精确参数，并将每相表示为二状态系统：

```text
dx/dt = F*x + K*v
i = H*x + J*v
```

梯形离散后：

```text
P = inverse(I - dt*F/2)*(I + dt*F/2)
Q = inverse(I - dt*F/2)*(dt*K/2)
```

再转换为节点法可直接使用的诺顿形式：

```text
i(k) = Y*v(k) + w(k)
Y = H*Q + J
w(k) = A*w(k-1) + B*v(k-1)
A = H*P*inverse(H)
B = H*(P*Q + Q)
```

实现位于：

```text
transformer.cpp
Transformer::computePhysicalTwoStateParameters()
Transformer::updateHistory()
Transformer::updateState()
```

这一修改统一了变压器的状态定义和伴随递推，是降低故障恢复误差的核心改进之一。

## 7. 装载精确的 t=-1 us 初始状态

### 7.1 原方法的问题

零状态启动会令线路电容、线路电感、变压器励磁支路、负载电感和电源内电感的历史源全部从零开始。负时间预热可以接近稳态，但不能保证每个内部状态与 Powergui 完全一致。

即使 `t=0` 的端口电压接近，只要隐藏历史状态不同，故障切除时储能释放过程仍会不同。

### 7.2 最终初始状态夹具

新增 `IEEE9PhysicalInitialState.h/.cpp`，保存并装载：

```text
63 个 t=-1 us 节点预测电压
99 个器件历史状态
9 个电源 Norton 历史量
线路、变压器和负载的分组历史量
```

数据采用 Q40 格式：

```text
double_value = q40_integer / 2^40
```

初始化时按设备类型明确映射：

```text
VoltageSource -> 电源内部电感历史量
Transformer   -> 每相两个历史源及上一拍端口电压
PI_line       -> 串联历史量、i 端电容历史量、j 端电容历史量
Load          -> 三相负载电感历史量
```

`Simulation::run(initial_predictor_voltage)` 在进入 `step 0` 前先用 `t=-dt` 电压更新器件上一拍状态，然后才计算第 0 步历史注入和节点电压。

这项修改解决的是“隐藏状态不一致”，而不仅是“初始波形看起来不一致”。从 kV 级进一步下降到几十 V 级，依赖模型方程与完整初值同时一致；只做其中一项不足以达到最终结果。

## 8. 故障支路改为自然过零开断

### 8.1 从 TwoValueResistor 切换到 CircuitBreakerPhase

原 IEEE9 故障支路由 `TwoValueResistor` 和 `FaultEvent` 直接改变电阻状态，在 `3.0 s` 命令清除时立即开路。该方式没有交流电流过零过程。

最终故障支路使用：

```cpp
CircuitBreakerPhase(
    fault_node, 0,
    1e-3, 1e9,
    false, freq);
```

动作过程为：

```text
step 2,000,001：断路器立即合闸，投入故障
step 3,000,000：发出 CloseToOpen 命令
step 3,000,000 之后：保持闭合并检查电流符号
step 3,003,332：故障电流自然过零，实际开断
```

过零判断为：

```text
i(k-1)*i(k) <= 0
或 abs(i(k)) <= I_EPS
```

检测到过零后改变断路器导纳并重建、重新分解节点导纳矩阵。

### 8.2 命令状态与物理拓扑分离

新增输出：

```text
fault_command_active
fault_topology_closed
fault_I_physical
fault_switch
```

因此在 `3.000000~3.003332 s`：

```text
fault_command_active  = 0
fault_topology_closed = 1
```

这表示外部命令已经撤销，但物理断路器仍在等待自然过零。比较 Simulink 协议电流时使用命令窗口；验证电气开断时使用物理电流和拓扑状态，避免把两个定义混为一谈。

最终最大 BUS7A 误差发生在 `step 3,005,938`，即实际开断后约 `2.606 ms`。这说明剩余峰值是开断后的恢复暂态差异，不是把 3.0 s 命令时刻误当作实际开断时刻。

## 9. 电源、负荷和输出合同对齐

除三类核心动态元件外，还进行了以下配套修正：

1. 系统频率统一为 `60 Hz`；
2. 三组电源使用精确峰值、初相角和内电感；
3. 负荷参数与目标 Simulink 模型一致；
4. 电源与负载增加 `setInitialHistory()`，使初值夹具能进入实际内部电感；
5. `Curve` 增加可访问当前时间和节点电压的派生电流输出；
6. C++ 电压输出按 kV、比较脚本按 `1000` 转换为 V；
7. C++ 时间轴按 `row_index*1 us` 重建，不使用精度不足的文本时间列。

这些修改保证最终比较是在同一单位、同一时间步和同一输出定义下进行。

## 10. 阶段性结果

现有阶段报告显示误差下降不是一次参数调整完成的：

| 阶段报告 | BUS7A RMSE | BUS7A 最大绝对误差 | 峰值步号 |
|---|---:|---:|---:|
| `emtp_sim_transformer_breaker_comparison` | `6.888619 kV` | `50.078991 kV` | `3,003,865` |
| `emtp_sim_physical_line_layout_comparison` | `0.608474 kV` | `29.020460 kV` | `868` |
| `emtp_final_candidate_20260715` | `0.887272 kV` | `7.202745 kV` | `3,006,769` |
| `emtp_low_error_final_20260715` | `0.005754 kV` | `0.075720 kV` | `3,005,948` |
| `emtp_zero_cross_final_20260715` | `0.005754 kV` | `0.075645 kV` | `3,005,938` |

这些阶段不是严格的单变量消融实验，中间同时存在参数、拓扑和初值合同变化，因此不能把每一行差值完全归因于某一个函数。但它们能证明以下事实：

- 仅替换断路器或单个元件模型，不能直接得到 80 V 级结果；
- 线路拓扑和精确离散化可以把几十 kV 误差显著降低；
- 变压器状态定义和完整初始历史量决定能否从数 kV 进一步降到几十 V；
- 最后的自然过零逻辑将 `75.720 V` 小幅优化为 `75.645 V`，它主要修正开断物理时序，而不是独自贡献全部百倍改善。

## 11. 为什么最终仍有约 75.6 V 峰值

最终模型仍不是 Simulink 内核的逐指令复刻，剩余误差来源包括：

1. C++ 和 Powergui 的矩阵装配、求解顺序及浮点舍入不同；
2. 过零点在 `1 us` 网格上检测，没有对零点进行子步插值；
3. 断路器在本步求解后检测过零，开路拓扑从下一求解步生效；
4. L2/L4/L5/L6 仍使用 aggregate 序参数模型；
5. 输出文本有效数字会贡献亚伏级量化误差；
6. 故障切除后储能网络对极小历史状态差异非常敏感。

但误差在开断后快速衰减，并回到约几 V 的稳定周期残差，没有持续增长。这表明最终模型是稳定的，剩余误差属于离散实现差异，而不是参数或拓扑仍存在数量级错误。

## 12. 本次改进没有采用的方法

为保证结果具有物理意义，本次没有使用：

- 根据 MATLAB golden 波形反向拟合自由参数；
- 对 BUS7A 输出直接加补偿量；
- 使用 POD、ARX 或其他数据驱动输出修正；
- 任意调整断路器电阻、snubber 或阈值来压低峰值；
- 将 `Rm`、`Lm` 错误串联；
- 混用中间 103 状态夹具和最终 99 状态夹具。

误差下降来自离散模型、物理参数、初始状态和事件合同对齐，而不是对最终输出做后处理。

## 13. 结论

BUS7A 暂态最大绝对误差从 `7.202745 kV` 降至 `75.645 V`，不是某一个公式单独产生的结果，而是以下闭环共同作用：

```text
1 us 时间轴和整数步事件
        +
PI-line 精确相域伴随模型
        +
L1/L3 浮动中性点电容
        +
变压器每相二状态 SPS 等效
        +
63 节点/99 历史量精确初值
        +
CircuitBreakerPhase 自然过零开断
        +
命令、拓扑、物理电流分离验证
```

其中，线路和变压器修改解决“方程不一致”，Q40 初始状态解决“隐藏储能状态不一致”，断路器修改解决“故障切除时刻不一致”。三类问题全部对齐后，故障恢复阶段才从 kV 级误差下降到约 `80 V`。

## 14. 复核文件

主要源码：

```text
main.cpp
Devices.h
Devices.cpp
transformer.h
transformer.cpp
Simulation.h
Simulation.cpp
Curve.h
Curve.cpp
IEEE9PhysicalInitialState.h
IEEE9PhysicalInitialState.cpp
```

模型修正说明：

```text
docs/IEEE9_model_correction_guide.md
```

改进前结果：

```text
reports/emtp_final_candidate_20260715/comparison/summary.json
```

最终结果：

```text
reports/emtp_zero_cross_final_20260715/comparison/summary.json
reports/emtp_zero_cross_final_20260715/curve_V.dat
reports/emtp_zero_cross_final_20260715/curve_I.dat
```
