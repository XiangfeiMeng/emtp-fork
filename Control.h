/* =============================================================================
 * 仿真控制与事件系统（Control）
 * -----------------------------------------------------------------------------
 * - 统一存放：时间步长 dt、仿真时长 t_end、录波区间，以及各类时序事件。
 * - 事件类型：
 *     FaultEvent：在指定节点施加/移除故障（通过双值电阻到地）
 *     VoltageRampEvent：对电压源施加幅值缩放（scale，-1 作用于全部源）
 *     FrequencyOffsetEvent：对电压源施加绝对频率偏移（df_Hz，-1 作用于全部源）
 * - 曲线：使用 (name, node) 记录电压；更多通道由 Curve::add*Trace 注册。
 * ============================================================================= */

#pragma once
#include <string>
#include <variant>
#include <vector>

using std::string;
/** 故障事件：与原义相同 **/
struct FaultEvent {
    double time; ///< 事件发生的时间 (s)
    int node; ///< 故障节点
    bool apply; ///< true=施加故障, false=移除
};

/** 线电压幅值按时间段线性偏移（百分比，正为抬升） **/
struct VoltageRampEvent {
    double time; ///< 动作时间(s)
    int idx; ///< 目标源索引（-1 表示全部）
    double scale; ///< 幅值缩放，例如 1.05
};

/** 频率按时间段偏移（绝对偏移，单位 Hz，负值即下降） **/
struct FrequencyOffsetEvent {
    double time; ///< 动作时间(s)
    int idx; ///< 目标源索引（-1 表示全部）
    double df_Hz; ///< 频偏(Hz)，如 -0.01
};

// 可扩展的事件集合
using Event = std::variant<FaultEvent, VoltageRampEvent, FrequencyOffsetEvent>;

class Control {
public:
    void addVoltageTrace(string name, int node); /// 追加一个电压录波通道
    ///  name 曲线名（列名）
    ///  node 节点号（1-based；0 代表地，将记为 0）
    // 基本仿真参数
    double dt = 50e-6;
    double t_end = 3.3;
    double plot_t_start = 1.9;
    double plot_t_end = 3.3;

    // 事件与曲线配置（由各算例在 main 内设置）
    std::vector<Event> events;
    std::vector<std::pair<std::string, int>> traces;
};
