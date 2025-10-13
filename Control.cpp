#include "Control.h"

// 追加电压通道到 traces；Curve::sample() 将按列输出到 curve_V.dat
void Control::addVoltageTrace(string name, int node)
{
    traces.emplace_back(name, node);
}