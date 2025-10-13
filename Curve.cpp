#include "Curve.h"
#include <iomanip>
#include <iostream>

Curve::Curve(const Control& ctrl, const TwoValueResistor* fault_switch)
    : ctrl(ctrl)
    , fault_switch(fault_switch)
{
    fileV.open("curve_V.dat");
    if (!fileV.is_open())
        std::cerr << "Error: Could not open curve_V.dat\n";
    fileI.open("curve_I.dat");
    if (!fileI.is_open())
        std::cerr << "Error: Could not open curve_I.dat\n";
    fileS.open("curve_S.dat");
    if (!fileS.is_open())
        std::cerr << "Error: Could not open curve_S.dat\n";
}

Curve::Curve(const Control& ctrl)
    : ctrl(ctrl)
{
    fileV.open("curve_V.dat");
    if (!fileV.is_open())
        std::cerr << "Error: Could not open curve_V.dat\n";
    fileI.open("curve_I.dat");
    if (!fileI.is_open())
        std::cerr << "Error: Could not open curve_I.dat\n";
    fileS.open("curve_S.dat");
    if (!fileS.is_open())
        std::cerr << "Error: Could not open curve_S.dat\n";
}

Curve::~Curve()
{
    if (fileV.is_open())
        fileV.close();
    if (fileI.is_open())
        fileI.close();
    if (fileS.is_open())
        fileS.close();
}

// 仅在 [plot_t_start, plot_t_end] 内记数据；P/Q 通道每 0.01 s 下采样一次（S_sample_T）
void Curve::sample(double t, const Eigen::VectorXd& V)
{
    // 电压
    if (t >= ctrl.plot_t_start && t <= ctrl.plot_t_end) {
        if (!ctrl.traces.empty() && fileV.is_open()) {
            if (!header_written_V) {
                fileV << "Time";
                for (const auto& tr : ctrl.traces)
                    fileV << "\t" << tr.first;
                fileV << "\n";
                header_written_V = true;
            }
            fileV << t;
            for (const auto& tr : ctrl.traces) {
                double v = 0.0;
                int node_idx = tr.second;
                if (node_idx > 0 && node_idx <= V.size())
                    v = V(node_idx - 1);
                fileV << "\t" << v;
            }
            fileV << "\n";
        }

        // 电流
        if (!curr_traces.empty() && fileI.is_open()) {
            if (!header_written_I) {
                fileI << "Time";
                for (const auto& ct : curr_traces)
                    fileI << "\t" << ct.name;
                fileI << "\n";
                header_written_I = true;
            }
            fileI << t;
            for (const auto& ct : curr_traces) {
                double i = 0.0;
                if (ct.reader)
                    i = ct.reader();
                fileI << "\t" << i;
            }
            fileI << "\n";
        }

        // 功率（每 0.01s）
        if (!pwr_traces.empty() && fileS.is_open()) {
            if (t - last_S_sample_t >= S_sample_T - 1e-12) {
                if (!header_written_S) {
                    fileS << "Time";
                    for (const auto& pt : pwr_traces) {
                        fileS << "\t" << pt.nameP << "\t" << pt.nameQ;
                    }
                    fileS << "\n";
                    header_written_S = true;
                }
                fileS << t;
                for (const auto& pt : pwr_traces) {
                    double P = pt.P_reader ? pt.P_reader() : 0.0;
                    double Q = pt.Q_reader ? pt.Q_reader() : 0.0;
                    fileS << "\t" << P << "\t" << Q;
                }
                fileS << "\n";
                last_S_sample_t = t;
            }
        }
    }
}
