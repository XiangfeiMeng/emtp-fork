/* =============================================================================
 * 曲线记录（Curve）
 * -----------------------------------------------------------------------------
 * - 将每步/定期采样的电压、电流、功率输出为 TSV 文本文件：
 *     curve_V.dat：Time + 若干节点电压[kV]
 *     curve_I.dat：Time + 用户注册的电流[kA]
 *     curve_S.dat：Time + 成对的 P/Q（MW/MVAr），默认每 0.01 s 采样一次
 * - 只负责 I/O，不参与电网络求解。
 * ============================================================================= */

#pragma once
#include "Control.h"
#include "Devices.h"
#include <Eigen/Dense>
#include <fstream>
#include <functional>
#include <string>
#include <vector>

class Curve {
public:
    Curve(const Control& ctrl, const TwoValueResistor* fault_switch);
    Curve(const Control& ctrl);
    ~Curve();

    void sample(double t, const Eigen::VectorXd& V);

    // 注册电流曲线
    void addCurrentTrace(const std::string& name, std::function<double()> reader)
    {
        curr_traces.push_back({ name, std::move(reader) });
    }

    void addDerivedCurrentTrace(const std::string& name,
        std::function<double(double, const Eigen::VectorXd&)> reader)
    {
        curr_traces.push_back({ name, {}, std::move(reader) });
    }

    // 功率采样（每 0.01s）
    void addPowerTrace(const std::string& nameP, std::function<double()> P_reader,
        const std::string& nameQ, std::function<double()> Q_reader)
    {
        pwr_traces.push_back({ nameP, std::move(P_reader), nameQ, std::move(Q_reader) });
    }

private:
    struct CurrTrace {
        std::string name;
        std::function<double()> reader;
        std::function<double(double, const Eigen::VectorXd&)> derived_reader;
    };
    struct PwrTrace {
        std::string nameP;
        std::function<double()> P_reader;
        std::string nameQ;
        std::function<double()> Q_reader;
    };

    const Control& ctrl;
    const TwoValueResistor* fault_switch;

    std::ofstream fileV, fileI, fileS; // 新增 fileS
    bool header_written_V { false };
    bool header_written_I { false };
    bool header_written_S { false };

    std::vector<CurrTrace> curr_traces;
    std::vector<PwrTrace> pwr_traces;

    // 采样节流：功率每 0.01s 采样一次
    const double S_sample_T { 0.01 };
    double last_S_sample_t { -1e9 };
};
