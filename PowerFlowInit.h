#pragma once

#include "PowerFlow.h"

#include <Eigen/Dense>
#include <vector>

namespace PFInit {

struct Source {
    int bus { 0 };
    int internal { 0 };
    double V_kV { 0.0 };
    double V_A { 0.0 };
    double R { 0.0 };
    double L { 0.0 };
};

struct Xfmr {
    int bus1 { 0 };
    int bus2 { 0 };
    int internal { 0 };
    double ratio { 1.0 };
    double R1 { 0.0 };
    double R2 { 0.0 };
    double L2 { 0.0 };
    double Rm { 0.0 };
    double Lm { 0.0 };
};

struct Load {
    int bus { 0 };
    double L { 0.0 };
};

struct Line {
    int from { 0 };
    int to { 0 };
    Eigen::Matrix3d R;
    Eigen::Matrix3d L;
    Eigen::Vector3d Cp;
    double Cg { 0.0 };
};

struct Inductor {
    int from { 0 };
    int to { 0 };
    double L { 0.0 };
};

struct Capacitor {
    int from { 0 };
    int to { 0 };
    double C { 0.0 };
};

struct Data {
    double freq { 50.0 };
    int nodeCount { 0 };
    std::vector<Source> sources;
    std::vector<Xfmr> xfmrs;
    std::vector<Load> loads;
    std::vector<Line> lines;
    std::vector<Inductor> inductors;
    std::vector<Capacitor> capacitors;
};

struct LineState {
    Eigen::Vector3d series;
    Eigen::Vector3d phase_i;
    Eigen::Vector3d phase_j;
    double ground_i { 0.0 };
    double ground_j { 0.0 };
};

struct State {
    Eigen::VectorXd V;
    std::vector<double> source;
    std::vector<Eigen::Vector3d> load;
    std::vector<LineState> line;
    std::vector<double> inductor;
    std::vector<double> capacitor;
};

Data makeData(double freq);
State init(const Result& flow, const Data& data, double t, double dt);

}