#pragma once

#include <Eigen/Dense>
#include <complex>
#include <vector>

struct Branch {
    int from { 0 };
    int to { 0 };
    std::complex<double> Y { 0.0, 0.0 };
};

struct Input {
    double freq { 50.0 };
    Eigen::MatrixXcd Y;
    Eigen::VectorXcd I;
    std::vector<Branch> lines;
};

struct Result {
    Eigen::VectorXcd V;
    Eigen::VectorXd V_abs;
    Eigen::VectorXd V_A;
    Eigen::VectorXd P;
    Eigen::VectorXd Q;
    std::vector<std::complex<double>> I_line;
    std::vector<std::complex<double>> S_line;
    double residual { 0.0 };
};

Input makeInput();
Result solve(const Input& in);