#include "PowerFlow.h"

#include <algorithm>
#include <stdexcept>
namespace {
void stamp(Eigen::MatrixXcd& Y, int from, int to, std::complex<double> y)
{
    Y(from, from) += y;
    Y(to, to) += y;
    Y(from, to) -= y;
    Y(to, from) -= y;
}
}


Input makeInput()
{
    Input in;
    in.freq = 60.0;
    in.Y = Eigen::MatrixXcd::Zero(9, 9);
    in.I = Eigen::VectorXcd::Zero(9);
    constexpr double pi = 3.14159265358979323846;
    const double w = 2.0 * pi * in.freq;
    const double vll[] = { 17.16, 18.45, 14.145 };
    const double ang[] = { 0.0, 9.28001, 4.66475 };
    const double Ls[] = { 0.000456, 0.000189, 0.000128 };
    for (int k = 0; k < 3; ++k) {
        const std::complex<double> z(0.001, w * Ls[k]);
        const auto y = 1.0 / z;
        const auto e = std::polar(vll[k] / std::sqrt(3.0), (ang[k] - 90.0) * pi / 180.0);
        in.Y(k, k) += y;
        in.I(k) += y * e;
    }
    const auto yl = 1.0 / std::complex<double>(1.2658353381428, w * 0.093107481093445);
    const std::complex<double> ys(0.0, w * 0.637898940812743e-6);
    const int lines[][2] = { {3,6}, {3,8}, {4,6}, {4,7}, {5,7}, {5,8} };
    for (const auto& line : lines) {
        stamp(in.Y, line[0], line[1], yl);
        in.Y(line[0], line[0]) += ys;
        in.Y(line[1], line[1]) += ys;
        in.lines.push_back({ line[0], line[1], yl });
    }
    const double R[] = { 423.2, 529.0, 587.777777777778 };
    const double L[] = { 2.80643216318709, 4.00918880455298, 4.67738693864515 };
    for (int k = 0; k < 3; ++k)
        in.Y(6 + k, 6 + k) += 1.0 / R[k] + 1.0 / std::complex<double>(0.0, w * L[k]);
    const double v1[] = { 16.5, 18.0, 13.8 };
    const double v2[] = { 230.0, 230.0, 230.0 };
    const double r1[] = { 2.7225e-6, 3.24e-6, 1.9044e-6 };
    const double r2[] = { 0.000529, 0.000529, 0.000529 };
    const double l2[] = { 0.080825, 0.087701, 0.082228 };
    const double rm[] = { 1361.2, 1620.0, 952.2 };
    const double lm[] = { 3.6108, 4.2972, 2.5258 };
    for (int k = 0; k < 3; ++k) {
        const double n = v1[k] / v2[k];
        const auto yp = 1.0 / std::complex<double>(r1[k], 0.0);
        const auto yt = 1.0 / std::complex<double>(r2[k], w * l2[k]);
        const auto ym = 1.0 / rm[k] + 1.0 / std::complex<double>(0.0, w * lm[k]);
        const auto d = yp + ym + yt / (n * n);
        const auto ap = yp / d;
        const auto as = yt / (n * d);
        in.Y(k, k) += yp * (1.0 - ap);
        in.Y(k, 3 + k) -= yp * as;
        in.Y(3 + k, k) -= yt * ap / n;
        in.Y(3 + k, 3 + k) += yt * (1.0 - as / n);
    }
    return in;
}

Result solve(const Input& in)
{
    if (in.Y.rows() == 0 || in.Y.rows() != in.Y.cols() || in.I.size() != in.Y.rows())
        throw std::invalid_argument("power-flow input dimensions are invalid");

    Eigen::FullPivLU<Eigen::MatrixXcd> lu(in.Y);
    if (!lu.isInvertible())
        throw std::runtime_error("power-flow matrix is singular");

    Result out;
    out.V = lu.solve(in.I);
    const int n = out.V.size();
    out.V_abs.resize(n);
    out.V_A.resize(n);
    out.P.resize(n);
    out.Q.resize(n);

    const Eigen::VectorXcd I = in.Y * out.V;
    for (int k = 0; k < n; ++k) {
        out.V_abs(k) = std::abs(out.V(k));
        out.V_A(k) = std::arg(out.V(k));
        const auto S = out.V(k) * std::conj(I(k));
        out.P(k) = S.real();
        out.Q(k) = S.imag();
    }

    out.I_line.reserve(in.lines.size());
    out.S_line.reserve(in.lines.size());
    for (const auto& line : in.lines) {
        if (line.from < 0 || line.to < 0 || line.from >= n || line.to >= n)
            throw std::invalid_argument("power-flow branch bus index is invalid");
        const auto I_line = line.Y * (out.V(line.from) - out.V(line.to));
        out.I_line.push_back(I_line);
        out.S_line.push_back(out.V(line.from) * std::conj(I_line));
    }

    out.residual = (in.Y * out.V - in.I).norm() / std::max(1.0, in.I.norm());
    return out;
}