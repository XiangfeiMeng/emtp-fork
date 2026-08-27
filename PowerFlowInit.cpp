#include "PowerFlowInit.h"

#include <cmath>
#include <complex>
#include <stdexcept>

namespace PFInit {
namespace {
using Complex = std::complex<double>;
constexpr double pi = 3.14159265358979323846;
constexpr double phaseShift = 2.0 * pi / 3.0;

Complex phaseFactor(int p)
{
    const double a = p == 0 ? 0.0 : (p == 1 ? -phaseShift : phaseShift);
    return std::polar(1.0, a);
}

double instant(Complex x, double w, double t)
{
    return std::sqrt(2.0) * std::real(x * std::polar(1.0, w * t));
}

Eigen::Vector3cd phases(Complex Va)
{
    return { Va, Va * phaseFactor(1), Va * phaseFactor(2) };
}

Complex busV(const Result& flow, int bus)
{
    if (bus < 0)
        return 0.0;
    if (bus >= flow.V.size())
        throw std::invalid_argument("initial-state bus index is invalid");
    return flow.V(bus);
}

void checkRange(int first, int count, int size)
{
    if (first < 0 || count < 0 || first + count > size)
        throw std::invalid_argument("initial-state internal node range is invalid");
}
}

Data makeData(double freq)
{
    Data d;
    d.freq = freq;
    d.nodeCount = 63;
    const double v[] = { 17.16, 18.45, 14.145 };
    const double a[] = { 0.0, 9.28001, 4.66475 };
    const double Ls[] = { 0.000456, 0.000189, 0.000128 };
    for (int k = 0; k < 3; ++k) d.sources.push_back({ k, 27 + 3 * k, v[k], a[k], 0.001, Ls[k] });
    const double n[] = { 16.5 / 230.0, 18.0 / 230.0, 13.8 / 230.0 };
    const double r1[] = { 2.7225e-6, 3.24e-6, 1.9044e-6 };
    const double l2[] = { 0.080825, 0.087701, 0.082228 };
    const double rm[] = { 1361.2, 1620.0, 952.2 };
    const double lm[] = { 3.6108, 4.2972, 2.5258 };
    for (int k = 0; k < 3; ++k) d.xfmrs.push_back({ k, 3 + k, 36 + 9 * k, n[k], r1[k], 0.000529, l2[k], rm[k], lm[k] });
    const double LL[] = { 2.80643216318709, 4.00918880455298, 4.67738693864515 };
    for (int k = 0; k < 3; ++k) d.loads.push_back({ 6 + k, LL[k] });
    Eigen::Matrix3d R, L;
    R << 13.5293831769454,12.2635478388026,12.2635478388026,12.2635478388026,13.5293831769454,12.2635478388026,12.2635478388026,12.2635478388026,13.5293831769454;
    L << 0.198642487743882,0.105535006650437,0.105535006650437,0.105535006650437,0.198642487743882,0.105535006650437,0.105535006650437,0.105535006650437,0.198642487743882;
    const int lines[][2] = {{3,6},{3,8},{4,6},{4,7},{5,7},{5,8}};
    for (const auto& line : lines) d.lines.push_back({ line[0], line[1], R, L, Eigen::Vector3d::Constant(0.637898940812743e-6), 2.99136421769046e-6 });
    return d;
}
State init(const Result& flow, const Data& data, double t, double dt)
{
    if (data.freq <= 0.0 || dt <= 0.0)
        throw std::invalid_argument("initial-state frequency and time step must be positive");
    if (data.nodeCount < 3 * flow.V.size())
        throw std::invalid_argument("initial-state node count is too small");

    const double w = 2.0 * pi * data.freq;
    State state;
    state.V = Eigen::VectorXd::Zero(data.nodeCount);

    for (int bus = 0; bus < flow.V.size(); ++bus) {
        const auto Vabc = phases(flow.V(bus));
        for (int p = 0; p < 3; ++p)
            state.V(3 * bus + p) = instant(Vabc(p), w, t);
    }

    state.source.reserve(3 * data.sources.size());
    for (const auto& src : data.sources) {
        if (src.L <= 0.0)
            throw std::invalid_argument("source inductance must be positive");
        checkRange(src.internal, 3, data.nodeCount);
        const Complex Z(src.R, w * src.L);
        const Complex Ea = std::polar(src.V_kV / std::sqrt(3.0),
            (src.V_A - 90.0) * pi / 180.0);
        const double G = dt / (2.0 * src.L);

        for (int p = 0; p < 3; ++p) {
            const Complex k = phaseFactor(p);
            const Complex E = Ea * k;
            const Complex V = busV(flow, src.bus) * k;
            const Complex I = (E - V) / Z;
            const Complex Vi = E - src.R * I;
            const double v = instant(V, w, t);
            const double vi = instant(Vi, w, t);

            state.V(src.internal + p) = vi;
            state.source.push_back(-instant(I, w, t) - G * (v - vi));
        }
    }

    for (const auto& xfmr : data.xfmrs) {
        if (xfmr.ratio == 0.0 || xfmr.R1 == 0.0 || xfmr.L2 <= 0.0
            || xfmr.Rm == 0.0 || xfmr.Lm <= 0.0)
            throw std::invalid_argument("transformer initial-state parameters are invalid");
        checkRange(xfmr.internal, 9, data.nodeCount);

        const Complex Yp = 1.0 / Complex(xfmr.R1, 0.0);
        const Complex Ys = 1.0 / Complex(xfmr.R2, w * xfmr.L2);
        const Complex Ym = 1.0 / xfmr.Rm + 1.0 / Complex(0.0, w * xfmr.Lm);
        const Complex V1 = busV(flow, xfmr.bus1);
        const Complex V2 = busV(flow, xfmr.bus2);
        const Complex Vi1 = (Yp * V1 + Ys * V2 / xfmr.ratio)
            / (Yp + Ym + Ys / (xfmr.ratio * xfmr.ratio));
        const Complex Vi2 = Vi1 / xfmr.ratio;
        const Complex Ic = Ys * (Vi2 - V2) / xfmr.ratio;

        for (int p = 0; p < 3; ++p) {
            const Complex k = phaseFactor(p);
            const int base = xfmr.internal + 3 * p;
            state.V(base) = instant(Ic * k, w, t);
            state.V(base + 1) = instant(Vi1 * k, w, t);
            state.V(base + 2) = instant(Vi2 * k, w, t);
        }
    }

    state.load.reserve(data.loads.size());
    for (const auto& load : data.loads) {
        if (load.L <= 0.0)
            throw std::invalid_argument("load inductance must be positive");
        const auto V = phases(busV(flow, load.bus));
        const double G = dt / (2.0 * load.L);
        Eigen::Vector3d hist;
        for (int p = 0; p < 3; ++p) {
            const Complex I = V(p) / Complex(0.0, w * load.L);
            hist(p) = instant(I, w, t) - G * instant(V(p), w, t);
        }
        state.load.push_back(hist);
    }

    state.line.reserve(data.lines.size());
    for (const auto& line : data.lines) {
        const auto Vi = phases(busV(flow, line.from));
        const auto Vj = phases(busV(flow, line.to));
        const Eigen::Matrix3cd Z = line.R.cast<Complex>()
            + Complex(0.0, w) * line.L.cast<Complex>();
        const Eigen::Vector3cd I = Z.fullPivLu().solve(Vi - Vj);
        const Eigen::Matrix3d G = (line.R + (2.0 / dt) * line.L).inverse();

        Eigen::Vector3d vi, vj, current, ci, cj;
        for (int p = 0; p < 3; ++p) {
            vi(p) = instant(Vi(p), w, t);
            vj(p) = instant(Vj(p), w, t);
            current(p) = instant(I(p), w, t);
            ci(p) = instant(Complex(0.0, w * line.Cp(p)) * Vi(p), w, t);
            cj(p) = instant(Complex(0.0, w * line.Cp(p)) * Vj(p), w, t);
        }

        LineState hist;
        hist.series = current - G * (vi - vj);
        hist.phase_i = ci - (2.0 / dt) * line.Cp.cwiseProduct(vi);
        hist.phase_j = cj - (2.0 / dt) * line.Cp.cwiseProduct(vj);

        const Complex V0i = Vi.sum() / 3.0;
        const Complex V0j = Vj.sum() / 3.0;
        const double v0i = instant(V0i, w, t);
        const double v0j = instant(V0j, w, t);
        hist.ground_i = instant(Complex(0.0, w * line.Cg) * V0i, w, t)
            - (2.0 * line.Cg / dt) * v0i;
        hist.ground_j = instant(Complex(0.0, w * line.Cg) * V0j, w, t)
            - (2.0 * line.Cg / dt) * v0j;
        state.line.push_back(hist);
    }

    state.inductor.reserve(data.inductors.size());
    for (const auto& item : data.inductors) {
        if (item.L <= 0.0)
            throw std::invalid_argument("inductance must be positive");
        const Complex V = busV(flow, item.from) - busV(flow, item.to);
        const Complex I = V / Complex(0.0, w * item.L);
        state.inductor.push_back(instant(I, w, t)
            - dt / (2.0 * item.L) * instant(V, w, t));
    }

    state.capacitor.reserve(data.capacitors.size());
    for (const auto& item : data.capacitors) {
        if (item.C <= 0.0)
            throw std::invalid_argument("capacitance must be positive");
        const Complex V = busV(flow, item.from) - busV(flow, item.to);
        const Complex I = Complex(0.0, w * item.C) * V;
        state.capacitor.push_back(instant(I, w, t)
            - 2.0 * item.C / dt * instant(V, w, t));
    }

    return state;
}

}