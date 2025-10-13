/* =============================================================================
 * KluSolverEigen：对 SuiteSparse/KLU 的极薄封装（与 Eigen 对接）
 * -----------------------------------------------------------------------------
 * 用法建议：
 *   - 首次/结构变化：调用 analyzePattern(A)（要求 A 已压缩：isCompressed() == true）
 *   - 结构不变，仅数值更新：反复调用 factorize(A)
 *   - 求解：solve(b)
 *
 * 注意事项：
 *   - A 必须为方阵，且 nonZeros 在结构阶段固定不变
 *   - analyzePattern() 内部复制 CSC 结构，factorize() 仅刷新 Ax
 *   - 该类在析构时会释放 KLU 的 symbolic/numeric 资源
 * ============================================================================= */

#pragma once
#include <Eigen/Sparse>
#include <stdexcept>
#include <suitesparse/klu.h>
#include <vector>

class KluSolverEigen {
public:
    KluSolverEigen() = default;
    ~KluSolverEigen() { release(); }

    // 结构分析（模式分析），仅当稀疏结构改变时需要重新调用
    void analyzePattern(const Eigen::SparseMatrix<double>& A)
    {
        release(); // 释放已存在的对象
        if (!A.isCompressed()) {
            // Eigen 构建后通常已压缩；如未压缩，用户需 setFromTriplets 后调用 makeCompressed()
            throw std::runtime_error("KluSolverEigen: matrix must be compressed.");
        }
        n_ = A.rows();
        if (A.cols() != n_)
            throw std::runtime_error("KluSolverEigen: matrix must be square.");
        nnz_ = A.nonZeros();

        // 拷贝 CSC 结构
        Ap_.assign(A.outerIndexPtr(), A.outerIndexPtr() + (n_ + 1));
        Ai_.assign(A.innerIndexPtr(), A.innerIndexPtr() + nnz_);
        Ax_.assign(A.valuePtr(), A.valuePtr() + nnz_);

        // KLU symbolic
        klu_common c;
        klu_defaults(&c);
        S_ = klu_analyze(n_, Ap_.data(), Ai_.data(), &c);
        if (!S_)
            throw std::runtime_error("KluSolverEigen: klu_analyze failed.");
        analyzed_ = true;
    }

    // 数值分解（结构不变时反复调用即可）
    void factorize(const Eigen::SparseMatrix<double>& A)
    {
        if (!analyzed_)
            analyzePattern(A); // 保险：若外部忘记调用，自动分析一次
        if (!A.isCompressed())
            throw std::runtime_error("KluSolverEigen: matrix must be compressed.");
        if (A.rows() != n_ || A.cols() != n_)
            throw std::runtime_error("KluSolverEigen: size changed; call analyzePattern().");
        if (A.nonZeros() != nnz_)
            throw std::runtime_error("KluSolverEigen: nonzero count changed; call analyzePattern().");

        // 更新数值
        std::copy(A.valuePtr(), A.valuePtr() + nnz_, Ax_.begin());

        klu_common c;
        klu_defaults(&c);
        // 若已有 N_，则使用 refactor 更高效；否则第一次用 factor
        if (N_) {
            if (!klu_refactor(Ap_.data(), Ai_.data(), Ax_.data(), S_, N_, &c))
                throw std::runtime_error("KluSolverEigen: klu_refactor failed.");
        } else {
            N_ = klu_factor(Ap_.data(), Ai_.data(), Ax_.data(), S_, &c);
            if (!N_)
                throw std::runtime_error("KluSolverEigen: klu_factor failed.");
        }
    }

    /// 解线性方程组 G·x=b；要求已完成 factorize()，RHS 维度一致
    Eigen::VectorXd solve(const Eigen::VectorXd& b) const
    {
        if (!N_ || !S_)
            throw std::runtime_error("KluSolverEigen: factorize() was not called.");
        if (b.size() != n_)
            throw std::runtime_error("KluSolverEigen: RHS size mismatch.");
        Eigen::VectorXd x = b; // KLU 就地覆盖
        klu_common c;
        klu_defaults(&c);
        int nrhs = 1;
        if (!klu_solve(S_, N_, n_, nrhs, x.data(), &c))
            throw std::runtime_error("KluSolverEigen: klu_solve failed.");
        return x;
    }

private:
    void release()
    {
        klu_common c;
        klu_defaults(&c);
        if (N_) {
            klu_free_numeric(&N_, &c);
            N_ = nullptr;
        }
        if (S_) {
            klu_free_symbolic(&S_, &c);
            S_ = nullptr;
        }
        analyzed_ = false;
        n_ = 0;
        nnz_ = 0;
        Ap_.clear();
        Ai_.clear();
        Ax_.clear();
    }

    // KLU 内部对象
    klu_symbolic* S_ = nullptr;
    klu_numeric* N_ = nullptr;

    // 维度与数据
    int n_ = 0;
    int nnz_ = 0;
    std::vector<int> Ap_; // 列指针
    std::vector<int> Ai_; // 行索引
    std::vector<double> Ax_; // 数值
    bool analyzed_ = false;
};
