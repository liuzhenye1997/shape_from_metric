#include <vector>
#include <variant>

#include <Eigen/Dense>
#include <Eigen/SparseCore>
#define EIGEN_USE_MKL_ALL
#define EIGEN_CHOLMOD_SUPPORT
#include <Eigen/CholmodSupport>
#include <Eigen/UmfPackSupport>
#include <Eigen/PardisoSupport>


using index_t = int64_t;
//using index_t = int;
using SpMat = Eigen::SparseMatrix<double, 0, index_t>;

enum SPLINEAR_SOLVER_DIRECT { SLS_PARDISO_LLT = 0, SLS_PARDISO_LDLT, SLS_PARDISO_LU, SLS_UMF_LU,
                              SLS_CHOLMOD, SLS_CHOLMOD_LLT_SIMPILICIAL, SLS_CHOLMOD_LLT_SUPER, SLS_CHOLMOD_LDLT, SLS_NUM };

struct SparseLinearSolver
{
private:
    const int SLS_INTEGRITY_CHECK = 0x12345678;
    int memIntFlag = SLS_INTEGRITY_CHECK;

public:
    
    //using solver_t = Eigen::CholmodSupernodalLLT<SpMat>; 
    using solver_t = std::variant<Eigen::PardisoLLT<SpMat>, Eigen::PardisoLDLT<SpMat>, Eigen::PardisoLU<SpMat>, Eigen::UmfPackLU<SpMat>, Eigen::CholmodDecomposition<SpMat>>;
    solver_t imp;

    SparseLinearSolver(int s = 0) { emplaceSolver(imp, s); }

    static void emplaceSolver(solver_t &solver, int s)
    {
        const Eigen::CholmodMode cholmodMode[] = { Eigen::CholmodAuto, Eigen::CholmodSimplicialLLt, Eigen::CholmodSupernodalLLt, Eigen::CholmodLDLt };
        switch (s) {
        case SLS_PARDISO_LLT: solver.emplace<SLS_PARDISO_LLT>(); break;  //imp.emplace<Eigen::PardisoLLT<SpMat> >();
        case SLS_PARDISO_LDLT: solver.emplace<SLS_PARDISO_LDLT>(); break;
        case SLS_PARDISO_LU: solver.emplace<SLS_PARDISO_LU>(); break;
        case SLS_UMF_LU: 
            solver.emplace<SLS_UMF_LU>(); break;
        case SLS_CHOLMOD: 
        case SLS_CHOLMOD_LLT_SIMPILICIAL: 
        case SLS_CHOLMOD_LLT_SUPER: 
        case SLS_CHOLMOD_LDLT: 
            solver.emplace<SLS_CHOLMOD>(); break;
            std::get<SLS_CHOLMOD>(solver).setMode(cholmodMode[s-SLS_CHOLMOD]); break;
        default:
            fprintf(stderr, "unexpected linear solver %d!", s);
        }
    }

    template<class mat>
    static bool solve(int s, const mat& M, const double* b, double *x, size_t n) {
        const size_t m = M.rows();
        solver_t solver;
        emplaceSolver(solver, s);
        return std::visit([&](auto& s) {
            s.compute(M);
            if (s.info() != Eigen::Success) return false;
            s._solve_impl(Eigen::Map<const Eigen::MatrixXd>(b, m, n), Eigen::Map<Eigen::MatrixXd>(x, m, n));
            return s.info() == Eigen::Success;
        }, solver);
    }

    inline bool solve(const double* b, double *x, size_t n) {
        return std::visit([&](auto& solver) {
            const size_t m = solver.cols();
            solver._solve_impl(Eigen::Map<const Eigen::MatrixXd>(b, m, n), Eigen::Map<Eigen::MatrixXd>(x, m, n));
            return solver.info() == Eigen::Success;
        }, imp);
    }

    // use template matrix type instead of SpMat, otherwise a Map<SpMat> will be implicitly converted to a SpMat, which will be destructed later and causing problem
    // not sure if the solver will access the nonzeros during symbolic factorization
    // inline void analyzePattern(const SpMat& M) 
    template<class mat>
    inline void analyzePattern(const mat& M) { std::visit([&](auto& solver) { solver.analyzePattern(M);}, imp); }

    template<class mat>
    inline void compute(const mat& M) { std::visit([&](auto& solver) { solver.compute(M);}, imp); }

    template<class mat>
    inline bool factorize(const mat& M) {
        return std::visit([&](auto& solver) {
            solver.factorize(M);
            return solver.info() == Eigen::Success;
        }, imp);
    }

    inline size_t rows() { return std::visit([&](const auto& solver) { return solver.rows();}, imp); }

    inline bool success() { return std::visit([&](const auto& solver) { return solver.info()==Eigen::Success;}, imp);}

    ~SparseLinearSolver() { memIntFlag = 0; /*cancel the integrity check*/ }
    bool isMemOk() const { return memIntFlag == SLS_INTEGRITY_CHECK; }
};

