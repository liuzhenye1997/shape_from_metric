#define _CRT_SECURE_NO_WARNINGS
#define _SCL_SECURE_NO_WARNINGS

#include <mex.h>
#include <string>
#include <vector>

#include "splsolver.h"

const mxClassID mxPOINTER_CLASS = mxINDEX_CLASS;
using pointer_t = size_t;

//////////////////////////////////////////////////////////////////////////
/* memory management - redirect to MATLAB memory manager */
void* operator new(size_t size)
{
    void *ptr = mxMalloc(size);
    mexMakeMemoryPersistent(ptr);
    return ptr;
}
void* operator new[](size_t size)
{
    void *ptr = mxMalloc(size);
    mexMakeMemoryPersistent(ptr);
    return ptr;
}

void operator delete(void* ptr) { mxFree(ptr); }
void operator delete[](void* ptr) { mxFree(ptr); }


//////////////////////////////////////////////////////////////////////////
// replace mxGetNzmax with nnz, as mxGetNzmax has over-allocated elements (zeros)
size_t nnz(const mxArray* m) { return *(mxGetJc(m)+mxGetN(m)); }

void verifyMatIsSquareSparse(const mxArray* m, bool checkDoublePrecision=true) 
{
    if (!mxIsSparse(m)) mexErrMsgTxt("Invalid input: matrix should be sparse.");
    if (mxGetM(m) != mxGetN(m)) mexErrMsgTxt("Invalid input: matrix should be square matrix (nxn).");

    if (checkDoublePrecision && !mxIsDouble(m)) mexErrMsgTxt("Invalid input: matrix should be in double precision.");
}

Eigen::Map<const SpMat> extractSpMat(const mxArray *m, double* replaceNonzeros = nullptr)
{
    static_assert(sizeof(index_t)==8, "only for 64bit");
    return Eigen::Map<const SpMat>(mxGetM(m), mxGetN(m), nnz(m), (const index_t*)mxGetJc(m), (const index_t*)mxGetIr(m), replaceNonzeros?replaceNonzeros:mxGetPr(m));
}

enum spl_solve_stage{
    LST_SYMFACT = 1,
    LST_NUMFACT = 2,
    LST_SOLVE = 3,
    LST_NUMFACT_AND_SOLVE = 4,
    LST_ALL = 5,
    LST_CLOSE_HANDLE = 6
};

//////////////////////////////////////////////////////////////////////////
void mexFunction(int nlhs, mxArray *plhs[],	int nrhs, const mxArray*prhs[])
{
    if (nrhs < 2) 
        mexErrMsgTxt("Invalid input: not enough input, x = SPLsolver(mode, handle, A, Av, b);");

    using MapMatC = Eigen::Map<const Eigen::MatrixXd>;
    using MapMat = Eigen::Map<Eigen::MatrixXd>;

    const int mode = int( mxGetScalar(prhs[0]) );

    if (mode == LST_SYMFACT) {
        // SPLsolver(mode, A, solver)
        int solvertype = nrhs > 2 ? int(mxGetScalar(prhs[2])) : SLS_PARDISO_LU;
        SparseLinearSolver *solver = new SparseLinearSolver(solvertype);
        plhs[0] = mxCreateNumericMatrix(1, 1, mxPOINTER_CLASS, mxREAL);
        pointer_t *ptr_handle = (pointer_t*)mxGetData(plhs[0]);
        *ptr_handle = (pointer_t)solver;

        const mxArray *mat_A = prhs[1];
        verifyMatIsSquareSparse(mat_A, false);
        solver->analyzePattern( extractSpMat(mat_A) );// not sure if SPLsolver will access the nonzeros during symbolic factorization
        return;
    }

    if (mode == LST_ALL) {
        // SPLsolver(mode, A, b, solvername)
        if( nrhs<3 )
            mexErrMsgTxt("Invalid input: not enough input, SPLsolver(mode, A, b, solvername)");

        const mxArray *mat_A = prhs[1];
        const mxArray *mat_b = prhs[2];
        verifyMatIsSquareSparse(mat_A);

        const size_t m = mxGetM(mat_A);
        const size_t n = mxGetN(mat_b);

        if (!mxIsDouble(mat_b)) mexErrMsgTxt("Invalid input: vector b should be in double precision.");
        if (mxGetM(mat_b) != m) mexErrMsgTxt("Invalid input: the size of vector b does not match the coefficient matrix.");

        plhs[0] = mxCreateDoubleMatrix(m, n, mxREAL);

        int solvertype = nrhs > 3 ? int(mxGetScalar(prhs[3])) : SLS_PARDISO_LU;
        bool suc = SparseLinearSolver::solve(solvertype, extractSpMat(mat_A), mxGetPr(mat_b), mxGetPr(plhs[0]), n);
        if (!suc) fprintf(stderr, "SparseLinearSolver failed to fullsolve\n");
        if (nlhs > 1) plhs[1] = mxCreateLogicalScalar(suc);
        return;
    }

    const pointer_t *ptr_handle = (pointer_t*)mxGetData(prhs[1]);
    SparseLinearSolver *solver = (SparseLinearSolver*)(*ptr_handle);
    if (!solver->isMemOk()) mexErrMsgTxt("bad memory! is the handle deleted already?");

    switch (mode) {
    case LST_CLOSE_HANDLE:
        delete solver;
        break;
    case LST_NUMFACT:
    case LST_NUMFACT_AND_SOLVE: 
    case LST_SOLVE:
    if (mode != LST_SOLVE) {
        // factorization only: SPLsolver(mode, handle, A, nonzeros)
        if( nrhs<4+(mode==LST_NUMFACT_AND_SOLVE) )
            mexErrMsgTxt("Invalid input: not enough input, missing new coefficient matrix (nonzero structure matrix + nonzero array) for factorization");

        const mxArray *mat_A = prhs[2];
        const mxArray *mat_nonzeros = prhs[3];

        verifyMatIsSquareSparse(mat_A, false);

        const size_t m = solver->rows();
        if (m != mxGetN(mat_A)) mexErrMsgTxt("Invalid input: the size of coefficient matrix does not match with the original input.");

        if (!mxIsDouble(mat_nonzeros)) mexErrMsgTxt("Invalid input: nonzero vector should be in double precision.");

        if (mxGetNumberOfElements(mat_nonzeros) < nnz(mat_A)) mexErrMsgTxt("Invalid input: size of nonzero vector does not match the sparse matrix.");

        bool suc = solver->factorize( extractSpMat(mat_A, mxGetPr(mat_nonzeros)) );
        if (!suc) fprintf(stderr, "SparseLinearSolver failed to factorize\n");
    }

    if (mode == LST_NUMFACT) {
        plhs[0] = mxCreateLogicalScalar(solver->success());
    }
    else{
        // solve only: SPLsolver(mode, handle, b);
        // numeric factorization and solve:  SPLsolver(mode, handle, A, nonzeros, b)
        int iargb = 2 + (mode == LST_NUMFACT_AND_SOLVE) * 2;
        if( nrhs<iargb+1 )
            mexErrMsgTxt("Invalid input: not enough input, missing right hand side vector for solve");

        const mxArray *mat_b = prhs[iargb];
        if (!mxIsDouble(mat_b)) mexErrMsgTxt("Invalid input: vector b should be in double precision.");

        const size_t m = solver->rows();
        if (mxGetM(mat_b) != m) mexErrMsgTxt("Invalid input: the size of vector b does not match the coefficient matrix.");

        const size_t n = mxGetN(mat_b);
        plhs[0] = mxCreateDoubleMatrix(m, n, mxREAL);

        bool suc = solver->success() && solver->solve( mxGetPr(mat_b), mxGetPr(plhs[0]), n ) ;

        if (!suc) fprintf(stderr, "SparseLinearSolver failed to solve\n");

        if(nlhs>1) plhs[1] = mxCreateLogicalScalar(suc);
    }
    break;

    default:
        mexErrMsgTxt("unrecognized mode");
    }
}
