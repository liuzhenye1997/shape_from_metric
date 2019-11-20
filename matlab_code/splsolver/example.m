% example code for using the splsolver

A = sprand(10, 10, 0.2);
A = A*A' + speye(10)*1e-2;
b = rand(10, 2);

%% init pardiso solver & perform symbfact with the nonzero pattern of A
solver = splsolver(A, 'llt');   % for SPD matrix A, llt (cholesky) is the fastest solver

%% option 1: perform numeical factorization and solve separately
% solve can be performed multiple times, i.e. reusing the symbolic and
% numerical factorization, if the matrix does not change
solver.refactorize(A);
x1 = solver.solve(b);

%% option 2: combine numerical factorization & solve Ax = b
x2 = solver.refactor_solve(A, b);

%% if the matrix changes (but with the SAME sparsity pattern), then the numerical factorization need to be perfomed again.
% refactorization & solve Bx = b
B = A + spdiags( rand(10,1), 0, 10, 10 );
x3 = solver.refactor_solve(B, b);
