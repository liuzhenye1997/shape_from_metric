classdef splsolver < handle
    properties(SetAccess = protected)
        handle = -1;
        A = [];  % matrix, nonzero pattern
    end
    
    methods(Static)
        function i = solverid(solvername)
            switch lower(solvername)
                case {'llt', 'llt_pardiso'}
                    i = 0;
                case {'ldlt', 'ldlt_pardiso'}
                    i = 1;
                case {'lu', 'lu_pardiso'}
                    i = 2;
                case {'lu_umf', 'umf'}
                    i = 3;
                case {'cholmod', 'cholmod_auto'}
                    i = 4;
                case {'cholmod_llt_simplicial', 'cholmod_simplicial'}
                    i = 5;
                case {'cholmod_llt_supernodal', 'cholmod_supernodal'}
                    i = 6;
                case 'cholmod_ldlt'
                    i = 7;
                otherwise
                    error('unknown solver %s', solvername);
            end
        end
        
        function print(varargin)
            verbose = false;
            if verbose, fprintf(varargin{:}); end
        end
        
        function x = fullsolve(A, b, solvername)
            if nargin<3, solvername = 'lu'; end
            tic
            [x, f] = splsolver_imp(5, A, b, splsolver.solverid(solvername));
            splsolver.print('\nsplsolver full      = %fs', toc);
            assert(f, 'linear solver failed');
        end
    end

    methods
        function s = splsolver(A, solvername)
            if nargin<2, solvername = 'lu'; end
            s.A = A;
            tic
            s.handle = splsolver_imp(1, A, splsolver.solverid(solvername));
            splsolver.print('\nsplsolver    symf   = %fs', toc);
        end

        function refactorize(s, nonzeros)
            % TODO: checke nonzeros is a vector
            % TODO: return complex if input right side is complex
            tic
            f = splsolver_imp(2, s.handle, s.A, nonzeros);
            assert(f, 'linear solver failed with factorization');
            splsolver.print('\nsplsolver    numf   = %fs', toc);
        end
        
        function x = solve(s, b)
            tic
            x = splsolver_imp(3, s.handle, b);
            splsolver.print('\nsplsolver    solve  = %fs', toc);
        end
        
        function x = mldivide(s, b)
            tic
            x = splsolver_imp(3, s.handle, b);
            fprintf('\nsplsolver    solve  = %fs', toc);
        end
        
        function [x, f] = refactor_solve(s, nonzeros, b)
            tic
            [x, f] = splsolver_imp(4, s.handle, s.A, nonzeros, b);
            splsolver.print('\nsplsolver    f&s    = %fs\n', toc);
            assert(f, 'linear solver failed');
        end
        
        function delete(s)
            if s.handle ~= -1
                tic
                splsolver_imp(6, s.handle);
                splsolver.print('\nsplsolver    clean up   = %fs\n', toc);
            end
        end
    end
end
