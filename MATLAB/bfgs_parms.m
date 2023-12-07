function [parms] = bfgs_parms()
% ------------------------------------------------------------------------
% Authors : Philip E. Gill
% Date    : November 16, 2019
% Purpose : Sets the control parameters
% ------------------------------------------------------------------------
% Updates
%   01/24/23, J.B., Addition of trust-region parameters 
%   11/17/23, J.B., Preparation for release
%   11/29/23, J.B., Further preparations
% ------------------------------------------------------------------------
  parms.solver        = 'LDLtr'  ;  % bfgsR, bfgsTR_MS
                                     
  parms.maxIterations = 1500      ;  % maximum iterates allowed
  parms.maxIterations = 3000      ;  % maximum iterates allowed
  parms.maxIterations = 6000      ;  % maximum iterates allowed

  parms.tiny          = 1.0e-12   ;  % tiny number
  parms.tolStat       = 1.0e-4    ;  % stationary tolerance

  parms.wolfeTols     = zeros(1,2);
  parms.wolfeTols(1)  = 0.9       ;  % Wolfe gradient tolerance for quasi-Newton
  parms.wolfeTols(2)  = 1.0e-4    ;  % Wolfe function tolerance
  parms.jfmax         = 20        ;  % max functions per line search
  parms.gammaC        = 0.5       ;  % contraction factor in backtracking
  parms.condMax       = 1.0e+16   ;  % maximum condition estimate
  parms.dxMax         = 100       ;  % maximum  change in x
  parms.dxInf         = 1.0e+10   ;  % infinite change in x

% parms.updateTol     = 1.0e-8    ;  % Tolerance for rejecting BFGS update
  parms.updateTol     = 0.99      ;  % Tolerance for rejecting BFGS update
  parms.ResetFreq     = NaN       ;  % Reset frequency
  parms.lines         = 15        ;  % lines between printing header
  parms.printSol      = 1         ;  % to print solution to screen
  parms.infinity      = 1.0e+15   ;  % user-defined +infinity
  parms.Assert        = 0         ;  % (0/1) check consistency
  parms.trials        = 1         ;  % Repeat each problem # times; metrics averaged
  parms.show_warnings = false     ;  % (0/1) print warnings
  parms.count_near    = true      ;  % count near optimal solutions
  parms.print_run_det = true      ;  % option to print run details
  parms.write_file    = true      ;  % write final run statistics to file
  parms.out_file      = false     ;  % write one line per iteration to file
  parms.searchType    = 'Armijo'  ;  % Armijo           line search
  parms.searchType    = 'Wolfe'   ;  % Wolfe            line search
  %parms.searchType    = 'MT'      ;  % More-Thuente     line search

  parms.unBoundedf    =-1.0e+9    ;  % definition of an unbounded objective value
  parms.debugItn      = NaN       ;  % iteration for setting dbstop
  parms.debugItn      = 125       ;  % iteration for setting dbstop
  
  % Trust-region parameters (proposed strategy)
  parms.c1            = 1e-5      ; % step acceptance
  parms.c2            = 0.75      ; % function improvement: if c2 < rho
  parms.c3            = 0.8       ; % fraction to radius: if nsk <= c3*Delta
  parms.c4            = 2         ; % radius increase
  parms.c5            = 0.1       ; % lower range for const. Delta: c5 <= rho <= c6
  parms.c6            = 0.75      ; % upper range for const. Delta
  parms.c7            = 0.5       ; % radius decrease

  parms.trTol         = 1e-7      ; % tolerance for the TR subproblem
  parms.trMax         = 30        ; % maximum iterations for the TR subproblem

  % CG subproblem tolerances (related to Alg. 3)
  parms.tolICG        = 1e-7      ; % Relative stopping for the inv. CG iteration
  parms.maxitICG      = 35        ; % Maximum iterations for CG
  parms.shrink        = 0.5       ; % Contraction for the CG iteration
  parms.maxitBack     = 3         ; % Maximum iterations for the backtrack
  
  % Branching to More-Sorensen
  parms.nmax          = 100       ; % Problem dimension to use MS for
  
  % TR-Medium paramters (More-Sorensen trust-region subproblem)
  parms.TRrelTol      = 1e-6      ;  % 1e-6 trust-region relative tolerance.
  parms.TRabsTol      = 1e-3      ;  % trust-region absolute tolerance.
  parms.TRlog         = 0         ;  % trust-region log off/on
  parms.TRitnMax      = 30        ;  % trust-region iteration limit

  parms.probsize      = 50000     ; % Only run problems below n=probsize
  
end