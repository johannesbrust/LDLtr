%---------------------- EXAMPLE.m ----------------------------------------%
%
% Example run of LDLtr on problem without CUTEst installation
% For details on how the CUETEst experiments are done one can
% look into runUCset.m, bfgs_parms.m, bfgsTest.m, and runProb.m,
% for instance
%
%-------------------------------------------------------------------------%
% 07/12/23, J.B., preparation for release

addpath('../');


% Load a problem
[prob] = getCUSTomProblem();

% Run the solver, specified in bfgs_parms.m
[x,status,Stats] = bfgsTest( prob );