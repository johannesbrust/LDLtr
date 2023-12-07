function [x,status,Stats] =  runProb( varargin )
%        [x,status,Stats] =  runProb( varargin )
% If 0 arguments are provided it will run the custom problem defined
% in getCUSTomProblem.m using whatever solver is already set in the bfgs_parms file:
% >> runProb();
%
% If 1 argument is provided it is assumed to be a cutest problem name, and the
% given cutest problem will be run using whatever solver is already set in the
% bfgs_parms file:
% >> runProb(‘AKIVA’);
%
% if 2 arguments are provided it assumes the latter is the desired solver name
% and will set the parms.solver value in the parms file and run cutest.
% >> runProb(‘AKIVA’, ‘bfgsZ’);
%
% Warning: calling with one argument that is not a problem name, e.g.,
% runProb(‘bfgsM’) will not work correctly.
%--------------------------------------------------------------------------
% Updates
% 12/02/21, J.B., deleting "*.d,*.o,*.dylib,*.f,mcutest.*" files on MacOS
% 04/03/22, J.B., Modification to run all problems that have <25000
% variables

if nargin == 2
   % Read in the parms file and set the solver variable to the provided solver
   solver = varargin{2};
   fname = 'bfgs_parms.m';
   ptxt  = fileread(fname);
   ptxt  = regexprep(ptxt, 'bfgs[CFHILMNRSYZ]+', solver);
   pid   = fopen(fname, 'w');
   fwrite(pid, ptxt);
   fclose(pid);
end

if nargin > 0 % Solver and problem provided, run solver on cutest problem
   prob_name = varargin{1}
   [prob] = getCUTEstProblem( prob_name );
   
   [x,status,Stats] = bfgsTest(prob);
   
   cutest_terminate
else
   [prob] = getCUSTomProblem();
   [x,status,Stats] = bfgsTest(prob);
end
if ismac==1
    delete( '*.d',...
            '*.o',...
            '*.dylib',...
            '*.f',...
            'mcutest.*');
end