function [prob] = getCUTEstProblem(pname)
%        [prob] = getCUTEstProblem(pname)
% Defines the data for a nonlinearly constrained problem of the form:
%
%         minimize   f(x)
%         subject to  cL <= c(x) <= cU
%                     xL <=   x  <= xU
%
% where cL and cU are vectors of lower and upper bounds on the constraint
% vector-valued function c(x), and xL and xU are vectors of lower and upper
% bounds on the primal variables x.
%
% The vector c(x) can be empty.
% cL and xL can be -inf, cU and xU can be +inf.
%
% -------------------------------------------------------------------------
% INPUT : pname   the name of the CUTEst problem to be solved. e.g., 'HS11'
%
% OUTPUT: prob    problem to be solved, a struct.
% -------------------------------------------------------------------------

%==========================================================================
% Authors : Philip E. Gill, Vyacheslav Kungurtsev and Daniel P. Robinson
% Date    : October, 2016
% Purpose : Solve general optimization problems.

% Intended for pedagogical use only. It's only a toy!
% Current : January 14, 2019.
% 
% Updates:
% 02/12/21, J.B., Addition of 'cutest2matlab_osx' function
%==========================================================================

format long  E
format compact

% Decode the cutest problem given as input parameter pname.

%path1 = getenv('PATH')
%path1 = [path1 ':/Users/peg/cutest/CUTEst/bin/']
%setenv('PATH', path1)

%a        = sprintf('%s','cutest2matlab');
if isunix == 1 && ismac == 0
    a = sprintf('%s','runcutest -p matlab -D ');
    c = sprintf('%s %s',a,pname);

    a = sprintf('%s','cutest2matlab $MASTSIF/');
    %a = sprintf('%s','cutest2matlab_osx $MASTSIFA/');
    c = sprintf('%s%s',a,strtrim(pname));

else
    a = sprintf('%s','cutest2matlab_osx $MASTSIF/');
    %a = sprintf('%s','cutest2matlab_osx $MASTSIFA/');
    c = sprintf('%s%s',a,strtrim(pname));
end

unix(c);
prob      = cutest_setup;
prob.name = pname;

% Function handles for cutest problem.

prob.cons = @(x)cutest_cons( x );
prob.obj  = @(x)cutest_obj( x );
if  prob.m > 0,
  prob.hess = @(x,y)cutest_hess(x, y);
else
  fprintf('Warning : Problem has no general constraints.\n\n')
  prob.hess = @(x)cutest_hess(x);
end

% Collect and print statistics.

infinity = 1.0e+15 ;  % definition of +infinity

bl       = prob.bl;
bu       = prob.bu;
cl       = prob.cl;
cu       = prob.cu;

equality = [ prob.equatn ];

% Variables and slacks with finite lower and upper bounds.

nFixed = length(find( bl ==  bu));
nFree  = length(find((bl < -infinity) & (bu >  infinity)));

nxLow  = length(find((bl > -infinity) & (bu >  infinity)));
ncLow  = length(find((cl > -infinity) & (cu >  infinity)));

nxUpp  = length(find((bu <  infinity) & (bl < -infinity)));
ncUpp  = length(find((cu <  infinity) & (cl < -infinity)));

nxRnge = length(find((bl > -infinity) & (bu <  infinity)));
ncRnge = length(find((cl > -infinity) & (cu <  infinity) & ~equality));

nEq    = length(find(equality));

% Print information about the problem.

fprintf('\nProblem');
fprintf('\n=======\n');
fprintf(' Name                            : %14s\n',prob.name);
fprintf(' Variables                       : %14d\n',prob.n);
fprintf(' General Constraints             : %14d\n',prob.m);

if nEq    > 0, fprintf(' Equality constraints            : %14d\n',nEq   ); end
if nFixed > 0, fprintf(' Fixed variables                 : %14d\n',nFixed); end
if nFree  > 0, fprintf(' Free  variables                 : %14d\n',nFree ); end
if nxLow  > 0, fprintf(' Variables with lower bounds     : %14d\n',nxLow ); end
if nxUpp  > 0, fprintf(' Variables with upper bounds     : %14d\n',nxUpp ); end
if nxRnge > 0, fprintf(' Variables with range constraints: %14d\n',nxRnge); end
if ncLow  > 0, fprintf(' Slacks    with lower bounds     : %14d\n',ncLow ); end
if ncUpp  > 0, fprintf(' Slacks    with upper bounds     : %14d\n',ncUpp ); end
if ncRnge > 0, fprintf(' Slacks    with range constraints: %14d\n',ncRnge); end
fprintf('\n');
