function [x,itn,nf,skipped,status] = bfgsM(prob,outfile)
%bfgsM   [x,itn,nf,skipped,status] = bfgsM(prob,outfile) finds
%        a local minimizer of a scalar-valued multivariate function
%        using a variant of the BFGS quasi-Newton method with a
%        Wolfe line search.
%
%        prob is a structure that provides the details of the
%        CUTEst problem being solved.
%
%        On successful termination, x is the first point found at
%        which the norm of the gradient is less than sqrt(eps).
%
%        Input arguments are set in the file bfgs_parms.m.
%
%        The meaning of each line of output is as follows
%
%        Itn          The current iteration number.
%
%        Nf           The cumulative number of times that the
%                     objective function f(x) has been evaluated.
%
%        Step         The steplength taken along the search direction.
%
%        f(x)         The current value of the objective function f(x).
%
%        Norm g       The current value of the gradient of f(x). In
%                     the neighborhood of a local minimizer, Norm g
%                     will be small.
%
%        cond H       condition number of the approximate Hessian
%
%        This version dated 29 Feb 2020.
%        Philip E. Gill, University of California, San Diego.

format compact

% -------------------------------
% Assign local control parameters
% -------------------------------
parms         = bfgs_parms();

solver        = parms.solver       ;  % Solver
maxIterations = parms.maxIterations;  % maximum iterates allowed
tolStny       = parms.tolStat      ;  % stationarity       tolerance
wolfeTols     = parms.wolfeTols    ;  % Wolfe function tolerance
dxMax         = parms.dxMax        ;  % maximum  change in x
dxInf         = parms.dxInf        ;  % infinite change in x
updateTol     = parms.updateTol    ;  % Tolerance for rejecting BFGS update
lines         = parms.lines        ;  % lines between printing header
ResetFreq     = parms.ResetFreq    ;  % Reset frequency
infinity      = parms.infinity     ;  % user-defined + infinity
unBoundedf    = parms.unBoundedf   ;  % unbounded objective value
debugItn      = parms.debugItn     ;  % iteration for setting dbstop
searchType    = parms.searchType   ;  % line search type
count_near    = parms.count_near   ;  % Option to count near-optimal results
print_run_det = parms.print_run_det;  % Option to print iteration details
% -----------------------
% End: control parameters
% -----------------------
if print_run_det
   str    = [ ' Solver' ];
   fprintf(        '%s%27s%15s\n', str, ':', solver);
   fprintf(outfile,'%s%27s%15s\n', str, ':', solver);

   str    = [ ' Line search' ];
   fprintf(        '%s%22s%15s\n', str, ':', searchType);
   fprintf(outfile,'%s%22s%15s\n', str, ':', searchType);
end
% Initialize the counters

itn      = 0;
nf       = 1;                 % number of function evaluations
skipped  = 0;

if isnan(ResetFreq), ResetFreq = maxIterations + 1; end
if isnan(debugItn ), debugItn  = maxIterations + 1; end

QN       =  0;                SD     = 1;

stepMax  = 10^5;
MjjMax   = 1.0d+4;            MjjMin = 1.0d-2;

x        = prob.x;
n        = prob.n;
[f,g]    = prob.obj(x);

normg    = norm(g);
normgMax = 1 + normg;
fMax     = 1 + abs(f);

Mjj    = 1/max(1,normg);      Mjj   =  min(max(Mjj,MjjMin),MjjMax);
M      = diag(Mjj*ones(1,n));

Info   = '     ';
% header = '  Itn    Nf    Step     jf    Objective      Norm g   cond H';
header = '  Itn    Nf    Step     jf    Objective      Norm g';

step   =  0;                  iExit = 1; % First line of output
dType  = SD;

Optimal  = false;   NearOptimal = false;   Unbounded = false;
ItnLimit = false;   BadSearch   = false;

%--------------------------------------------------------------------------
% Main loop.
%--------------------------------------------------------------------------
while (true)

  if      normg    <= tolStny
    Optimal     = true;
  elseif  f        <= unBoundedf
    Unbounded   = true;
  elseif count_near && BadSearch && (normg  <= eps^(2/3)*normgMax ||  ...
                                     abs(f) <= eps^(2/3)*fMax)
    NearOptimal = true;
  end

  ItnLimit  = itn > maxIterations - 1;

  if BadSearch && dType == QN
    Mjj     = 1/max(1,normg);    Mjj   =  min(max(Mjj,MjjMin),MjjMax);
    M       = diag(Mjj*ones(1,n));
    dType   = SD;              Info(5:5) = 'R';
    BadSearch = false;
  end

  %------------------------------------------------------------------------
  % Print iteration details.
  %------------------------------------------------------------------------
  if print_run_det
     nLine    = rem( itn, lines );
     if nLine == 0
       fprintf('\n%s\n', header);  fprintf(outfile, '\n%s\n', header);
     end

     if itn == 0
       cstep = '    --    ';
       cjf    =    '-- ';
     else
       cstep = sprintf ( ' %9.2e', step );
       cjf   = sprintf ( '%2g '  ,   jf );
     end
     if iExit ~= 1
       cExit = sprintf ( '(%1g)', iExit );
     else
       cExit = '   ';
     end

     str1  = sprintf ( '%5g %5g',        itn,      nf );
     str2  = sprintf ( '%s%s%s',               cstep,  cExit, cjf );
     str3  = sprintf ( '%14.7e %9.2e', f, norm(g) );
     % str3  = sprintf ( '%14.7e %9.2e %8.1e', f, norm(g), cond(M) );
     str4  = Info;
     str   = [ str1 str2 str3 str4 ];
     fprintf('%s\n', str);    fprintf(outfile, '%s\n', str);
     Info  = '     ';
  end

  if  Optimal || Unbounded || NearOptimal || ItnLimit || BadSearch
    break;
  end

  p         = -M*g;
  gtp       =  g'*p;     normp = norm(p);

  stepMax   = dxInf/(1+normp);
  stepLimit = dxMax/(1+normp);
  step      = min([1 stepLimit stepMax]);
  if step  == stepLimit, Info(5:5) = 'l'; end

  if     strcmp(searchType,'Wolfe' )
    [step,xnew,fnew,gnew,jf,iExit] ...
        = wolfeLSfg(prob,wolfeTols,step,stepMax,f,g,p,x);
  elseif strcmp(searchType,'Armijo')
    [step,xnew,fnew,gnew,jf,iExit] ...
        = armijoLS(prob,step,stepMax,f,g,p,x);
  end

  nf    = nf + jf;

  if  iExit > 3
    BadSearch = true; % Line search failure
  else

    % BFGS update.

    d      = step* p;
    y      = gnew - g;
    ytd    = y'*d;
    My     = M*y;
    ytMy   = y'*My;
    sigma  = ytd/(ytMy);

    normg  = norm(gnew);  normgMax = max( normgMax, normg );

    ResetH = itn > 0  &&  mod( itn, ResetFreq ) == 0;

    if ResetH
      Mjj  = sigma/max(1,normg);      Mjj =  min(max(Mjj,MjjMin),MjjMax);
      M    = diag(Mjj*ones(1,n));
      Info(4:4) = 'R';

    elseif ytd > - updateTol*(1 - wolfeTols(1))*step*gtp,
      %--------------------------------------------------
      % Update the inverse Hessian approximation.
      %--------------------------------------------------
      sd   = (1/ytd)*d;        sMy = (1/ytMy)*My;
      w    = sd - sMy;
      M    = M - sMy*My' + sd*d' + (ytMy*w)*w';
    else
      Info(3:3) = 'n';
      skipped   = skipped + 1;
    end
    x    = xnew;  f = fnew;  g = gnew;
  end
  itn    = itn + 1;
  if itn >= debugItn,
    jj   = 1; % Set break here
  end
end

if      Optimal
  status = 0;
elseif  Unbounded
  status = 1;
elseif  NearOptimal
  status = 2;
elseif  ItnLimit
  status = 3;
elseif  BadSearch
  status = 4;
end
