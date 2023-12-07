function [x,itn,nf,skipped,status] = bfgsR(prob,outfile)
%bfgsR   [x,itn,nf,skipped,status] = bfgsR(prob,outfile) finds
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
%        This version dated Jan 22 2020.
%        Philip E. Gill, University of California, San Diego.

  format compact

% -------------------------------
% Assign local control parameters
% -------------------------------
  parms         = bfgs_parms();

  solver        = parms.solver       ;  % Solver
  maxIterations = parms.maxIterations;  % maximum iterates allowed
  tolStny       = parms.tolStat      ;  % stationarity tolerance
  wolfeTols     = parms.wolfeTols    ;  % Wolfe function tolerance
  condMax       = parms.condMax      ;  % maximum condiion estimate
  dxMax         = parms.dxMax        ;  % maximum change in x
  dxInf         = parms.dxInf        ;  % infinite change in x
  updateTol     = parms.updateTol    ;  % Tolerance for rejecting BFGS update
  lines         = parms.lines        ;  % lines between printing header
  ResetFreq     = parms.ResetFreq    ;  % Reset frequency
  infinity      = parms.infinity     ;  % user-defined + infinity
  Assert        = parms.Assert       ;  % (0/1) check consistency
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
  RjjMax   = 1.0d+1;            RjjMin = 1.0d-2;

  x        = prob.x;
  n        = prob.n;
  [f,g]    = prob.obj(x);

  normg    = norm(g);
  normgMax = 1 + normg;
  fMax     = 1 + abs(f);

  Rjj    = sqrt(max(1,normg));  Rjj   =  min(max(Rjj,RjjMin),RjjMax);
  R      = diag(Rjj*ones(1,n)); condR = 1;

  Info   = '     ';
  header = '  Itn    Nf    Step     jf    Objective      Norm g';
  step   =  0;                   iExit = 1; % First line of output
  dType  = SD;

  Optimal  = false;   NearOptimal = false;   Unbounded = false;
  ItnLimit = false;   BadSearch   = false;
  optsU.UT = true;    optsL.LT    = true;    % For linsolve

% --------------------------------------------------------------------------
% Main loop.
% --------------------------------------------------------------------------
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
      Rjj     = sqrt(max(1,normg));  Rjj   =  min(max(Rjj,RjjMin),RjjMax);
      R       = diag(Rjj*ones(1,n)); condR = 1;
      dType   = SD;              Info(5:5) = 'R';
      BadSearch = false;
    end

  % ------------------------------------------------------------------------
  % Print iteration details.
  % ------------------------------------------------------------------------
    if print_run_det
      nLine    = rem( itn, lines );
      if nLine == 0
        fprintf('\n%s\n', header);  fprintf(outfile, '\n%s\n', header);
      end

      if itn == 0
        cstep = '    --    ';
        cjf   =    '-- ';
      else
        cstep = sprintf ( ' %9.2e', step );
        cjf   = sprintf ( '%2g '  ,   jf );
      end
      if iExit ~= 1
        cExit = sprintf ( '(%1g)', iExit );
      else
        cExit = '   ';
      end
      str1  = sprintf ( '%5g %5g',            itn,      nf );
      str2  = sprintf ( '%s%s%s',             cstep,  cExit, cjf );
      str3  = sprintf ( '%14.7e %9.2e', f, norm(g));
      str4  = Info;
      str   = [ str1 str2 str3 str4 ];
      fprintf('%s\n', str);    fprintf(outfile, '%s\n', str);
      Info  = '     ';
    end

    if  Optimal || Unbounded || NearOptimal || ItnLimit || BadSearch
      break;
    end

    q   = -linsolve(R', g, optsL);  p     =  linsolve(R, q, optsU);
    gtp =  g'*p;                    normp = norm(p);

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

      d      = step*p;
      y      = gnew - g;
      ytd    = y'*d;

      normg  = norm(gnew);  normgMax = max( normgMax, normg );

%     if ytd > - updateTol*(1 - wolfeTols(1))*step*gtp,
      if ytd > 0 %- (1 - updateTol)*step^2*gtp,
    %   --------------------------------------------------
    %   Update the Cholesky factor R such that H = R'*R;
    %   --------------------------------------------------
        normq = norm(q);
        u     = q/normq;
        if  y'*g > 0,
          v   = -y/sqrt(ytd) + g/normq;
        else
          v   =  y/sqrt(ytd) + g/normq;
        end

        k  = find(u, 1, 'last');
        v_ = linsolve(R', v, optsL);
        R  = modR(R, u, v_, k);
        Rjj   = abs(diag(R));             condR = max(Rjj)/min(Rjj);
        if condR*condR > condMax
          gamma  = sqrt(ytd)/(step*normq);
          Rjj    = gamma*sqrt(max(1,normg)); Rjj   =  min(max(Rjj,RjjMin),RjjMax);
          R      = diag(Rjj*ones(1,n));      condR = 1;
          dType  = SD;                   Info(4:4) = 'M';
        else
          dType = QN;
        end
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
end
%--------------------------------------------------------------------------
function   [R] = modR(R, u, v_, k)
  [bet_, s]       = rec1(u, k);
  [lam, bet, gam] = rec2(bet_, s, u, v_, k);
  R               = rec3(R, bet, u, gam, v_, lam);
end
%--------------------------------------------------------------------------
function [bet_, s] = rec1(u, k)
  n    = length(u);
  bet_ = zeros(n, 1);
  s    = zeros(n, 1);
  bet_(k) = 1/u(k);
  for j = k-1:-1:1
    rj        = bet_(j + 1)*u(j);
    s(j)      = 1/sqrt(rj^2 + 1);
    bet_(j)   =    s(j)*bet_(j + 1);
    bet_(j+1) = rj*s(j)*bet_(j + 1);
  end
end
%--------------------------------------------------------------------------
function [lam, bet, gam] = rec2(bet_, s, u, v, k)
  n    = length(u);

  lam  =  ones(n, 1);
  bet  = zeros(n, 1);
  gam  = zeros(n, 1);

  gam_ = 1/(bet_(1));
  lam_ = bet_(1)*u(1) + gam_*v(1);

  for j = 1:n-1
    lam(j)      = sqrt(lam_^2 + s(j)^2);
    c_          = lam_/lam(j);             s_ = s(j)/lam(j);
    bet(j)      = c_*bet_(j) - s_*bet_(j+1);
    gam(j)      = c_*gam_;
    gam_        = s_*gam_;
    bet_(j + 1) = s_*bet_(j) + c_*bet_(j + 1);
    lam_        = bet_(j + 1)*u(j + 1) + gam_*v(j + 1);
    if j == k
      return;
    end
  end
  lam(n) = lam_;
end
%--------------------------------------------------------------------------
function [R] = rec3(R, alph, q, bet, p, lam)
  n = length(bet);
  w = zeros(n, 1);
  z = zeros(n, 1);

  for i = n:-1:1
    w(i)    =   p(i)*R(i, i);
    z(i)    =   q(i)*R(i, i);
    R(i, i) = lam(i)*R(i, i);
    for j = i + 1:n
      r        = R(i, j);
      R(i , j) = lam(i)*r + alph(i)*z(j) + bet(i)*w(j);
      w(j)     = w(j) + p(i)*r;
      z(j)     = z(j) + q(i)*r;
    end
  end
end
