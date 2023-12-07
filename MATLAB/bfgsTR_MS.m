function [x,itn,nf,skipped,status] = bfgsTR_MS(prob,outfile)
%bfgsTR_MS  [x,itn,nf,skipped,status] = bfgsTR_MS(prob,outfile) finds
%        a local minimizer of a scalar-valued multivariate function
%        using a variant of the BFGS quasi-Newton method with a
%        trust-region subproblem.
%
%        prob is a structure that provides the details of the
%        CUTEst problem being solved.
%
%        On successful termination, x is the first point found at
%        which the norm of the gradient is less than a tolerance.
%
%        Input arguments are set in the file bfgs_parms.m.
%
%        The meaning of each line of output is as follows
%
%        Itn          The current iteration number.
%
%        f(x)         The current value of the objective function f(x).
%
%        Norm g       The current value of the gradient of f(x). In
%                     the neighborhood of a local minimizer, Norm g
%                     will be small.
%
%        Delta        The radius of the trust-region.
%
%        Rho          Sufficient decrease metric: (fk1-fk)/Q(sk)
%
%        TRit         Trust-region subproblem iterations
%
%        Phi          Subproblem constraint: |1/Delta - 1/norm(sk)|
%
%--------------------------------------------------------------------------
% This version extends an implementation from 29 Feb 2020 by Philip E. Gill
%
% The new version is first implemented 24 Jan 2023 by Johannes J. Brust
% University of California, San Diego.
%
% Updates: 

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

c1            = 1e-5               ; % step acceptance
c2            = 0.75               ; 
c3            = 0.8                ; 
c4            = 2                  ; % radius increase
c5            = 0.1                ; 
c6            = 0.75               ; 
c7            = 0.5                ; % radius decrease
  
trTol         = 1e-3               ; % tolerance for the TR subproblem
trMax         = 20                 ; % maximum iterations for the TR subproblem

optsSPD.SYM   = true               ; % Solving with a spd system
optsSPD.POSDEF= true               ;

optsUT.UT      = true              ; % Upper triangular system
optsUTT.UT     = true              ; % Transposed upper triangular
optsUTT.TRANSA = true              ;

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
trIt     = 0;                 % Trust-region iterations
rho      = 0;                 % Sufficient decrease condition
phi      = 0;                 % Trust-region subproblem

numAccept = 0;
numTRInc  = 0;
numTRDec  = 0;

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
B      = diag((1/Mjj)*ones(1,n));
In     = eye(n);

Info   = '     ';
% header = '  Itn    Nf    Step     jf    Objective      Norm g   cond H';
header1 = '  Itn    Nf    Delta     jf    Objective      Norm g';
header = '  Itn    Objective    Norm g     Delta    Rho      TRit   Phi';


step   =  0;                  iExit = 1; % First line of output
dType  = SD;

Optimal  = false;   NearOptimal = false;   Unbounded = false;
ItnLimit = false;   BadSearch   = false;

% Start the quasi-Newton method based on an initial line-search
p         = -Mjj*g;
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

% quasi-Newton update
sk    = xnew - x;
yk    = gnew - g; 
sy    = sk'*yk;

if sy > 0
    Bs = B*sk;
    B  = B - Bs*(Bs'/(sk'*Bs)) + yk*(yk'/sy);
end

nsk     = norm(sk);
Delta   = 2*nsk;

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
            fprintf('\n%s\n', header);  fprintf(outfile, '\n%s\n', header1);
        end
        
        %header = '  Itn    Objective    Norm g     Delta    Rho      TRit   Phi';
        
        if itn == 0
            cRho = '    --    ';
            cPhi = '    --    ';
            cjf    =    '-- ';
            ctrIt  =    '-- ';
        else
            cRho  = sprintf ( ' %9.2e', rho );
            cPhi  = sprintf ( ' %9.2e', abs(phi) );
            cjf   = sprintf ( '%2g '  ,   jf );
            ctrIt = sprintf ( '%2g '  ,   trIt );
        end
        if iExit ~= 1
            cExit = sprintf ( '(%1g)', iExit );
        else
            cExit = '   ';
        end
        
        cstep = sprintf ( ' %9.2e', Delta );
        
        str1  = sprintf ( '%5g %5g',        itn,      nf );
        str2  = sprintf ( '%s%s%s',               cstep,  cExit, cjf );
        str3  = sprintf ( '%14.7e %9.2e', f, norm(g) );
        % str3  = sprintf ( '%14.7e %9.2e %8.1e', f, norm(g), cond(M) );
        str4  = Info;
        
        % Trust-region infos
        citn  = sprintf ( '%5g', itn );
        str5  = [ citn, str3, cstep, cRho, ctrIt, cPhi ];
        
        str   = [ str1 str2 str3 str4 ];
        fprintf('%s\n', str5);    fprintf(outfile, '%s\n', str);
        Info  = '     ';
    end
    
    if  Optimal || Unbounded || NearOptimal || ItnLimit || BadSearch
        break;
    end
    
    % Trust-region subproblem solution
    % Try unconstrained minimizer first
    
    sk = linsolve(B,-g,optsSPD);
    nsk = norm(sk);
    
    if nsk - Delta > 0; solveTRSUB = 1; else solveTRSUB = 0; end;
    
    sig = 0;
    
    trIt  = 0;
    if solveTRSUB == 1
        
        [sk,~,~,~,trIt,~]= solveTrMedium(B,g,Delta,sig,parms);
        
        nsk = norm(sk);
        
        phi = 1/nsk - 1/Delta;
        
%         [R,p1]  = chol(B);
%         
%         % Ensure the matrix is positive definite
%         %       ii = 0; fac = 1e-3;
%         %       while p > 0
%         %           [R,p]  = chol(B+fac*In);
%         %           fac    = fac * 10;
%         %           ii     = ii + 1;
%         %       end
%         %
%         phi       = 1/nsk - 1/Delta;
%         q         = linsolve(R,sk,optsUTT);
%         nq        = norm(q);
%         
%         while (abs(phi) > trTol) && (trIt < trMax)
%             
%             sig        = sig + (nsk/nq)^2*((nsk-Delta)/Delta);
%             
%             [R,p]      = chol(B+sig.*In);
%             
%             sk         = linsolve(R,-g,optsUTT);
%             sk         = linsolve(R,sk,optsUT);
%             q          = linsolve(R,sk,optsUTT);
%             
%             nsk        = norm(sk);
%             nq         = norm(q);
%             
%             phi        = 1/nsk - 1/Delta;
%             trIt       = trIt + 1;
%             
%         end
        
    end
    
    % Compute a trial step
    xk1t        = x+sk;
    [fk1t,gk1t] = prob.obj(xk1t);    
    yk          = gk1t - g;
    
    jf          = jf + 1;
    
    ared        = f - fk1t;
    Bs          = B*sk;
    sBs         = sk'*Bs;
    gs          = g'*sk;
    
    pred        = -(gs + 0.5*sBs);
    
    rho         = ared/pred;
    
    if pred < 0
        sk = 1e-12.*ones(n,1);
        rho = -1;
    end
    
    nsk = norm(sk);
    
    if c1 < rho
        xnew = xk1t;
        gnew = gk1t;
        fnew = fk1t;
        
        numAccept = numAccept + 1;
        
        %isAccept = 1;
        
    else
        xnew = x;
        gnew = g;
        fnew = f;
        
        %isAccept = 0;
    end
    
    if c2 < rho
        
        if nsk <= c3*Delta %nsk <= c3*Deltak
            Delta = 1*Delta;
        else
            Delta = c4*Delta;
            numTRInc = numTRInc + 1;
        end
        
    elseif (c5 <= rho) && (rho <= c6)
        Delta = 1*Delta;
    else
        Delta = c7*Delta;
        numTRDec = numTRDec + 1;
    end
    
    % BFGS update.
    
    sy = sk'*yk;
    
    normg  = norm(gnew);  normgMax = max( normgMax, normg );
    
    ResetH = itn > 0  &&  mod( itn, ResetFreq ) == 0;
    
    if ResetH
        Mjj  = sigma/max(1,normg);      Mjj =  min(max(Mjj,MjjMin),MjjMax);
        B    = diag((1/Mjj)*ones(1,n));
        Info(4:4) = 'R';
        
    elseif sy > 0 && rho > c1
        %--------------------------------------------------
        % Update the Hessian approximation.
        %--------------------------------------------------
        
        B = B - Bs*(Bs'/sBs) + yk*(yk'/sy);
        
    else
        Info(3:3) = 'n';
        skipped   = skipped + 1;
    end
    x    = xnew;  f = fnew;  g = gnew;
    itn    = itn + 1;
    if itn >= debugItn,
        jj   = 1; % Set break here
    end
    if Delta < 1e-22
        BadSearch = true;
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
