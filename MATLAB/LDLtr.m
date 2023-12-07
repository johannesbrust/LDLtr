function [x,itn,nf,skipped,status] = LDLtr(prob,outfile)
%LDLtr [x,itn,nf,skipped,status] = LDLtr(prob,outfile) finds
%        a local minimizer of a scalar-valued multivariate function
%        using a variant of the BFGS quasi-Newton method. 
%        The algorithm uses two phases. The first phase computes an
%        approximate optimal shift, while the second phase applies
%        a specific conjugate gradient iteration.
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
% 01/26/23, J.B., Initial version of WTRL2
% 02/10/23, J.B., Scaled trust-region subproblem
% 02/14/23, J.B., Updated scaled problem
% 02/16/23, J.B., Remove normalization
% 02/23/23, J.B., Add diagonal shift in subproblem
% 02/24/23, J.B., Low rank approximation in subproblem
% 03/01/23, J.B., Removing function evaluations
% 03/02/23, J.B., Use of inverse update, store and update only inverse
% 05/09/23, J.B., Function implements a CG solver in the subproblem
% 05/16/23, J.B., Further correction of CG in Newton solver
% 05/23/23, J.B., Use of a banded approximation for the TR subproblem
% 11/22/23, J.B., Preparation for release
% 11/27/23, J.B., Further updates for release
% 11/29/23, J.B., Updating for release
% 12/07/23, J.B., Preparation for release

format compact

% -------------------------------
% Assign local control parameters
% -------------------------------
parms         = bfgs_parms();

solver        = parms.solver       ;  % Solver
maxIterations = parms.maxIterations;  % maximum iterates allowed
tolStny       = parms.tolStat      ;  % stationarity       tolerance
lines         = parms.lines        ;  % lines between printing header
unBoundedf    = parms.unBoundedf   ;  % unbounded objective value
count_near    = parms.count_near   ;  % Option to count near-optimal results
print_run_det = parms.print_run_det;  % Option to print iteration details

%
% Initial line-search parameters
%
wolfeTols     = parms.wolfeTols    ;  % Wolfe function tolerance
dxMax         = parms.dxMax        ;  % maximum  change in x
dxInf         = parms.dxInf        ;  % infinite change in x
searchType    = parms.searchType   ;  % line search type

%
% Trust-region parameters
%
c1            = parms.c1           ; % step acceptance
c2            = parms.c2           ; 
c3            = parms.c3           ; 
c4            = parms.c4           ; % radius increase
c5            = parms.c5           ; 
c6            = parms.c6           ; 
c7            = parms.c7           ; % radius decrease
trTol         = parms.trTol        ; % tolerance for the TR subproblem
trMax         = parms.trMax        ; % maximum iterations for the TR subproblem

%
% Parameters for solving with triangular systems
%
optsUT.UT      = true              ; % Upper triangular system
optsUTT.UT     = true              ; % Transposed upper triangular
optsUTT.TRANSA = true              ;

% Parameter to branch to the MS subproblem solver
nmax            = parms.nmax       ; 

%
% CG subproblem tolerances (related to Alg. 3)
%
tolICG         = parms.tolICG      ; %1e-10; 1e-7
maxitICG       = parms.maxitICG    ; %40;
maxitBack      = parms.maxitBack   ; 
shrink         = parms.shrink      ;

% -----------------------
% End: control parameters
% -----------------------
if print_run_det
   str    = ' Solver' ;
   fprintf(        '%s%27s%15s\n', str, ':', solver);
   fprintf(outfile,'%s%27s%15s\n', str, ':', solver);

   str    = ' Line search' ;
   fprintf(        '%s%22s%15s\n', str, ':', searchType);
   fprintf(outfile,'%s%22s%15s\n', str, ':', searchType);
end
%
% Initialize the counters
%
itn      = 0;
nf       = 1;                 % number of function evaluations
skipped  = 0;
trIt     = 0;                 % Trust-region iterations
rho      = 0;                 % Sufficient decrease condition
phi      = 0;                 % Trust-region subproblem

numAccept = 0;
numTRInc  = 0;
numTRDec  = 0;

%
% Initialize the QN matrix
%

stepMax  = 10^5; %#ok<NASGU>
MjjMax   = 1.0d+4;            MjjMin = 1.0d-2;

x        = prob.x;
n        = prob.n;
[f,g]    = prob.obj(x);

normg    = norm(g);
normgMax = 1 + normg;
fMax     = 1 + abs(f);

Mjj      = 1/max(1,normg);      Mjj   =  min(max(Mjj,MjjMin),MjjMax);

% Representation of the LDLT factorization
ons     = ones(n,1);
In      = eye(n);
D       = (1/Mjj)*ons; %#ok<NASGU>

% Inverse quantities
R       = In;
DI      = Mjj*ons;
rr      = 1;
GI      = zeros(n,rr);
dn      = max(DI);

Info   = '     ';
header1 = '  Itn    Nf    Delta     jf    Objective      Norm g';
header = '  Itn    Objective    Norm g     Delta    Norm s    Rho      TRit   Phi     min(D)';


step   =  0;                  iExit = 1; %#ok<NASGU> % First line of output

Optimal  = false;   NearOptimal = false;   Unbounded = false;
ItnLimit = false;   BadSearch   = false;

% Start the quasi-Newton method based on an initial line-search
p         = -Mjj*g;
normp     = norm(p);

stepMax   = dxInf/(1+normp);
stepLimit = dxMax/(1+normp);
step      = min([1 stepLimit stepMax]);
if step  == stepLimit, Info(5:5) = 'l'; end

if     strcmp(searchType,'Wolfe' )
   [step,xnew,fnew,gnew,jf,iExit] ...
        = wolfeLSfg(prob,wolfeTols,step,stepMax,f,g,p,x); %#ok<*ASGLU>
elseif strcmp(searchType,'Armijo')
   [step,xnew,fnew,gnew,jf,iExit] ...
        = armijoLS(prob,step,stepMax,f,g,p,x);
else
    %
    % For testing the More-Thuente line-search
    % (Using the parameters from the Fortran implementation)
    %
    gtol             = 0.9;
    ftol             = 1.0e-4;
    xtol             = 1e-12;
    stpmin           = 1e-14;
    stpmax           = 20;
    maxfev           = 200; 
    [xnew,fnew,gnew,step,iExit,jf] = cvsrch(prob.obj,n,x,f,g,p,step,ftol,gtol,xtol, ...
                     stpmin,stpmax,maxfev);
end

nf    = nf + jf;

% quasi-Newton update
sk    = xnew - x;
yk    = gnew - g; 
nsk   = norm(sk);
nyk   = norm(yk);
sy    = (sk/nsk)'*(yk/nyk)*nsk*nyk;

% For testing
% x       = xnew;
% g       = gnew;
% f       = fnew;
% normg   = norm(g);

if sy > 0
    
    Hy = R*(DI.*(R'*yk));
    
    [ R, DI, GI ] = bfgsRDRUpdate_E1( R, DI, GI, sk, yk, Hy, nsk, nyk, sy );
    
end

% Radius initialization
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
    
    %----------------------------------------------------------------------
    % Print iteration details.
    %----------------------------------------------------------------------
    if print_run_det
        nLine    = rem( itn, lines );
        if nLine == 0
            fprintf('\n%s\n', header);  fprintf(outfile, '\n%s\n', header1);
        end
        
        %header = '  Itn    Objective    Norm g     Delta    Rho      TRit   Phi';
        
        if itn == 0
            cRho = '    --    ';
            cPhi = '    --    ';
            cDn  = '    --    ';
            cjf    =    '-- ';
            ctrIt  =    '-- ';
        else
            cRho  = sprintf ( ' %9.2e', rho );
            cPhi  = sprintf ( ' %9.2e', abs(phi) );
            cDn   = sprintf ( ' %2.2e', dn );
            cjf   = sprintf ( '%2g '  ,   jf );
            ctrIt = sprintf ( '%2g '  ,   trIt );
        end
        if iExit ~= 1
            cExit = sprintf ( '(%1g)', iExit );
        else
            cExit = '   ';
        end
        
        cstep = sprintf ( ' %9.2e %9.2e', Delta, nsk );
        
        str1  = sprintf ( '%5g %5g',        itn,      nf );
        str2  = sprintf ( '%s%s%s',               cstep,  cExit, cjf );
        str3  = sprintf ( '%14.7e %9.2e', f, norm(g) );
        % str3  = sprintf ( '%14.7e %9.2e %8.1e', f, norm(g), cond(M) );
        str4  = Info;
        
        % Trust-region infos
        citn  = sprintf ( '%5g', itn );
        str5  = [ citn, str3, cstep, cRho, ctrIt, cPhi, cDn ];
        
        str   = [ str1 str2 str3 str4 ];
        fprintf('%s\n', str5);    fprintf(outfile, '%s\n', str);
        Info  = '     ';
    end
    
    if  Optimal || Unbounded || NearOptimal || ItnLimit || BadSearch
        break;
    end
    

    %
    % Computation of a full quasi-Newton step sk, and check
    % whether a trust-region subproblem needs to be solved.
    % Trial function values for later comparisons are also computed,
    % as is the smallest value in the diagonal d stored.
    %    
    D       = 1./DI;    
    hk      = R'*g;     
    sig     = 0;
    vk      = -(hk ./ (D+sig));    
    sk      = R*(DI.*(-hk)); %linsolve(L,-hk./D,optsLTT);    
    nsk     = norm(sk);    
    dn      = min(D);
   
    % Check if the subproblem needs to be solved
    if nsk - Delta > 0; solveTRSUB = 1; else; solveTRSUB = 0; end
    
    % Store a copy of the quasi-Newton step for later comparisons
    sk1     = sk;
    sk1_    = sk;
    
    [fk1t1,gk1t1] = prob.obj(x+sk1);
    ared1   = f - fk1t1;            
    ngk1t1  = norm(gk1t1);
    
    if ngk1t1 <= tolStny && Delta < tolStny
        solveTRSUB = 0;
    end
    
    %
    % Solve the trust-region subproblem using Alg. 2
    % This is Newton's method on the 1 dimensional secular equation
    %
    %
    trIt = 0;
    if solveTRSUB == 1 && nmax < n
        
        phi      = 1/nsk - 1/Delta;
        relphi   = (Delta-nsk)/Delta;
        
        %
        % Sparse matrix, and solve (an alternative way to compute the step)
        %
        % SD      = spdiags(D,0,n,n);        
        % AA1     = spdiags(GI(:,1),0,n,n);        
        % vkp     = (SD + sig.*AA1)\(-(AA1*vk));
        
        % Vector operation
        vkp     = (-GI(:,1).*vk) ./ (D + sig.* GI(:,1));

        skp     = R*(vkp);
        
        while (abs(relphi) > trTol) && (trIt < trMax)
        
            skskp      = (sk'*skp);
            sig        = sig + ((nsk*nsk)/skskp)*((Delta-nsk)/Delta);
            
            % Sparse matrix computation
            % vk       = (SD + sig.*AA1)\(-hk);
        
            vk         = -hk ./ (D + sig.*GI(:,1));            
            sk         = R*(vk);%linsolve(L,(vk),optsLTT);

            % Evaluate the function to check progress
            if normg < 1e0
                [fk1t2,gk1t2] = prob.obj(x+sk);
                ared1_ = f - fk1t2;

                if ared1_ > ared1 && sig > 0

                    ared1 = ared1_;
                    sk1 = sk;
                    fk1t1 = fk1t2;
                    gk1t1 = gk1t2;

                end

                ngk1t1 = norm(gk1t1);

                if ngk1t1 <= tolStny && Delta < tolStny

                    ared1 = ared1_;
                    sk1 = sk;
                    fk1t1 = fk1t2;
                    gk1t1 = gk1t2;

                    break;

                end
            end
            
            nsk        = norm(sk);
                  
            %vkp        = (SD + sig.*AA1)\(-(AA1*vk));
            vkp        = (-GI(:,1).*vk) ./ (D + sig.* GI(:,1));
                        
            skp        = R*(vkp);%linsolve(L,(vkp),optsLTT);                        
            phi        = 1/nsk - 1/Delta;            
            relphi     = (Delta-nsk)/Delta;            
            trIt       = trIt + 1;
            
        end
        
        % (For safeguarding)
        % Bi-section method if Newton iteration 
        % does not find a minimizer
        skN = sk;
        relphiN = relphi;
        sigN = sig;
        
        if sig < 0 || trIt == trMax
            sk   = sk1_;
            nsk  = norm(sk); 

            sigb    = 10*max(D);
            
            %vkb          = (SD + sigb.*AA1)\(-hk);
            vkb     = -hk ./ (D + sig.*GI(:,1));
            
            skb     = R*(vkb);%linsolve(L,(vkb),optsLTT);            
            phib    = 1/norm(skb) - 1/Delta;

            itPhib = 2;
            while (phib < 0) && itPhib < 10

                sigb        = 10^(itPhib)*max(D);
                
                %vkb          = (SD + sigb.*AA1)\(-hk);
                vkb          = -hk ./ (D + sig.*GI(:,1));
                
                skb         = R*(vkb); %linsolve(L,(vkb),optsLTT);        
                phib        = 1/norm(skb) - 1/Delta;

                itPhib      = itPhib + 1; 

            end

            siga        = 0;
            phia        = 1/nsk - 1/Delta;                        
            phi         = phia;
            
            relphi      = (Delta-nsk)/Delta;

            while (abs(relphi) > trTol) && (trIt < 2*trMax)

                sig     = (siga+sigb)/2;

                %vk      = (SD + sig.*AA1)\(-hk);
                vk       = -hk ./ (D + sig.*GI(:,1));
                
                sk      = R*(vk); %linsolve(L,(vk),optsLTT);
                nsk     = norm(sk);

                phi     = 1/nsk - 1/Delta;                
                relphi  = (Delta-nsk)/Delta;
                
                % Check intermediate improvements
                % Evaluate the function to check progress
                if normg < 1e0
                    [fk1t2,gk1t2] = prob.obj(x+sk);
                    ared1_ = f - fk1t2;
                    %Ared2 = ared1_;

                    if ared1_ > ared1 || ...
                    abs(ared1_)/abs(f) < eps && normg > norm(gk1t2) && Delta < tolStny%&& ared1_ > 0 

                        ared1 = ared1_;
                        sk1 = sk;
                        fk1t1 = fk1t2;
                        gk1t1 = gk1t2;

                    end

                    ngk1t1 = norm(gk1t2);

                    if ngk1t1 <= tolStny && Delta < tolStny

                        ared1 = ared1_;
                        sk1 = sk;
                        fk1t1 = fk1t2;
                        gk1t1 = gk1t2;

                        break;

                    end
                end

                if phi < 0
                    siga = sig;
                else
                    sigb = sig;
                end

                trIt = trIt + 1;

            end
        end
        
        % Pick Newton step if Bisection does not improve
        if sigN >= 0 && abs(relphiN) < abs(relphi)

            sk = skN;
            sig = sigN;

        end
    
    elseif solveTRSUB == 1 && n <= nmax
        
        %
        % More-Sorensen trust-region algorithm
        %

        B          = linsolve(R,In,optsUT);
        B          = linsolve(R,diag(D)*B,optsUTT);
        
        [sk,sig,~,~,trIt,~]= solveTrMedium(B,g,Delta,sig,parms);
            
    end
    
    %    
    % CG iteration with backtracking strategy
    %
    sigCG           = sig;
    [fktCG,gktCG]   = prob.obj(x+sk);  nf = nf + 1;  
    ngktCG          = norm(gktCG);
    itMin           = 0;
    hasGRADIM       = 0;

    if solveTRSUB == 1
        for i = 1:maxitBack

            %
            % CG iteration for (sig * R'R + D) vk = -hk
            %
            vkCG           = inv_ldl_CG1( sigCG, -hk, R, D, tolICG, maxitICG );                        
            skCG           = R*(vkCG); %linsolve(L,(vk),optsLTT);

            [fk1t2,gk1t2]  = prob.obj(x+skCG); nf = nf + 1;            
            ngk1t2         = norm(gk1t2);

            % If the function decreases or roundoff errors are present
            % and gradient norm decreases the step is saved
            if fk1t2 < fktCG && abs(fk1t2-f)/abs(f) > 1e-12
                sk = skCG;
                fktCG = fk1t2;
                itMin= i;
            elseif abs(fk1t2-f)/abs(f) < 1e-12 && ngk1t2 < ngktCG && ...
                Delta < 1e-7
                sk = skCG;
                fktCG = fk1t2;
                ngktCG = ngk1t2;
                itMin= i;
                hasGRADIM = 1;
            end

            sigCG = shrink*sigCG;

        end
    end

    %
    % Evaluate the quadratic model and find rho for 
    % the step from computing sigma (ie., Alg. 1)
    %
    Bs1          = linsolve(R,sk1,optsUT);
    Bs1          = linsolve(R,D.*Bs1,optsUTT);
    sBs1         = sk1'*Bs1;
    gs1          = g'*sk1;
    
    pred1        = -(gs1 + 0.5*sBs1);    
    rho1         = ared1/pred1;
        
    if ngk1t1 <= tolStny && Delta < tolStny
        rho1 = 1;
    end
        
    % Safeguard in case the model has a negative optimal value
    if pred1 < 0        
        sk1 = zeros(n,1);
        rho1 = -1;
    end
    
    %
    % Compute a trial step
    %
    xk1t        = x+sk;
    [fk1t,gk1t] = prob.obj(xk1t);        
    nf          = nf + 1;    
    ared        = f - fk1t;
    
    %
    % Evaluate the quadratic model and find rho
    %
    Bs          = linsolve(R,sk,optsUT);
    Bs          = linsolve(R,D.*Bs,optsUTT);
    sBs         = sk'*Bs;
    gs          = g'*sk;    
    pred        = -(gs + 0.5*sBs);
    
    rho         = ared/pred;

    %
    % Additional acceptance for small changes
    % This can accept steps that reduce the gradient norm
    % when the computed values are close to round-off errors
    %
    if abs(ared)/abs(f) < 1e-12 && Delta < 1e-7 && ...
        rho1 ~= 1 && hasGRADIM == 1 % && norm(g) < 1e0

        if norm(gk1t) < norm(g)
            rho = c5;
        end        
    end
    
    if pred < 0 %&& ared < 0
        sk = zeros(n,1);        
        rho = -1;
    end
        
    %
    % Use the best of availabe steps in terms of rho
    %
    if rho1 > rho        
        rho = rho1;
        gk1t = gk1t1;
        fk1t = fk1t1;
        sk   = sk1;
        xk1t = x+sk;                
    end
    
    yk  = gk1t - g;
    
    nsk = norm(sk);
    nyk = norm(yk);
    
    %
    % Trust-region parameter updates
    % Parameters c1, c2, c3, c4, c5, c6 and c7 are all 
    % related to the trust-region step
    % The updates to the shrink parameter are adjustments to 
    % the value gamma in Alg. 3
    %

    if c1 < rho

        % The trial step is accepted
        xnew = xk1t;
        gnew = gk1t;
        fnew = fk1t;
        
        numAccept = numAccept + 1;

        %
        % For acceptable trial steps, the shrink (gamma) scaling
        % for the CG iteration can be updated. This
        % quantity is doubled or halved depending on the outcome
        %
        if itMin == maxitBack
            shrink = 0.5*shrink;        %0.25*shrink;
        elseif itMin <= 1
            shrink = (1/0.5)*shrink;    %(1/0.25)*shrink;
        end
        shrink = min(shrink,0.25);
        shrink = max(shrink,0.25^10);
        
    else
        xnew = x;
        gnew = g;
        fnew = f;
    end
    
    %
    % Depending on the improvement through the trial step
    % the radius Delta is updated
    %
    if c2 < rho
        
        if nsk <= c3*Delta %eps %c3*Delta %nsk <= c3*Deltak
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
    
    %
    % BFGS update.
    %
    sy = sk'*yk;
    
    normg  = norm(gnew);  normgMax = max( normgMax, normg );
        
    if sy > 0 && rho >= 0 %> c1
        
        %
        % Update the Hessian approximation.
        %        
        Hy = R*(DI.*(R'*yk));
        
        [ R, DI, GI ] = bfgsRDRUpdate_E1( R, DI, GI, sk, yk, Hy, nsk, nyk, sy );
        
    else
        Info(3:3) = 'n';
        skipped   = skipped + 1;
    end

    %
    % Update the iteration
    %
    x       = xnew;  f = fnew;  g = gnew;
    itn     = itn + 1;
    
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

end
