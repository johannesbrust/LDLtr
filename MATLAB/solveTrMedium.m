function [s,sigma,Qs,Solved,itns,factored]=solveTrMedium(H,g,delta,sigma,parms)
%        [s,sigma,Qs,Solved,itns,factored]=solveTrMedium(H,g,delta,sigma,parms)
% =========================================================================
% Uses the More'-Sorenson algorithm to find an approximate solution of the
% quadratically constrained trust-region subproblem:
%
%    min Q(s) = g'*s + (1/2)s'*H*s  subject to twoNorm(s) <= delta.

% On input:
%    H              is the Hessian of the quadratic model
%    g              is the linear term of the quadratic model
%    delta          is the positive trust-region radius
%    sigma          is an initial guess for the optimal sigma
%    relTol, absTol are relative and absolute tolerances controlling the
%                   accuracy of the approximate solution.
%
% On output:
%    s              is an approximate minimizer of Q(s)
%    sigma          approximates the optimal sigma
%    Solved         = 3 if stopped when width is small (p+tau*z)
%                   = 2 if stopped when width is small (p)
%                   = 1 if an approximate solution has been found
%                   = 0 if no solution could be found (giving s=0)
%    Qs             is the final value of Q(s)
%    itns           is the number of iterations
%    factored       is the number of Cholesky factorizations needed

% The method finds an approximate zero of  psi(sigma), where
%              psi(sigma) = 1/delta - 1/twonorm(s(sigma)) with
%   (H + sigma I)s(sigma) = - g.

%  Local variables:
%  ---------------
%   sigmaL, sigmaU  define the current interval of uncertainty for sigma.
%   sigmaS          is a lower bound on -lambda_n, the least eigenvalue of H.
%   sigmaN          is the Newton estimate of sigma.

% The optimality conditions are:
%    sigma (twonorm(s) - delta) = 0, with sigma >=0,   twoNorm(s) <= delta
%   (H + sigma I)s =-g               with (H + sigma I) positive semidefinite

% First version dated 09-Jun-2008
% This  version dated 22-Oct-2022
%
% Philip Gill and Joey Reed
% University of California, San Diego.
% =========================================================================

  tic % Start the matlab timer

% -------------------------------
% Assign local control parameters
% -------------------------------
  relTol =  parms.TRrelTol  ;  % trust-region relative tolerance.
  absTol =  parms.TRabsTol  ;  % trust-region absolute tolerance.
  TRlog  =  parms.TRlog     ;  % trust-region log off/on
  itnMax =  parms.TRitnMax  ;  % trust-region iteration limit

  tolg       = sqrt(eps); % Input absolute termination tol for g
  if nargin > 4
    tolrel   = relTol;  tolabs = absTol;
  else
    tolrel   = 1e-6;    tolabs = 0;
  end
  nullVecTol = tolrel*(2 - tolrel);
  tolmin     = min(tolrel, 1e-9 );
  lowBias    = 0.999;

  factored   = 0;
  Solved     = 0;

  delsq      = delta^2;

  [n,m]      = size(H);
  e          = ones(n,1);

  normg      = norm(g);        gtg = normg^2;     normH = norm(H,1);

  eigMin     = min(diag(H));
  sigmaS     = max([ -eigMin, (normg/delta)-normH ]);

% Compute the initial value for sigmaL and sigmaU.
% sigmaL is set negative because any sigma <= sigmaL is rejected without
% factoring, but we want to factor if the first sigma is 0.
% If the inertia for a zero sigma is wrong, we set sigmaL = sigma = 0.
% If sigma is 0 on subsequent iterations, the matrix will not be refactored
% (sigma is always nonnegative).

  if sigmaS > 0
    sigmaL = sigmaS;
  else
    sigmaL = -1;
  end
  sigmaU   = (normg/delta) + normH;

  if sigma >= sigmaL  &&  sigma <= sigmaU
    sigmaN  = sigma;
  else
    sigmaN = min( sigmaU, max([ sigmaL, 0 ]));
  end

% sigmaN   = max([ sigmaL, sigma 0 ]);
% sigmaN   = min([ sigmaU, sigma   ]);

  relphi   = inf;
  s        = zeros(n,1);
  if TRlog
    header1 = sprintf ('    TRitn     sigma      sigmaN      sigmaS');
    header2 = sprintf ('      sigmaL      sigmaU   pos  TR resid\n');
    header  = [ header1 header2 ];
    fprintf('\n%s', header);
  end

  for  itns = 1:itnMax
  % --------------------------------
  % Safeguard sigmaN
  % --------------------------------
    if sigmaN <= sigmaL

    % We are outside the bounds of our safeguarding algorithm.
    % sigmaN = lowBias*sigmaL + (1 - lowBias)*sigmaU is another option

      sigmaN = max([(1 - lowBias)*sigmaU,sqrt(sigmaL*sigmaU)]);
    end

  % sigmaN cannot be too close to sigmaL;

    tolL = tolmin*(1 + abs(sigmaL));
    if sigmaN < sigmaL +  tolL
      sigmaN  = sigmaL +  tolL;
    end

  % --------------------------------
  % Factor H + sigma I
  % --------------------------------
    Hsigma    = H + diag(sigmaN*e);
    [R,Rtype] = chol(Hsigma);         factored = factored + 1;
    posDef    = Rtype == 0;

    if TRlog
      fprintf ( '      %3g %11.4e %11.4e %11.4e %11.4e %11.4e %3g', ...
                itns, sigma, sigmaN, sigmaS, sigmaL, sigmaU, posDef );
      if relphi ~= inf
        fprintf ( ' %9.2e\n', relphi );
      else
        fprintf ( '\n' );
      end
    end

    if ~posDef
      neg    =  n     - Rtype;         pos     =  Rtype - 1;

      r      =  R'\Hsigma(1:pos,pos+1);
      d      = -(Hsigma(pos+1,pos+1) - r'*r);
      u      = -R\r;

      sigmaS = max([sigmaS,sigmaN + d/(u'*u + 1)]);
      sigmaL = max([sigmaL,sigmaS]);
      sigmaN = max((1-lowBias)*sigmaU,sqrt(sigmaL*sigmaU));
    % sigmaN = (sigmaL + sigmaU)/2 is another option

    else
      sigma  = sigmaN;
      w      = R'\(-g);                Rp = w;         normRp = sqrt(Rp'*Rp);
      p      = R\w;                 normp = norm(p);      ptp = normp^2;
      relphi = (normp - delta)/delta;

      if  abs(relphi) <= tolrel
        s = p;                     Solved = 1;          break;
      else
        if  relphi < 0
          if  sigma == 0
            s = p;                 Solved = 1;          break;
          else
            [z]    = approxNull(R);    ptz = p'*z;     negphi = delsq - ptp;
            if  ptz >= 0
              tau  = negphi/(ptz + sqrt(ptz^2 + negphi));
            else
              tau  = negphi/(ptz - sqrt(ptz^2 + negphi));
            end
            Rz     = R*z;                              normRz = norm(Rz);
            sigmaS = max(sigmaS, sigma-normRz^2);
            if  (tau*normRz)^2 <= nullVecTol*max([tolabs,normRp^2+sigma*delsq])
              s = p + tau*z;        Solved = 1;         break;
            end

          % The width is too small

            if sigmaU - sigmaL <= 2*tolL
              p_z = p + tau*z;
              if p'*g + p'*H*p/2 > p_z'*g +  p_z'*H*p_z/2
                Solved = 3;
                s      = p_z;
              else
                Solved = 2;
                s      = p;
              end
              break;
            end
          end
          sigmaU = sigma;
        else
          sigmaL = sigma;
        end
      % ------------------------------------------
      % Compute a new sigma using Newton's method.
      % ------------------------------------------
        if  normg >= tolg
          q      = R'\p;   normq = norm(q);  qtq = normq^2;
          sigmaN = sigma + relphi*ptp/qtq;
        else
          sigmaN = max([ 0 sigmaS ]);
        end
        sigmaN   = max([ 0 sigmaN ]);
      end
    end
  end

  Qs     = s'*g + s'*H*s/2;
  TRtime = toc;
end
%--------------------------------------------------------------------------
function [z] = approxNull(R)
%        [z] = approxNull(R)
%               returns an approximate null vector of A based on the
%               LINPACK algorithm.  The steps of this algorithm are the
%               following:
%       1.  Solve R'w = e where e has components 1 or -1.  We
%           choose e to make w maximal the sense described below.
%       2.  Solve Rv=w
%       3.  z = v/||v||, and this is the approximate null vector.

  [n,m] = size(R);
  p     = zeros(n,1);   w = p;

  for k = 1:n,
    wPos = (1-p(k))/R(k,k);  wNeg = -(1+p(k))/R(k,k);

    sumPos = 0;  sumNeg = 0;

    for j=k+1:n
      sumPos = sumPos + abs(p(j)+R(k,j)*wPos);
      sumNeg = sumNeg + abs(p(j)+R(k,j)*wNeg);
    end

    sumPos = abs(wPos) + sumPos;
    sumNeg = abs(wNeg) + sumNeg;

    if sumPos > sumNeg
      w(k) = wPos;
    else
      w(k) = wNeg;
    end

    % Update p

    p = p + w(k)*R(k,:)';
  end

  v = R\w;
  z = v/norm(v);
end

%--------------------------------------------------------------------------
function [z] = approxNull1(R)
%        [z] = approxNull1(R)
%               returns an approximate null vector of A based on the
%               LINPACK algorithm.  The steps of this algorithm are the
%               following:
%       1.  Solve R'w = e where e has components 1 or -1.  We
%           choose e to make w maximal the sense described below.
%       2.  Solve Rv=w
%       3.  z = v/||v||, and this is the approximate null vector.
%--------------------------------------------------------------------------
  [n,m] = size(R);
  w     = zeros(n,1); e = ones(n,1); I  = diag(e);
  Rd    = diag(R);
% -------------------------------------------------------------------------
% The basic process is the following:
%   1. Preallocate a zero vector w of the appropriate dimension.
%   2. At each iteration k, compute the two possibilities for the next
%      entry of the vector w.  Base the choice on which entry gives you the
%      maximal one norm upon multiplying R'*w^, where w^ =
%      [w(1:k);w_new;zeros].  The algorithm is choosing entries dynamically
%      so that we have maximal norm at each step.
% -------------------------------------------------------------------------
  for k = 1:n
    Rw     = R'*w;
    wPos   =  (1 - Rw(k))/Rd(k);
    wNeg   = -(1 + Rw(k))/Rd(k);
    wPos   = w + wPos*I(:,k);    wNeg = w + wNeg*I(:,k);

    RwPos  = R'*wPos;           RwNeg = R'*wNeg;
    accPos = RwPos(k+1:n);     accNeg = RwNeg(k+1:n);

    sumPos = norm(accPos,1) + abs(wPos(k));
    sumNeg = norm(accNeg,1) + abs(wNeg(k));

    if sumPos > sumNeg
      w = wPos;
    else
      w = wNeg;
    end
  end
  v = R\w;
  z = v/norm(v);
end
