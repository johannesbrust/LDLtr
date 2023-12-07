function [ R, DI, GI ] = bfgsRDRUpdate_E1( R, DI, GI, sk, yk, Hy, nsk, nyk, sy )
% bfgsRDRUpdate_E Updating the inverse LDLT factorization of the BFGS matrix
% and its inverse using eigenvalues of the rank-2 update
%   bfgsRDRUpdate( R, DI, GI, sk, yk, Hy, nsk, nyk, sy )
%
%   Updates the inverse of a LDL' factorization by a rank-2 matrix
%
%   L1 D1 L1' = L D L' - Bs*(Bs/sk'*Bs)' + yk*(yk/sk'*yk)'
%
%   where R = inv(L1), and DI = inv(D1)
%
%   On output D and R are overwritten by DI and R1.
%   The vector GI contains the column norms of R
%
%   INPUTS:
%   R : Upper triangular
%   DI: Diagonal (stored as vector)
%   GI : Column lengths of R (stored as vector)
%   sk : Vector for the BFGS update
%   yk : Vector for the BFGS update
%   Hy : Vector for the inverse BFGS update
%   nsk : norm(sk) (scalar)
%   nyk : norm(yk) (scalar)
%   sy : sk'*yk (scalar)
%
%   OUTPUTS:
%   R : Updated upper triangular
%   DI: Updated diagonal (stored as vector)
%   GI : Diagonals of R'*R
%
%--------------------------------------------------------------------------
% 01/26/23, J.B., Initial implementation
% 03/02/23, J.B., Inverse updating
% 03/10/23, J.B., Eigenvalue update
% 05/23/23, J.B., Storing of diagonals

n   = length(sk);
rr  = size(GI,2);

% Based on eigenvalues
nHy = norm(Hy);
w   = yk/nyk;
w1  = sk/nsk;
nsy = (w1'*w);
a1  = (nsk*(nsy/nyk)+w'*(Hy/nyk))/nsy^2;
a2  = -(nHy)/(nsy*nyk);

%a1 = (sy +yk'*Hy)/sy^2;
%a2 = -1/sy;

% Eigenvalues and eigenvectors
rt   = sqrt(a1^2+4*a2^2);
lam1 = 0.5*(a1+rt);
lam2 = -a2^2/lam1;%0.5*(a1-rt);

v1   = [1;a2/lam1];
nv1  = norm(v1);
v1   = v1/nv1;
%lam1 = lam1*nv1;

%lam1a2 = a2/lam1;
v2      = [-v1(2); v1(1)];
%v2     = [0;(1/sqrt(1+(lam1a2)^2))];
%v2(1)  = -lam1a2*v2(2);

% Computing errors
% BD = R*diag(DI)*R';
% B  = BD;

% Inverse quantities
% w   = yk/nyk;
% nsy = nsk*((sk/nsk)'*w);
% gam = 1 + ((w)'*Hy)/nsy;
% wI   = sk - Hy./gam;
% nwI  = norm(wI);
% wI   = wI/nwI;
% alpI = (gam*nwI*nwI)/sy;

wI   = w1.*v1(1) + (Hy./nHy).*v1(2);
%nwI  = norm(wI);
%wI   = wI./nwI;
alpI = lam1; %*nwI*nwI;

% Error comput
%BD = BD + (alpI.*wI)*wI';

% Algorithm C1
% For updating a symmetric pd matrix

% First update
for j = 1:n
    
    % Inverse updates    
    jI = n-j+1;
    
    pjI      = wI(jI,1);
    djI      = DI(jI,1) + alpI*pjI^2;
    betI     = (pjI*alpI)/djI;
    alpI     = (DI(jI,1)*alpI)/djI;
    
    DI(jI,1)  = djI;
    
    wI(1:(n-j),1) = wI(1:(n-j),1) - pjI*R(1:(n-j),jI);
    R(1:(n-j),jI) = R(1:(n-j),jI) + betI*wI(1:(n-j),1);
   
end

% Inverse quantities
% wI   = Hy;
% nwI  = norm(wI);
% wI   = wI/nwI;
% alpI = -(nwI*nwI)/(sy*gam);

wI   = w1.*v2(1) + (Hy./nHy).*v2(2);
% nwI  = norm(wI);
% wI   = wI./nwI;
alpI = lam2; % *nwI*nwI;

%alpI = lam2;

% Error comput
%BD = BD + (alpI.*wI)*wI';

for j = 1:n
        
    % Inverse updates    
    jI = n-j+1;
    
    pjI      = wI(jI,1);
    djI      = DI(jI,1) + alpI*pjI^2;
    betI     = (pjI*alpI)/djI;
    alpI     = (DI(jI,1)*alpI)/djI;
    
    DI(jI,1)  = djI;
    
    wI(1:(n-j),1) = wI(1:(n-j),1) - pjI*R(1:(n-j),jI);
    R(1:(n-j),jI) = R(1:(n-j),jI) + betI*wI(1:(n-j),1);
    
    cGI = min(j,rr); %max(rr - j + 1, 1);
    cGE = min(jI + rr -1, n);
    
    GI(jI,1:cGI) = R(1:jI,jI)'*R(1:jI,jI:cGE);
    
    %GI(jI)        = norm(R(1:jI,jI));
    
end

% Make matrix positive definite
[md, idx] = min(DI);
if md < 0; DI(idx) = abs(DI(idx)); end;

% Error comput
% V = [w1, Hy./nHy]*[v1,v2];
% V2  = [sk Hy]*[v1,v2];
% Lam = diag([lam1,lam2]);
% %ER = norm(BD - B + V2*Lam*V2','fro');
% ER1 = norm(BD - B + V*Lam*V','fro');
% 
% aa1 = (sy +yk'*Hy)/sy^2;
% aa2 = -1/sy;
% A2  = [aa1 aa2; aa2 0];
% RK2 = V2*A2*V2';
% 
% DD = diag([nsk, nHy]);
% 
% ER2 = norm(RK2-V*Lam*V','fro');

%er = norm(R*(DI.*(R'*yk))-sk);


end

