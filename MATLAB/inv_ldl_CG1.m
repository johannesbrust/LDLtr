function [ x ] = inv_ldl_CG1( sig, g, R, D, tol, maxit )
%inv_ldl_CG1 Approximate solve with a LDL factorization
%
% The system to solve is
%
%       (sig.inv(L)inv(L')+D)x = g
%
%-----------------------------------------------------------------%
% 05/09/2023, J.B., Initial version

x       = zeros(size(g,1),1);
r       = g;
nr0     = norm(r); %max(norm(r),1);
p       = r;
for icg = 1:maxit

    %Ap      = R'*(R*p)+D.*(p./sig);
    Ap      = sig.*(R'*(R*p)) + D.*p;
    %sig.*p +  linsolve(R,(D.*linsolve(R,p,optsUT)),optsUTT);
    %Ap      = aApply(p);
    rr      = r'*r;
    alpc    = rr/(p'*Ap);
    x       = x + alpc.*p;
    r       = r - alpc.*Ap;
    nr      = norm(r);
    if nr/nr0 < tol
        break;
    end
    betc    = (r'*r)/rr;
    p       = r + betc.*p;

end

end
