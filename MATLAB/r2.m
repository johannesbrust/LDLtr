%--------------------------------------------------------------------------
function [lam, bet, gam] = r2(bet_, s, w, z, Djj)
   n    = length(w);
   k    = find(abs(w) > eps*norm(w), 1, 'last');

   lam  =  ones(n, 1);
   bet  = zeros(n, 1);
   gam  = zeros(n, 1);

   gam_ = 1/(Djj(1)*bet_(1));
   lam_ = bet_(1)*w(1) + gam_(1)*z(1);

   for j = 1:n-1
      q           = Djj(j)/Djj(j + 1);
      lam(j)      = sqrt(lam_^2 + q*s(j)^2);
      c_          = lam_/lam(j);             s_ = s(j)/lam(j);
      bet(j)      = c_*bet_(j) - s_*bet_(j+1);
      gam(j)      = c_*gam_;
      gam_        = q*s_*gam_;
      bet_(j + 1) = q*s_*bet_(j) + c_*bet_(j + 1);
      lam_        = bet_(j + 1)*w(j + 1) + gam_*z(j + 1);
      if j == k
         return;
      end
   end
   lam(n) = lam_;
