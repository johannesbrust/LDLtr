function [bet_, s] = r1(w, Djj)
   n    = length(w);
   k    = find(abs(w) > eps*norm(w), 1, 'last');

   bet_ = zeros(n, 1);
   s    = zeros(n, 1);

   bet_(k) = 1/w(k);

   for j = k-1:-1:1
      rj        = bet_(j + 1)*w(j);
      s(j)      = 1/sqrt(rj^2 + Djj(j)/Djj(j + 1));

      bet_(j)   =    s(j)*bet_(j + 1);
      bet_(j+1) = rj*s(j)*bet_(j + 1);
   end
