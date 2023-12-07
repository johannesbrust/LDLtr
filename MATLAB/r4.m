function[R] = r4(R, alph, q, bet_, p, lam)
n = length(bet_);
w = zeros(n, 1);
z = zeros(n, 1);

for i = n:-1:1
   w(i)    =   p(i)*R(i, i);
   z(i)    =   q(i)*R(i, i);
   R(i, i) = lam(i)*R(i, i);
   for j = i+1:n
      r = R(i,j);
      R(i , j) = lam(i)*r + alph(i)*z(j) + bet_(i)*w(j);
      w(j)     = w(j) + p(i)*r;
      z(j)     = z(j) + q(i)*r;
   end
end