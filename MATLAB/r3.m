function[H_] = r3(R, bet_, u, s)
    n = length(bet_);
    H_ = zeros(n);
    w  = R'*u;
    H_(1, :) = bet_(1)*w';
    for i = 2:n
        H_(i, i-1) = s(i-1)*R(i-1, i-1);
        w(i:n) = w(i:n) - u(i-1)*R(i-1,i:n)';
        H_(i, i:n) = s(i-1)*R(i-1,i:n)' + bet_(i)*w(i:n);
    end
end