function[prob] = getCUSTomProblem()
   prob      = struct;
   prob.n    = 5;
   prob.m    = 0;
   prob.x    = [1; 2; 3; 4; 5];
   prob.name = "EXPERIM";
   prob.obj  = @(x)custom_obj(x);
   global B;
   A = [0.3596    0.1249    0.9569    0.7593    0.4633;
        0.5583    0.0244    0.9357    0.7406    0.2122;
        0.7425    0.2902    0.4579    0.7437    0.0985;
        0.4243    0.3175    0.2405    0.1059    0.8236;
        0.4294    0.6537    0.7639    0.6816    0.1750];
   B = fix(100 * (1/2*(A + A') + eye(5)));
   prob.B = B;
end

% function[f, g] = custom_obj(x)
%    H = toeplitz([16, 8, 4, 2, 1]);
%    c = ones(5, 1);
%    f = c'*x + 1/2 *x'*H*x;
%    g = c + H*x;
% end

function[f, g] = custom_obj(x)
   global B;
   c = -B*ones(5,1);
   f = c'*x + 1/2 * x'*B*x;
   g = c + B*x;
end

% function[prob] = getCUSTomProblem()
%    prob      = struct;
%    prob.n    = 4;
%    prob.m    = 0;
%    prob.x    = ones(4, 1);
%    prob.name = "POWELL";
%    prob.obj  = @(x)custom_obj(x);
% end

% function[f, g] = custom_obj(x)
%    G = [0.1, 0.1, 0.1, 0.1;
%         0.1, 0.2, 0.1, 0.1;
%         0.1, 0.1, 0.3, 0.1;
%         0.1, 0.1, 0.1, 0.4];
%    f = 1/2 *x'*G*x;
%    g = G*x;
% end