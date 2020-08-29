pkg load optim
gmax = 3;
t0 = 0;
T = 200;
x0 = [20; 280; 650];
h = 0.1;

constants.lambda1 = 0.192;
constants.lambda2 = 0.192;
constants.mu = 0.0;
constants.b1 = 5.85;
constants.b2 = 5.85;
constants.d = 0.00873;
constants.beta1 = 0.15;
constants.beta2 = 0.1;
constants.beta = 0.05;
%constants.alpha12 = 0.1;
%constants.alpha21 = 0.15;
constants.alpha12 = 0.5;
constants.alpha21 = 0.75;
constants.epsilon = 0.01;
%constants.omega = 1000;
constants.omega = 2000;

constants1 = constants;
constants1.alpha12 = 0.1;
constants1.alpha21 = 0.15;
constants1.omega = 1000;

points = t0 : 4 : T;
N = length(points);
discr_fun = @const_discr;

start = ones(N,1);

test_grad(start, @(g) objf(g, points, x0, h, discr_fun, constants));

test_grad(x0, @(v) dynamics(0, v, @(t) discr_fun(t, start, points), constants));

function [res, dg] = dynamics_dg(t, V, g, constants)
  if nargout == 1
    res = dynamics(t, V, g, constants);
  else
    [res, dv, dg] = dynamics(t, V, g, constants);
  endif
endfunction
test_grad(start, @(g) dynamics_dg(0, x0, @(t) discr_fun(t, g, points), constants));

function [v1, dv1] = test_v1(t, V, g, constants)
  lambda1 = constants.lambda1;
  beta1 = constants.beta1;
  alpha12 = constants.alpha12;
  if nargout == 1
    ggt = g(t);
  else
    [ggt, dgt] = g(t);
  endif
  function y = F(x)
    y = -log(x);
  endfunction
  function y = dF(x)
    y = -1 ./ x;
  endfunction

  V1 = V(1);
  V2 = V(2);
  K = V(3);
  frac1 = (V1 + alpha12 * V2) ./ K;
  v1 = lambda1 * V1 .* F(frac1) - beta1 * V1 .* ggt;
  if nargout > 1
    dv1dv1 = lambda1 * (F(frac1) + V1 .* dF(frac1) ./ K) - beta1 .* ggt;
    dv1dv2 = lambda1 * V1 .* dF(frac1) .* (alpha12 / K);
    dv1dk = lambda1 * V1 .* dF(frac1) .* ( -1 * frac1 ./ K);
    dv1 = [dv1dv1, dv1dv2, dv1dk];
  endif
endfunction
test_grad(x0, @(V) test_v1(0, V, @(t) discr_fun(t, start, points), constants));

test_grad(start, @(g) integrate(@(x) integral_func(x, control_to_res(g, points, x0, h, discr_fun, constants), constants), t0, T, h)); %BAD

test_grad(start, @(g) integral_func(50, control_to_res(g, points, x0, h, discr_fun, constants), constants)); %BAD

test_grad(start, @(g) control_to_res(g, points, x0, h, discr_fun, constants)(50)); %BAD

test_grad(start, @(g) discr_fun(50, g, points)); %GOOD

test_grad(start, @(g) rk4(@(t, V, dV = 0) dynamics(t, V, @(t) discr_fun(t, g, points), constants, dV), t0, 50, x0, h, length(g))); %BAD

function [k1, dk1] =  test_k1(t, x, dx, f)
  if nargout == 1
    k1 = f(t, x);
  else
    [k1, dfdx, dfdg] = f(t, x, dx);
    dk1 = dfdx * dx + dfdg;
  endif
endfunction

%x50 = [186.68; 392.93; 881.77];
%dx50 = [-53.02656; -72.45939; -84.42554]; %CHYBA
test_grad(start, @(g) test_k1(0, x0, zeros(length(x0), length(g)), @(t, V, dV = 0) dynamics(t, V, @(t) discr_fun(t, g, points), constants, dV)));

function [k2, dk2] = test_k4(t, x, dx, h, f)
  if nargout == 1
    k1 = f(t, x);
    k2 = f(t + h/2, x + h/2*k1);
    k3 = f(t + h/2, x + h/2*k2);
    k4 = f(t + h, x + h*k3);
  else
    [k1, dfdx, dfdg] = f(t, x, dx);
    dk1 = dfdx * dx + dfdg;
    [k2, dfdx, dfdg] = f(t + h/2, x + h/2*k1, dx + h/2*dk1);
    dk2 = dfdx * (dx + h/2*dk1) + dfdg;
%    [k3, dfdx, dfdg] = f(t + h/2, x + h/2*k2, dx + h/2*dk2);
%    dk3 = dfdx * (dx n+ h/2*dk2) + dfdg;
%    [k4, dfdx, dfdg] = f(t + h, x + h*k3, dx + h*dk3);
%    dk4 = dfdx * (dx + h*dk3) + dfdg;
%
%    dx = dx + h/6*(dk1 + 2*dk2 + 2*dk3 + dk4);
%
%    disp(dx);
  endif
endfunction
test_grad(start, @(g) test_k4(0, x0, zeros(length(x0), length(g)), h, @(t, V, dV = 0) dynamics(t, V, @(t) discr_fun(t, g, points), constants, dV))); %BAD

function [v, dv] = test_dv(t, x, dx, h, f)
  if nargout == 1
    k1 = f(t, x);
    v = x + h/2*k1;
  else
    [k1, dfdx, dfdg] = f(t, x, dx);
    dk1 = dfdx * dx + dfdg;
    v = x + h/2*k1;
    dv = dx + h/2*dk1;
  endif
endfunction
test_grad(start, @(g) test_dv(0, x0, zeros(length(x0), length(g)), h, @(t, V, dV = 0) dynamics(t, V, @(t) discr_fun(t, g, points), constants, dV)));

function [v, dg] = test_dyn_dg(t, V, h, g, constants, dV)
  lambda1 = constants.lambda1;
  lambda2 = constants.lambda2;
  mu = constants.mu;
  b1 = constants.b1;
  b2 = constants.b2;
  d = constants.d;
  beta1 = constants.beta1;
  beta2 = constants.beta2;
  beta = constants.beta;
  alpha12 = constants.alpha12;
  alpha21 = constants.alpha21;
  function y = F(x)
    y = -log(x);
  endfunction

  function y = dF(x)
    y = -1 ./ x;
  endfunction

  if nargout == 1
    ggt = g(t);
  else
    [ggt, dgt] = g(t);
  endif

  [k1, dfdx, dfdg] = dynamics(t, V, g, constants, dV);
  dk1 = dfdx * dV + dfdg;
  t = t + h/2;
  V = V + h/2*k1;
  dV = dV + h/2*dk1;
  
  V1 = V(1);
  V2 = V(2);
  K = V(3);
  frac1 = (V1 + alpha12 * V2) ./ K;
  frac2 = (V2 + alpha21 * V1) ./ K;
  pow = (V1 + V2).^(2/3);
  v1 = lambda1 * V1 .* F(frac1) - beta1 * V1 .* ggt;
  v2 = lambda2 * V2 .* F(frac2) - beta2 * V2 .* ggt;
  k = -mu * K + (V1 + V2) - d * pow .* K - beta * K .* ggt;
  v = [v1; v2; k];

  if nargout > 1 %gradient required
    dpow = 2/3 * pow ./ (V1 + V2);
    dV1 = dV(1,:);
    dV2 = dV(2,:);
    dK = dV(3,:);
    dfrac1 = ((dV1 + alpha12 * dV2) .* K - (V1 + alpha12 * V2) .* dK) ./ (K .* K);
    dfrac2 = ((dV2 + alpha21 * dV1) .* K - (V2 + alpha21 * V1) .* dK) ./ (K .* K);
    dpowdv = dpow .* (dV1 + dV2);
    dv1dg = lambda1 * (dV1 .* F(frac1) + V1 .* dF(frac1) .* dfrac1) - beta1 * (dV1 .* ggt + V1 .* dgt);    
    dv2dg = lambda2 * (dV2 .* F(frac2) + V2 .* dF(frac2) .* dfrac2) - beta2 * (dV2 .* ggt + V2 .* dgt);
    dkdg = -mu * dK + dV1 + dV2 - d * (dpowdv .* K + pow .* dK) - beta * (dK .* ggt + K .* dgt);
    dg = [dv1dg; dv2dg; dkdg];
  endif
endfunction

function [v, dv] = test_dyn_dv(t, V, h, g, constants, dV)
  lambda1 = constants.lambda1;
  lambda2 = constants.lambda2;
  mu = constants.mu;
  b1 = constants.b1;
  b2 = constants.b2;
  d = constants.d;
  beta1 = constants.beta1;
  beta2 = constants.beta2;
  beta = constants.beta;
  alpha12 = constants.alpha12;
  alpha21 = constants.alpha21;
  function y = F(x)
    y = -log(x);
  endfunction

  function y = dF(x)
    y = -1 ./ x;
  endfunction

  if nargout == 1
    ggt = g(t);
  else
    [ggt, dgt] = g(t);
  endif

  [k1, dfdx, dfdg] = dynamics(t, V, g, constants, dV);
  dk1 = dfdx * dV + dfdg;
  t = t + h/2;
  V = V + h/2*k1;
  dV = dV + h/2*dk1;
  
  V1 = V(1);
  V2 = V(2);
  K = V(3);
  frac1 = (V1 + alpha12 * V2) ./ K;
  frac2 = (V2 + alpha21 * V1) ./ K;
  pow = (V1 + V2).^(2/3);
  v1 = lambda1 * V1 .* F(frac1) - beta1 * V1 .* ggt;
  v2 = lambda2 * V2 .* F(frac2) - beta2 * V2 .* ggt;
  k = -mu * K + (V1 + V2) - d * pow .* K - beta * K .* ggt;
  v = [v1; v2; k];

  if nargout > 1 %gradient required
    dpow = 2/3 * pow ./ (V1 + V2);
    dv1dv1 = lambda1 * (F(frac1) + V1 .* dF(frac1) ./ K) - beta1 .* ggt;
    dv1dv2 = lambda1 * V1 .* dF(frac1) .* (alpha12 / K);
    dv1dk = lambda1 * V1 .* dF(frac1) .* ( -1 * frac1 ./ K);
    dv2dv1 = lambda2 * V2 .* dF(frac2) .* (alpha21 / K);
    dv2dv2 = lambda2 * (F(frac2) + V2 .* dF(frac2) ./ K) - beta2 .* ggt;
    dv2dk = lambda2 * V2 .* dF(frac2) .* ( -1 * frac2 ./ K);
    dkdv1 = 1 - d * dpow .* K;
    dkdv2 = 1 - d * dpow .* K;
    dkdk = -mu - d * pow - beta * ggt;
    dv = [dv1dv1, dv1dv2, dv1dk; dv2dv1, dv2dv2, dv2dk; dkdv1, dkdv2, dkdk];
  endif
endfunction

test_grad(x0, @(V) test_dyn_dv(0, V, h, @(t) discr_fun(t, start, points), constants, zeros(length(x0), length(start))));

load "bug_log"; %g dy test_grad
clear test_grad;
test_grad(g, @(g) objf(g, points, x0, h, discr_fun, constants)); %BAD

test_grad(x0, @(v) dynamics(0, v, @(t) discr_fun(t, start, points), constants)); %GOOD

test_grad(g, @(g) dynamics_dg(0, x0, @(t) discr_fun(t, g, points), constants)); %GOOD

function [v1, dv1] = test_v1(t, V, g, constants)
  lambda1 = constants.lambda1;
  beta1 = constants.beta1;
  alpha12 = constants.alpha12;
  if nargout == 1
    ggt = g(t);
  else
    [ggt, dgt] = g(t);
  endif
  function y = F(x)
    y = -log(x);
  endfunction
  function y = dF(x)
    y = -1 ./ x;
  endfunction

  V1 = V(1);
  V2 = V(2);
  K = V(3);
  frac1 = (V1 + alpha12 * V2) ./ K;
  v1 = lambda1 * V1 .* F(frac1) - beta1 * V1 .* ggt;
  if nargout > 1
    dv1dv1 = lambda1 * (F(frac1) + V1 .* dF(frac1) ./ K) - beta1 .* ggt;
    dv1dv2 = lambda1 * V1 .* dF(frac1) .* (alpha12 / K);
    dv1dk = lambda1 * V1 .* dF(frac1) .* ( -1 * frac1 ./ K);
    dv1 = [dv1dv1, dv1dv2, dv1dk];
  endif
endfunction
test_grad(x0, @(V) test_v1(0, V, @(t) discr_fun(t, start, points), constants)); %GOOD

test_grad(g, @(g) integrate(@(x) integral_func(x, control_to_res(g, points, x0, h, discr_fun, constants), constants), t0, T, h)); %BAD

test_grad(g, @(g) integral_func(1, control_to_res(g, points, x0, h, discr_fun, constants), constants)); %GOOD

test_grad(g, @(g) control_to_res(g, points, x0, h, discr_fun, constants)(1)); %GOOD

test_grad(g, @(g) discr_fun(1, g, points)); %GOOD

test_grad(g, @(g) rk4(@(t, V) dynamics(t, V, @(t) discr_fun(t, g, points), constants), t0, 1, x0, h, length(g))); %GOOD

test_grad(g, @(g) integral_func(50, control_to_res(g, points, x0, h, discr_fun, constants), constants));

function [y,dy] = test_integrate(f, t0, T, h)
  eval_points = t0:h:T;
  if nargout == 1
    y = f(eval_points);
  else
    [y, dy] = f(eval_points);
    dy = dy';
  endif
endfunction
test_grad(g, @(g) test_integrate(@(x) integral_func(x, control_to_res(g, points, x0, h, discr_fun, constants), constants), t0, T, h)); %BAD
%result
% 23.791
%max difference
% 85.294
%on position
% 676
test_grad(g, @(g) test_integrate(@(x) integral_func(x, control_to_res(g, points, x0, h, discr_fun, constants), constants), 650*h, 700*h, h)); %BAD

test_grad(g, @(g) integral_func(675*h, control_to_res(g, points, x0, h, discr_fun, constants), constants)); %BAD

test_grad(g, @(g) control_to_res(g, points, x0, h, discr_fun, constants)(675*h)); %GOOD

function [x, dx] = test_integral(t, dyn, constants)
  epsilon = constants.epsilon;
  omega = constants.omega;
  function y = G(x)
    y = (1 + tanh(x))/2;
  endfunction

  function y = dG(x)
    y = sech(x) .^ 2 / 2;
  endfunction
  
  if nargout == 1
    V = dyn(t);
  else
    [V, dV] = dyn(t);
  endif
%  x = V(:,1);
  x = G((V(:, 2) - V(:, 1)) / epsilon);
  if nargout > 1
%    dx = dV(:, 1, :);
    disp(squeeze(V))
    disp("")
    disp(squeeze(dV)')
    disp("")
    disp(squeeze(dG((V(:, 2) - V(:, 1)) / epsilon))')
    disp("")
    disp(((dV(:, 2, :) - dV(:, 1, :)) / epsilon))

    dx = dG((V(:, 2) - V(:, 1)) / epsilon) .* ((dV(:, 2, :) - dV(:, 1, :)) / epsilon);
    dx = squeeze(dx)';

    disp("")
    disp(dx')
  endif
endfunction
test_grad(g, @(g) test_integral(675*h, control_to_res(g, points, x0, h, discr_fun, constants), constants));

