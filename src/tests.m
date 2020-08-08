
t0 = 0;
T = 200;
x0 = [20; 280; 650];
STEP = 0.1;
h = STEP;
points = t0 : 4 : T;
N = length(points)

grad = dfxpdp(1, zeros(N,1), @(g,t) rk4(@(t, V) dynamics(t, V, @(t) ster(t, g, points)), t0, t, x0, STEP));
[x, dx] = rk4(@(t, V, dV) dynamics(t, V, @(t) ster(t, zeros(N,1), points), dV), t0, 1, x0, STEP, N);
dx = squeeze(dx);

f = @(t, V, dV) dynamics(t, V, @(t) ster(t, zeros(N,1), points), dV);
[x, dfdx, dfdg] = f(0, x0, [0;0;0]);

grad_x = dfxpdp(0, x0, @(p,t) f(t, p, [0;0;0]));
grad_g = dfxpdp(0, zeros(N,1), @(p,t) dynamics(t, x0, @(t) ster(t, p, points)));

function [x,dx, dk1,dk2,dk3,dk4] = test_rk(f, t, x0, h, grad_dim=0)
  x = x0;
  dx = zeros(length(x), grad_dim);
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
    [k3, dfdx, dfdg] = f(t + h/2, x + h/2*k2, dx + h/2*dk2);
    dk3 = dfdx * (dx + h/2*dk2) + dfdg;
    [k4, dfdx, dfdg] = f(t + h, x + h*k3, dx + h*dk3);
    dk4 = dfdx * (dx + h*dk3) + dfdg;

    dx = dx + h/6*(dk1 + 2*dk2 + 2*dk3 + dk4);
  endif
  x = x + h/6*(k1 + 2*k2 + 2*k3 + k4);
endfunction

grad = dfxpdp(0, zeros(N,1), @(p,t) test_rk(@(t, V) dynamics(t, V, @(t) ster(t, p, points)), t, x0, STEP));
[x,dx,dk1,dk2,dk3,dk4] = test_rk(f, 0, x0, STEP, N);

function k2 = test_k2(f, t, x0, h)
  x = x0;
  k1 = f(t, x);
  k2 = f(t + h/2, x + h/2*k1);
endfunction

grad_k2 = dfxpdp(0, zeros(N,1), @(p,t) test_k2(@(t, V) dynamics(t, V, @(t) ster(t, p, points)), t, x0, STEP));

dx = zeros(3, N);
[k1, dfdx, dfdg] = f(0, x0, [0;0;0]);
dk1 = dfdx * dx + dfdg;
[k2, dfdx, dfdg] = f(h/2, h/2*k1, dx + h/2*dk1);
dk2 = dfdx * (dx + h/2*dk1) + dfdg;
