%% Rungge-Kutta 4th order
% x' = f(t, x), x(t0) = x0
% x \in [t0, T]
% N steps
function x = rk4(f, t0, x0, T, N)
  t = t0;
  x = x0;
  h = (T - t0) ./ N;
  for i = 1 : N - 1
    K1 = f(t, x);
    K2 = f(t + h/2, x + h/2*K1);
    K3 = f(t + h/2, x + h/2*K2);
    K4 = f(t + h, x + h*K3);
    x = x + h/6*(K1 + 2*K2 + 2*K3 + K4);
    t = t + h;
  end
end
