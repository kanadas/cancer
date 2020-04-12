%% Euler scheme
% solves x' = f(t, x), x(t0) = x0
% on [t0, T]
% using N steps
function x = euler(f, t0, x0, T, N)
  h = (T - t0)/N;
  t = t0;
  x(1) = x0;
  for i = 1 : N - 1
    x(i + 1) = x(i) + h*f(t, x(i));
    t = t + h;
  endfor
endfunction