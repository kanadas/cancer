function x = adamus_bashforth(f, t0, x0, x1, T, N)
  h = (T - t0) / N;
  t = t0;
  x(1) = x0;
  x(2) = x1;
  for i = 2 : N - 1
    t += h;
    x(i + 1) = x(i) + h/2*(3*f(t,x(i)) - f(t - h, x(i - 1)));
  end
end