function x = rk2(f, t0, x0, T, N)
  t = t0;
  x(1) = x0;
  h = (T - t0) / N;
  for i = 1 : N - 1
    K1 = f(t, x(i));
    K2 = f(t + h, x(i) + h*K1);
    x(i + 1) = x(i) + h/2*(K1 + K2);
    t = t + h;
  end
end