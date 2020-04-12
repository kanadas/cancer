t0 = pi;
x0 = pi;
T = 6*pi;
N = 100;
f = @(t,x)x/t + t*cos(t);
x = euler(f, t0, x0, T, N);
x2 = rk2(f, t0, x0, T, N);
x4 = rk4(f, t0, x0, T, N);
xa = T*(1+sin(T));
printf("Blad:\n euler: %e\n rk2: %e\n rk4: %e\n", sum(abs(x - xa)), sum(abs(x2 - xa)), sum(abs(x4 - xa)));
