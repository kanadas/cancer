t0 = 0;
x0 = 0;
T = 10;
N = 1;
x = euler(@(t,x)1, t0, x0, T, N);
xa = T;
printf("Blad: %e\n", abs(x - xa));
