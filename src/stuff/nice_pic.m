t0 = 1.3;
x0 = -1.889;
T = 10;
N = 10000;
t = linspace(t0, T, N);
f = @(t,x)t*x + exp(-x);
#x = euler(f, t0, x0, T, N);
#x2 = rk2(f, t0, x0, T, N);
x4 = rk4(f, t0, x0, T, N);
xm = midpoint(f, t0, x0, x4(2), T, N);
xab = adamus_bashforth(f, t0, x0, x4(2), T, N);
plot(t, x4, "y");
hold on;
#plot(t, x, "g");
#plot(t, x2, "r");
plot(t, xm, "g");
plot(t, xab, "r");
ezplot(@(t,x) t.*x + exp(-x), [1.3, 10, -4, 4]);