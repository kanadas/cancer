t0 = 0;
x0 = 1;
T = 4;
N = 200;
t = linspace(t0, T, N);
f = @(t,x)-2*x;
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
plot(t, exp(-2*t), "k");
#ezplot(@(t,x) t.*x + exp(-x), [1.3, 10, -4, 4]);
