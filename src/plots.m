x = t0:1:T;

v = control_to_res(zeros(N,1), points, x0, (T - t0) / 2000);
plt = plot(x, v(x));
leg = legend(plt, {'V1', 'V2', 'K'});
legend boxoff;
set(leg, "fontsize", 15);
print plt "-djpg" "-color" "-S800,800" "plot_zero.jpg";

v = control_to_res(3*ones(N,1), points, x0, (T - t0) / 2000);
control = @(t) ster(t, 3*ones(N,1), points)
subplot(2, 1, 1);
plot(x, control(x));
axis([0, 200, 0, 3.5]);
leg = legend('control');
legend boxoff;
set(leg, "fontsize", 15);
subplot(2, 1, 2);
plt = plot(x, v(x));
leg = legend(plt, {'V1', 'V2', 'K'});
legend boxoff;
set(leg, "fontsize", 15);
print plt "-djpg" "-color" "-S800,800" "plot_max.jpg";

load "sol_150.mat"
load "sol_bang.mat"

v = control_to_res(p, points, x0, (T - t0) / 2000);
control = @(t) ster(t, p, points)
subplot(2, 1, 1);
plot(x, control(x));
axis([0, 200, 0, 3.5]);
leg = legend('control');
legend boxoff;
set(leg, "fontsize", 15);
subplot(2, 1, 2);
plot(x, v(x));
leg = legend({'V1', 'V2', 'K'});
legend boxoff;
set(leg, "fontsize", 15);

print plt "-djpg" "-color" "-S800,800" "plot_weird.jpg";
print plt "-djpg" "-color" "-S800,800" "plot_bang.jpg";

