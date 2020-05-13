pkg load optim

gmax = 3;
t0 = 0;
T = 200;
x0 = [20; 280; 650];

step = 4;
points = t0 : step : T;
N = (T - t0) / step;
disp(N);

load "solution0.mat";

#[p0, objf0, cvg, outp] = nonlin_min(@(g) objf(g, points, x0), zeros(N, 1), optimset("lbound", zeros(N, 1), "ubound", gmax * ones(N, 1)));

[p15, objf15, cvg, outp] = nonlin_min(@(g) objf(g, points, x0), 1.5 * ones(N, 1), optimset("lbound", zeros(N, 1), "ubound", gmax * ones(N, 1)));

[pmax, objfmax, cvg, outp] = nonlin_min(@(g) objf(g, points, x0), gmax * ones(N, 1), optimset("lbound", zeros(N, 1), "ubound", gmax * ones(N, 1)));

#{
starting_points = {zeros(N, 1),
		   ones(N,1),
		   2 * ones(N, 1),
		   gmax * ones(N, 1)}

for i = 1:length(starting_points)
  [p, objf, cvg, outp] = nonlin_min(@J, starting_points{i}, optimset("lbound", zeros(N, 1), "ubound", gmax * ones(N, 1)));
  printf("%d %f\n", i, objf)
endfor
