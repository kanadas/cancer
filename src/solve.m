pkg load optim

gmax = 3;
t0 = 0;
T = 200;
x0 = [20; 280; 650];

step = 4;
points = t0 : step : T;
N = (T - t0) / step + 1;
disp(N);

load "sol0.mat";

[p0, objf0, cvg, outp] = fmincon(OBJF = @(g) objf(g, points, x0),
				 X0   = zeros(N, 1),
				 A    = [],
				 B    = [],
				 AEQ  = [],
				 BEQ  = [],
				 LB   = zeros(N, 1),
				 UB   = gmax * ones(N, 1),
				 NONLCON = [],
				 OPTIONS=optimset("GradObj", "on"));

save "sol0.mat" p0;

step = 1;
points = t0 : step : T;
N = (T - t0) / step + 1;
disp(N);

[p, res, cvg, outp] = fmincon(OBJF = @(g) objf(g, points, x0),
			       X0   = zeros(N, 1),
			       A    = [],
			       B    = [],
			       AEQ  = [],
			       BEQ  = [],
			       LB   = zeros(N, 1),
			       UB   = gmax * ones(N, 1),
			       NONLCON = [],
			       OPTIONS=optimset("GradObj", "on"));

[p_max, res_max, cvg, outp] = fmincon(OBJF = @(g) objf(g, points, x0),
			       X0   = gmax * ones(N, 1),
			       A    = [],
			       B    = [],
			       AEQ  = [],
			       BEQ  = [],
			       LB   = zeros(N, 1),
			       UB   = gmax * ones(N, 1),
			       NONLCON = [],
			       OPTIONS=optimset("GradObj", "on"));

[p_mid, res_mid, cvg, outp] = fmincon(OBJF = @(g) objf(g, points, x0),
			       X0   = gmax/2 * ones(N, 1),
			       A    = [],
			       B    = [],
			       AEQ  = [],
			       BEQ  = [],
			       LB   = zeros(N, 1),
			       UB   = gmax * ones(N, 1),
			       NONLCON = [],
			       OPTIONS=optimset("GradObj", "on"));

%load "old_sol0.mat";
%[old_p0, old_objf0, cvg, outp] = nonlin_min(@(g) objf(g, points, x0), zeros(N, 1), optimset("lbound", zeros(N, 1), "ubound", gmax * ones(N, 1)));

%[p15, objf15, cvg, outp] = nonlin_min(@(g) objf(g, points, x0), 1.5 * ones(N, 1), optimset("lbound", zeros(N, 1), "ubound", gmax * ones(N, 1)));

%[pmax, objfmax, cvg, outp] = nonlin_min(@(g) objf(g, points, x0), gmax * ones(N, 1), optimset("lbound", zeros(N, 1), "ubound", gmax * ones(N, 1)));

#{
starting_points = {zeros(N, 1),
		   ones(N,1),
		   2 * ones(N, 1),
		   gmax * ones(N, 1)}

for i = 1:length(starting_points)
  [p, objf, cvg, outp] = nonlin_min(@J, starting_points{i}, optimset("lbound", zeros(N, 1), "ubound", gmax * ones(N, 1)));
  printf("%d %f\n", i, objf)
endfor
