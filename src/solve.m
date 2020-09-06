pkg load optim

%global DEBUG = 2;

gmax = 3;
t0 = 0;
T = 200;
x0 = [20; 280; 650];
h = 0.1;

constants.lambda1 = 0.192;
constants.lambda2 = 0.192;
constants.mu = 0.0;
constants.b1 = 5.85;
constants.b2 = 5.85;
constants.d = 0.00873;
constants.beta1 = 0.15;
constants.beta2 = 0.1;
constants.beta = 0.05;
%constants.alpha12 = 0.1;
%constants.alpha21 = 0.15;
constants.alpha12 = 0.5;
constants.alpha21 = 0.75;
constants.epsilon = 0.01;
%constants.omega = 1000;
constants.omega = 2000;

constants1 = constants;
constants1.alpha12 = 0.1;
constants1.alpha21 = 0.15;
constants1.omega = 1000;

points = t0 : 4 : T;

points = t0 : 1 : T;

points = t0 : 0.5 : T;

points = [(0 : 0.5 : 50), (51 : 1 : 149), (150 : 0.5 : 200)];

points = [(0: 1 : 24), (25 : 0.1 : 75), (76 : 1 : 200)];

N = length(points);

start_bang = 1:N;
start_bang = start_bang < N/2;

[y0, dy0] = objf(zeros(N,1), points, x0, h, @const_discr, constants);

[p, res, cvg, outp] = run_opt(points, x0, zeros(N,1), h, @const_discr, constants);
[y, dy] = objf(p, points, x0, h, @const_discr, constants);

#no gradient calculation experimment
points = t0 : 4 : T;
N = length(points);
[p_tstng, res_tstng, cvg_tstng, outp_tstng] = run_opt_bak(points, x0, zeros(N,1), h, @const_discr, constants);
[p_ng, res_ng, cvg_ng, outp_ng] = run_opt_bak(points, x0, zeros(N,1), h, @const_discr, constants, "lm_feasible", "off");

save "sol_tstng.mat" p_tstng res_tstng cvg_tstng outp_tstng;
save "sol_nograd.mat" p_ng res_ng outp_ng nograd_err;
load "data/sol_nograd.mat";

[res, grad_tstng] = objf(p_tstng, points, x0, h, @const_discr, constants);
hook.diffp =.001 * ones (size (p_tstng)); %Default
grads = zeros(8, size(p_tstng));
for i = 1:8
  grads(i,:) = dfpdp(p_tstng, @(g) objf(g, points, x0, h, @const_discr, constants), hook);
  diff(i) = norm(grad_tstng - grads(i,:), 1);
  disp(strcat("Diff(", num2str(i), ") = ", num2str(diff(i))));
  hook.diffp = 0.1*hook.diffp;
endfor
save "grads_findif.mat" grads;

for i = 1:8
  diff(i) = norm(grad_tstng - grad(i,:), 1);
endfor

[p, res, cvg, outp] = run_opt(points, x0, ones(N,1), h);
[p_max, res_max, cvg, outp] = run_opt(points, x0, gmax*ones(N,1), h);
[p, res, cvg, outp] = run_opt(points, x0, ((N-1):-1:0)*3/(N-1), h);
[p, res, cvg, outp] = run_opt(points, x0, start_bang, h);

[p_ng, res_ng, cvg, outp] = run_opt(points, x0, zeros(N,1), h, "off");

%assume N = 401
start1 = [(400:-2:0), (0:3:100), (100:-3:0), zeros(1,132)]*3/400;
start2 = [(300:-3:0), (0:4:100), (100:-4:0), zeros(1,248)] / 100;
start3 = [300*ones(1,50), zeros(1,351)]/100;
start4 = [300*ones(1,25), zeros(1,50), 300*ones(1,25), zeros(1,50), 300*ones(1,25) ,zeros(1,226)]/100;
start5 = [300*ones(1,50), 16*ones(1,351)]/100;
start6 = [zeros(1,100), 300*ones(1,10), zeros(1,291)] / 100;
start7 = [300*ones(1,50), zeros(1, 50), 16*ones(1,20), zeros(1,20), 30*ones(1,20), zeros(1,241)] / 100;
start8 = [300*ones(1,50), zeros(1, 50), 16*ones(1,40), 30*ones(1,10), 16*ones(1,20), 30*ones(1,10), 16*ones(1,20), 30*ones(1,10), 16*ones(1,191)] / 100;

%start2 = [300*ones(1,100), (300:-3:0), (0:4:100), (100:-4:0), (0:4:100), (100:-4:0), zeros(1,96)]/100;

[p0, res0, cvg, outp] = run_opt(points, x0, zeros(N,1), h);

start_lin = [0:3:1200]/400;

start = [zeros(1,100), 40*ones(1,301)]/100;
start = [15*ones(1,100), 55*ones(1,301)]/100;
start = [zeros(1,85), 40*ones(1,316)]/100;
start = [zeros(1,85), 55*ones(1,316)]/100;
start = [300*ones(1,20), zeros(1, 80), 50*ones(1,301)]/100;

[p, res, cvg, outp] = run_opt(points, x0, start, h, @const_discr, constants);
[p, res, cvg, outp] = run_opt(points, x0, start, h, @const_discr, constants, "active-set");
[p, res, cvg, outp] = run_opt(points, x0, start_lin, h, "active-set");

save "sol_bang.mat" p;

[p_c1, res_c1, cvg_c1, outp_c1] = run_opt(points, x0, zeros(N,1), h, @const_discr, constants1, "active-set");
save "data/sol_test_c1.mat" p_c1 res_c1 outp_c1;

%Test sqp siatka S_1 start40 (na testach wyszedł absurdalny wynik)
points = t0 : 1 : T;
N = length(points);
start = [zeros(1,42), 40*ones(1,159)]/100;
[p_sqs1, res_sqs1, cvg_sqs1, outp_sqs1] = run_opt(points, x0, start, h, @const_discr, constants, "active-set");
save "data/sol_test_sqs1.mat" p_sqs1 res_sqs1 outp_sqs1;

%Eksperymenty z gęstą siatką
points = t0 : 0.1 : T;
N = length(points);
h = 0.1;
start40 = [zeros(1,85*5), 40*ones(1,315*5+1)]/100;
[p_s1, res_s1, cvg_s1, outp_s1] = run_opt(points, x0, zeros(N,1), h, @const_discr, constants, "active-set");
disp("DONE 1")
[p_s2, res_s2, cvg_s2, outp_s2] = run_opt(points, x0, start40, h, @const_discr, constants, "active-set");

%BEST
% 2. & {\it sqp\/} & stała & $S_{0.5}$ & 0.1 & $g_{0,0.55}$ & 2.71 & 10 & 130 \\
points = t0 : 0.5 : T;
N = length(points);
start = [zeros(1,85), 55*ones(1,316)]/100;
[p_min, res_min, cvg_min, outp_min] = run_opt(points, x0, start, h, @const_discr, constants, "active-set");
save "data/sol_min.mat" p_min, res_min, outp_min;

p_best = p0;
res_best = res0;
noncont_best = 400;
for i = 50:5:150
  disp(strcat("Compting i=", num2str(i)));
  start = [zeros(1,i), 40*ones(1,401-i)]/100;
  [p, res, cvg, outp] = run_opt(points, x0, start, h);
  if(res < res_best)
    p_best = p;
    res_best = res;
    noncont_best = i;
  endif
  disp(strcat("Result for i=", num2str(i), " is: ", num2str(res)));
endfor

save "sol_best.mat" p_best;
save "sol_150.mat" p;

%Dla siatki niejednorodnej
points = [(0: 1 : 24), (25 : 0.1 : 75), (76 : 1 : 200)];
start = [zeros(1,200), 40*ones(1,451)]/100;

%test
t0 = 0;
T = 1;
x0 = [0; 0];
h = 0.001;

points = t0 : 0.01 : T;
N = length(points)

cons = [(N-1)/2, (N - 3/2 : -1 : 1/2)]*step^2;
fin_pos = 0.9993465;
starting = ones(N,1)*2;
#starting(1) = 1;

pstarting = [(16/3 + 0.07)*ones((N-1)/4, 1); zeros((N+1)/2, 1); (-16/3 - 0.07)*ones((N-1)/4, 1)];

[p0, objf0, cvg, outp] = fmincon(OBJF = @(g) objf(g, points, x0, h),
				 X0   = pstarting,
				 A    = [],
				 B    = [],
				 AEQ  = cons',
				 BEQ  = fin_pos,
				 LB   = [],
				 UB   = [],
				 NONLCON = [],
				 OPTIONS=optimset("GradObj", "on")); %, "Algorithm",  "active-set"));

%TESTSSSTSTTST
load "res/res_grid_solutions"; %solutions

function s = start_bang(grid, point, val1, val2)
  pos = lookup(grid, point);
  s = [val1*ones(1, pos - 1), val2*ones(1, length(grid) - pos + 1)];
endfunction
start0 = @(grid) zeros(1, length(grid));
start_max = @(grid) gmax*ones(1,length(grid));
start40 = @(grid) start_bang(grid, 42.5, 0, 0.4);
start55 = @(grid) start_bang(grid, 42.5, 0, 0.55);
start3050 = @(grid) start_bang(grid, 10, 3, 0) + start_bang(grid, 50, 0, 0.5);

points = t0 : 1 : T;
points = t0 : 0.5 : T;
%points = t0 : 0.1 : T;
points = [(0: 1 : 24), (25 : 0.1 : 75), (76 : 1 : 200)];
start = start40(points);
start = start0(points);

discr = @const_discr;

p = solutions{1};
[y, dy] = objf(p, points, x0, 0.1, discr, constants);
[y0, dy0] = objf(start, points, x0, 0.1, discr, constants);

points = t0 : 0.5 : T;
start = [zeros(1,85), 55*ones(1,316)]/100;
[p, res, cvg, outp] = run_opt(points, x0, start, h, @const_discr, constants, 1e-9);
save "data/sol_big_prec.mat" p, res, cvg, outp;

[p, res, cvg, outp] = run_opt(points, x0, zeros(N,1), h, @const_discr, constants, 1e-9);

start = 2.9*ones(N,1);
[p_cc, res_cc, cvg_cc, outp_cc] = run_opt(points, x0, start, h, @const_discr, constants1, 1e-9);
[y_cc, dy_cc] = objf(p_cc, points, x0, h, @const_discr, constants1);
[y0_cc, dy0_cc] = objf(start, points, x0, h, @const_discr, constants1);


[p, res, cvg, outp] = run_opt(points, x0, zeros(N,1), h, @const_discr, constants, 1e-9);
%save "data/res_zero_prec.mat" p res cvg outp;
load "data/res_zero_prec.mat";
[y_zprc, dy_zprc] = objf(p, points, x0, h, @const_discr, constants);
[y0_zprc, dy0_zprc] = objf(zeros(N,1), points, x0, h, @const_discr, constants);

t0 = 0;
T = 200;
h = 0.05;
g0 = [0; 0.55; 42];
g0 = [0; 0.4; 50];
[p_bang, res_bang, cvg_bang, outp_bang] = run_opt_bang(t0, T, x0, g0, h, constants);
save "data/bang_res.mat" p_bang, res_bang, cvg_bang, outp_bang;

points = t0 : 0.5 : T;
N = length(points);
pt = lookup(points, p_bang(3));
start = [p_bang(1)*ones(1,pt), p_bang(2)*ones(1,N - pt)];
h = 0.1;
[p_ns, res_ns, cvg_ns, outp_ns] = run_opt(points, x0, start, h, @const_discr, constants, 1e-9);

lag = zeros(length(dy), 1);
for i = 1:length(dy)
  if p(i) < 1e-9
    disp(dy(i))
  endif
  if p(i) < 1e-8 && dy(i) < 0
    lag(i) = 0;
  elseif p(i) > 3 - (1e-9) && dy(i) > 0
    lag(i) = 0;
  else
    lag(i) = dy(i);
  endif
endfor

%Loading backups
function s = start_bang(grid, point, val1, val2)
  pos = lookup(grid, point);
  s = [val1*ones(1, pos - 1), val2*ones(1, length(grid) - pos + 1)];
endfunction
load "data/bang_res.mat"; %p_bang res_bang cvg_bang outp_bang
start_comp = @(grid) start_bang(grid, p_bang(3), p_bang(1), p_bang(2));
discr = @const_discr;
grid = t0 : 0.5 : T;
h = 0.1;
start = start_comp(grid);
[y0, dy0] = objf(start, grid, x0, h, discr, constants);
results = zeros(8, 5);
for i = 1:8
  load(["bak/res/res_prec_sol_" num2str(i)]); %p outp
  [y, dy] = objf(p, grid, x0, h, discr, constants);
  results(i,:) = [y, outp.niter, outp.nobjf, norm(dy,1), norm(dy0,1)];
endfor
save("-ascii", "res/res_prec", "results");

[p_scaled, res_scaled, cvg_scaled, outp_scaled] = run_opt(points, x0, start, h, @const_discr, constants, 1e-9);

				%SOME IMPORTANT RESULTS

save "data/p_50_6.mat" p_55_6; %start g_55, tol=10^-6
save "data/p_50_9.mat" p_55_9; %start g_55, tol=10^-9
save "data/p_comp_6.mat" p_comp_6; %start g_computed, tol=10^-6
save "data/p_comp_9.mat" p_comp_9; %start g_computed, tol=10^-9

				%Some old, probably obsolete stuff
load "sol0.mat";
load "sol_omega2k_1.mat";
load "sol_omega2k.mat";
save "sol0.mat" p0;
