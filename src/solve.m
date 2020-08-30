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

load "sol0.mat";
load "sol_omega2k_1.mat";
load "sol_omega2k.mat";
save "sol0.mat" p0;

points = t0 : 4 : T;

points = t0 : 1 : T;

points = t0 : 0.5 : T;

points = [(0 : 0.5 : 50), (51 : 1 : 149), (150 : 0.5 : 200)];

points = [(0: 1 : 24), (25 : 0.1 : 75), (76 : 1 : 200)];

N = length(points);

start_bang = 1:N;
start_bang = start_bang < N/2;

[p, res, cvg, outp] = run_opt(points, x0, zeros(N,1), h, @const_discr, constants);
#no gradient calculation experimment
[p_ng, res_ng, cvg_ng, outp_ng] = run_opt(points, x0, zeros(N,1), h, @const_discr, constants, "lm_feasible", "off");

save "sol_nograd.mat" p_ng res_ng outp_ng;
load "data/sol_nograd.mat";

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
start = [zeros(1,90), 55*ones(1,311)]/100;
start = [zeros(1,85), 55*ones(1,316)]/100;
start = [300*ones(1,20), zeros(1, 80), 50*ones(1,301)]/100;

[p, res, cvg, outp] = run_opt(points, x0, start, h, @const_discr, constants);
[p, res, cvg, outp] = run_opt(points, x0, start, h, @const_discr, constants, "active-set");
[p, res, cvg, outp] = run_opt(points, x0, start_lin, h, "active-set");

save "sol_bang.mat" p;

[p_c1, res_c1, cvg_c1, outp_c1] = run_opt(points, x0, zeros(N,1), h, @const_discr, constants1, "active-set");
save "data/sol_test_c1.mat" p_c1, res_c1, outp_c1;

%Test sqp siatka S_1 start40 (na testach wyszedł absurdalny wynik)
points = t0 : 1 : T;
N = length(points);
start = [zeros(1,42), 40*ones(1,159)]/100;
[p_sqs1, res_sqs1, cvg_sqs1, outp_sqs1] = run_opt(points, x0, start, h, @const_discr, constants, "active-set");
save "data/sol_test_sqs1.mat" p_sqs1, res_sqs1, outp_sqs1;

%BEST
% 2. & {\it sqp\/} & stała & $S_{0.5}$ & 0.1 & $g_{0,0.55}$ & 2.71 & 10 & 130 \\
points = t0 : 0.5 : T;
N = length(points);
start = [zeros(1,85), 55*ones(1,316)]/100;
[p_min, res_min, cvg_min, outp_min] = run_opt(points, x0, start, h, @const_discr, constants, "active-set");

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



%{
Statystyki:
nowe parametry:
	Wynik dla:
	zero control: 381831.74497
	85x0 + 316x40:	301226.27469
	85x0 + 316x55:  277858.64599
	Eksperymenty:
	     start 0:
	     niter =  6
	     nobjf =  17
	     res =  324132.85096
	     start bang:
	     niter =  10
	     nobjf =  26
	     res =  283519.85873
	     start best:
	     niter =  19
	     nobjf =  42
	     res =  272466.20561
	     active_set zero start:
	     niter =  3
	     nobjf =  57
	     res =  340595.93819
	     active set best start:
	     niter =  7
             nobjf =  110
	     res =  272502.70105
	Siatka niejednorodna:
	       start 0:
	       niter =  20
	       nobjf =  42
	       res =  295303.82581
	       active set zero start:
	       niter =  2
	       nobjf =  53
	       res =  348091.22831
	       start bang:
	       niter =  14
	       nobjf =  35
	       res =  277775.51745
	       start best:
    	       niter =  19
	       nobjf =  42
	       res =  272848.18541
stare parametry:
      Wynik dla	 
      zero control: 567688.31755
      max control:  213575.68897
      Eksperymenty:
             start 0:
	     niter =  6
	     nobjf =  15
	     res =  234437.61671
	     max start:
	     niter =  1
	     nobjf =  2
    	     res =  213575.68897
	     active set zero start:
	     niter =  8
	     nobjf =  9
	     res =  213575.68897

