
function run_prec_exp()
  pkg load optim;
  gmax = 3;
  t0 = 0;
  T = 200;
  x0 = [20; 280; 650];

  c1.lambda1 = 0.192;
  c1.lambda2 = 0.192;
  c1.mu = 0.0;
  c1.b1 = 5.85;
  c1.b2 = 5.85;
  c1.d = 0.00873;
  c1.beta1 = 0.15;
  c1.beta2 = 0.1;
  c1.beta = 0.05;
  c1.alpha12 = 0.1;
  c1.alpha21 = 0.15;
  c1.epsilon = 0.01;
  c1.omega = 1000;

  c2 = c1;
  c2.alpha12 = 0.5;
  c2.alpha21 = 0.75;
  c2.omega = 2000;

  function s = start_bang(grid, point, val1, val2)
    pos = lookup(grid, point);
    s = [val1*ones(1, pos - 1), val2*ones(1, length(grid) - pos + 1)];
  endfunction
  load "data/bang_res.mat"; %p_bang res_bang cvg_bang outp_bang
  start_comp = @(grid) start_bang(grid, p_bang(3), p_bang(1), p_bang(2));

  constants = c2;
  backends = {"lm_feasible", "active-set"};
  discr = @const_discr;
  grid = t0 : 0.5 : T;
  h = 0.1;
  start = start_comp(grid);
  precisions = {
		1e-6,
		1e-9,
		1e-11,
		1e-13};

  results = zeros(length(backends)*length(precisions), 5);
  solutins = cell(length(backends)*length(precisions));
  idx = 1;
  for i = 1:length(backends)
    [y0, dy0] = objf(start, grid, x0, h, discr, constants);
    for j = 1:length(precisions)
      disp(strcat("Processing experiment  ", int2str(idx) ));
      [p, res, cvg, outp] = run_opt(grid,
				    x0,
				    start,
				    h,
				    discr,
				    constants,
				    precisions{j},
				    backends{i});
      [y, dy] = objf(p, grid, x0, h, discr, constants);
      save(["bak/res/res_prec_sol_" num2str(idx)], "p", "outp");
      solutions{idx} = p;
      results(idx,:) = [res, outp.niter, outp.nobjf, norm(dy,1), norm(dy0,1)];
      idx = idx + 1;
    endfor
  endfor
  save("res/res_prec_solutions", "solutions");
  save("-ascii", "res/res_prec");

endfunction
