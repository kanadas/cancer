%Run all experiments
function run_experiments()
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

  constants = {c1, c2};
  backends = {"lm_feasible", "active-set"};
  discretizations = {@const_discr, @linear_discr};
  grids = {t0 : 1 : T,
	   t0 : 0.5 : T,
	   [(0 : 0.5 : 50), (51 : 1 : 149), (150 : 0.5 : 200)],
	   [(0: 1 : 24), (25 : 0.1 : 75), (76 : 1 : 200)]
	  };
  steps = {0.5, 0.1, 0.02};

  function s = start_bang(grid, val)
    pos = lookup(grid, 42.5);
    s = [zeros(1, pos - 1), val*ones(1, length(grid) - pos + 1)];
  endfunction
  start0 = @(grid) zeros(1, length(grid));
  start_max = @(grid) gmax*ones(1,length(grid));
  start40 = @(grid) start_bang(grid, 0.4);
  start55 = @(grid) start_bang(grid, 0.55);
  starts = {start0, start_max, start40, start55};

  results = [];
  for const = constants
    for back = backends
      for discr = discretizations
	for grid = grids
	  for h = steps
	    for startf = starts
	      disp(strcat("Processing experiment ", num2str(length(results))));
	      start = startf{1}(grid{1});
	      [p, res, cvg, outp] = run_opt(grid{1}, x0, start, h{1}, discr{1}, const{1}, back{1});
	      results(end + 1) = res;
	    endfor
	  endfor
	endfor
      endfor
    endfor
  endfor

  save -ascii "results" results;
endfunction

