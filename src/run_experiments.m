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

%  factors = {backends, discretizations, grids, steps, starts};
%  disp(length(factors));
%  function exps = product(factors, acc)
%    len = length(factors);
%    if(len == 0)
%      exps = acc;
%    else
%      factor = factors{end};
%      nacc = {};
%      for exp = acc
%	for fc = factor
%	  exp{end+1} = fc;
%	  nacc{end+1} = exp;
%	endfor
%      endfor
%      disp(length(nacc))
%      exps = product(factors(1:(len - 1)), nacc);
%    endif
%  endfunction
  
  results = zeros(length(constants)*length(backends)*length(discretizations)*length(grids)*length(steps)*length(starts), 3);
  idx = 1;
  for i = 1:length(constants)
    for j = 1:length(backends)
      for k = 1:length(discretizations)
	for l = 1:length(grids)
	  for m = 1:length(steps)
            for n = 1:length(starts)
	      disp(strcat("Processing experiment  ", int2str(idx) ));
              start = starts{n}(grids{l});
              [p, res, cvg, outp] = run_opt(grids{l},
					    x0,
					    start,
					    steps{m},
					    discretizations{k},
					    constants{i},
					    backends{j});
              results(idx,:) = [res, outp.niter, outp.nobjf];
	      idx = idx + 1;
            endfor
	  endfor
	endfor
      endfor
    endfor
  endfor  

  disp(length(experiments));

 results = [];


save -ascii "results_all" results;
endfunction

