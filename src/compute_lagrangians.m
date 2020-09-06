
function compute_lagrangians()
  pkg load optim;
  gmax = 3;
  t0 = 0;
  T = 200;
  x0 = [20; 280; 650];
  function run(constants, backends, discretizetions, grids, steps, starts, precisions, fname)

    function L = lagrangian(dy, p)
      L = zeros(length(dy), 1);
      for zz = 1:length(dy)
	if p(zz) < 1e-8 && dy(zz) > 0
	  L(zz) = 0;
	elseif p(zz) > 3 - (1e-9) && dy(zz) < 0
	  L(zz) = 0;
	else
	  L(zz) = dy(zz);
	endif
      endfor
    endfunction
    
    load([fname "_solutions"]);
    results = zeros(length(solutions), 2);
    idx = 1;
    for i = 1:length(constants)
      for j = 1:length(backends)
	for k = 1:length(discretizations)
	  for l = 1:length(grids)
	    for m = 1:length(steps)
              for n = 1:length(starts)
		for o = 1:length(precisions)
		  disp(strcat("Processing experiment  ", int2str(idx) ));
		  start = starts{n}(grids{l});
		  p = solutions{idx};
		  [y, dy] = objf(p, grids{l}, x0, steps{m}, discretizations{k}, constants{i});
		  L = lagrangian(dy, p);
		  [y0, dy0] = objf(start, grids{l}, x0, steps{m}, discretizations{k}, constants{i});
		  L0 = lagrangian(dy0, start);
		  results(idx,:) = [norm(L,1), norm(L0,1)];
		  idx = idx + 1;
		endfor
	      endfor
	    endfor
	  endfor
	endfor
      endfor
    endfor

    save("-ascii", [fname "_lagrangians"], "results");
  endfunction

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
  start0 = @(grid) zeros(1, length(grid));
  start_max = @(grid) gmax*ones(1,length(grid));
  start40 = @(grid) start_bang(grid, 42.5, 0, 0.4);
  start55 = @(grid) start_bang(grid, 42.5, 0, 0.55);
  load "data/bang_res.mat"; %p_bang res_bang cvg_bang outp_bang
  start_comp = @(grid) start_bang(grid, p_bang(3), p_bang(1), p_bang(2));
%  start3050 = @(grid) start_bang(grid, 10, 3, 0) + start_bang(grid, 50, 0, 0.5);

% All experiments  
%  constants = {c1, c2};
%  backends = {"lm_feasible", "active-set"};
%  discretizations = {@const_discr, @linear_discr};
%  grids = {t0 : 1 : T,
%	   t0 : 0.5 : T,
%  	   t0 : 0.1 : T,
%	   [(0 : 0.5 : 50), (51 : 1 : 149), (150 : 0.5 : 200)],
%	   [(0: 1 : 24), (25 : 0.1 : 75), (76 : 1 : 200)]
%	  };
%  steps = {0.5, 0.1, 0.02};
%  starts = {start0, start_max, start40, start55, start3050};
%  precisions = {1e-6, 1e-9, 1e-10, 1e-11, 1e-12}

% start test
%  disp("TEST START")
%  constants = {c2};
%  backends = {"lm_feasible", "active-set"};
%  discretizations = {@const_discr};
%  grids = {t0 : 0.5 : T};
%  steps = {0.1};
%  starts = {start0,
%	    start_max,
%	    start55,
%	    start_comp,
%	   };
%  precisions = {1e-6, 1e-9};
%  run(constants, backends, discretizations, grids, steps, starts, precisions, "res/res_start")
%  
%% CC test  
%  disp("TEST (CC)")
%  constants = {c1};
%  backends = {"lm_feasible"};
%  discretizations = {@const_discr, @linear_discr};
%  grids = {t0 : 1 : T,
%	   t0 : 0.5 : T,
%	   [(0 : 0.5 : 50), (51 : 1 : 149), (150 : 0.5 : 200)],
%	  };
%  steps = {0.1};
%  starts = {start0, start_max};
%  precisions = {1e-6, 1e-9};
%  run(constants, backends, discretizations, grids, steps, starts, precisions, "res/res_paramCC")
%
%% Discretization test
%  disp("TEST DISCRETIZATION");
%  constants = {c2};
%  backends = {"lm_feasible", "active-set"};
%  discretizations = {@linear_discr}; %const computed in start test
%  grids = {t0 : 0.5 : T};
%  steps = {0.1};
%  starts = {start_comp};
%  precisions = {1e-9};
%  run(constants, backends, discretizations, grids, steps, starts, precisions, "res/res_discr")
%
%% Grid test
%  disp("TEST GRID")
%  constants = {c2};
%  backends = {"lm_feasible", "active-set"};
%  discretizations = {@const_discr};
%  grids = {t0 : 1 : T,
%%	   t0 : 0.5 : T, computed in start test
%%	   t0 : 0.1 : T,
%	   [(0: 1 : 24), (25 : 0.1 : 75), (76 : 1 : 200)]
%	  };
%  steps = {0.1};
%  starts = {start_comp};
%  precisions = {1e-9};
%  run(constants, backends, discretizations, grids, steps, starts, precisions, "res/res_grid")
%
%% h test
%  disp("TEST h")
%  constants = {c2};
%  backends = {"lm_feasible", "active-set"};
%  discretizations = {@const_discr};
%  grids = {t0 : 0.5 : T};
%  steps = {0.5,
%%	   0.1, computed in tart test
%	   0.02};
%  starts = {start_comp};
%  precisions = {1e-9};
%  run(constants, backends, discretizations, grids, steps, starts, precisions, "res/res_h")
%
  constants = {c2};
  backends = {"lm_feasible", "active-set"};
  discretizations = {@const_discr};
  grids = {t0 : 0.5 : T};
  steps = {0.1};
  starts = {start_comp};
  precisions = {
		1e-6,
		1e-9,
		1e-11,
		1e-13};
  run(constants, backends, discretizations, grids, steps, starts, precisions, "res/res_prec")
  
endfunction
