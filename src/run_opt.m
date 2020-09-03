function [p, res, cvg, outp] = run_opt(points, x0, g0, h, discr_fun, constants, algorithm = "lm_feasible", do_grad = "on")
  gmax = 3;
  N = length(points);
  maxIter = 0;
  if (algorithm == "active-set")
    maxIter = 200;
  endif
  [p, res, cvg, outp] = fmincon(OBJF = @(g) objf(g, points, x0, h, discr_fun, constants),
			       X0   = g0,
			       A    = [],
			       B    = [],
			       AEQ  = [],
			       BEQ  = [],
			       LB   = zeros(N, 1),
			       UB   = gmax * ones(N, 1),
			       NONLCON = [],
			       OPTIONS=optimset("GradObj", do_grad,
						"Algorithm", algorithm,
					        "MaxIter", maxIter,
					        "TolFun", 1e-6));
endfunction
