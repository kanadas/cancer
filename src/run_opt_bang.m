
function [p, res, cvg, outp] = run_opt_bang(t0, T, x0, g0, h, constants, TolFun=1e-9, algorithm = "lm_feasible")
  gmax = 3;  
  [p, res, cvg, outp] = fmincon(OBJF = @(x) objf_bang(x, x0, t0, T, h, constants),
				X0   = g0,
				A    = [],
				B    = [],
				AEQ  = [],
				BEQ  = [],
				LB   = [0; 0; t0],
				UB   = [gmax; gmax; T],
				NONLCON = [],
				OPTIONS=optimset("Algorithm", algorithm,
					         "MaxIter", 0,
					         "TolFun", TolFun,
						 "FinDiffRelStep", (1e-7)*ones(3,1)));
endfunction
