%Simple function to easily change algorithm to solve ode
function x = ode(f, t0, T, x0)
  %Constant length of step for now
  STEP = 0.1;
  N = ceil((T - t0)/STEP);
  x = rk4(f, t0, x0, T, N);
%  [t,v] = ode45(f, [t0, T], x0);
%  x = v(end, :);
endfunction
