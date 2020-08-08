% Numerical integration with constant step length function
% h has to be divisor of T - t0
function [x, dx] = integrate(f, t0, T, h)
  N = (T - t0) / h;
  %Trapezoidal quadrature for now
  trap_quad = ones(N + 1, 1) * h;
  trap_quad(1) /= 2;
  trap_quad(end) /= 2;
  
  eval_points = t0:h:T;
  if nargout == 1
    y = f(eval_points);
  else
    [y, dy] = f(eval_points);
  endif

  x = dot(trap_quad, y, 1);
  if nargout > 1
    dx = dy' * trap_quad;
  endif
endfunction
