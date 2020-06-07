%% Rungge-Kutta 4th order
% x' = f(t, x), x(t0) = x0
% pts: points of evaluation, sorted
% h: length of step, has to be a divisor of length of every interval
% TODO: gradients may be sparse, check if can use that
function [res, grad] = rk4(f, t0, pts, x0, h, grad_dim = 0)
  t = t0;
  x = x0;
  dx = zeros(length(x0), grad_dim);
  res = zeros(length(pts), length(x0));
  grad = zeros(length(pts), length(x0), grad_dim);
  i = 1;
  for p = pts
    while (t < p)
      if nargout == 1
	k1 = f(t, x);
	k2 = f(t + h/2, x + h/2*k1);
	k3 = f(t + h/2, x + h/2*k2);
	k4 = f(t + h, x + h*k3);
      else
	[k1, dfdx, dfdg] = f(t, x, dx);
	dk1 = dfdx * dx + dfdg;
	[k2, dfdx, dfdg] = f(t + h/2, x + h/2*k1, dx + h/2*dk1);
	dk2 = dfdx * (dx + h/2*dk1) + dfdg;
	[k3, dfdx, dfdg] = f(t + h/2, x + h/2*k2, dx + h/2*dk2);
	dk3 = dfdx * (dx + h/2*dk2) + dfdg;
	[k4, dfdx, dfdg] = f(t + h, x + h*k3, dx + h*dk3);
	dk4 = dfdx * (dx + h*dk3) + dfdg;

	dx = dx + h/6*(dk1 + 2*dk2 + 2*dk3 + dk4);
      endif
      x = x + h/6*(k1 + 2*k2 + 2*k3 + k4);
      t = t + h;
    endwhile
    res(i,:) = x;
    if nargout > 1
      grad(i, :, :) = dx;
    endif
    i += 1;
  end
end
