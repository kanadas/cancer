#{
function res = control_to_res(g, points, x0)

  function x = approx_dyn(t, g, points, x0)
    t = t';
    t0 = points(1);
    ghat = @(t) ster(t, g, points)
    if st(1) == t0
      is_t0 = true;
    else
      is_t0 = false;
      st = [t0, st];
      id = [0, id];
    endif
    [rt, y] = ode45(@(t, V) dynamics(t, V, ster), st, x0);
    y = [id' y];
    y = sortrows(y, [1]);
    y = y(:, 2:end);
    if !is_t0
      y = y(2:end, :);
    endif
    x = y;
  endfunction
  
  res = @(t) approx_dyn(t, g, points, x0);  
endfunction
#}

function res = control_to_res(g, points, x0, h)

  function [x, dx] = approx_dyn(t, g, points, x0, h)
    t0 = points(1);
    ghat = @(t) ster(t, g, points);
    if nargout == 1
      x = rk4(@(t, V) dynamics(t, V, ghat), t0, t, x0, h);
    else
      [x, dx] = rk4(@(t, V, dV) dynamics(t, V, ghat, dV), t0, t, x0, h, length(g));
    endif
  endfunction
  
  res = @(t) approx_dyn(t, g, points, x0, h);
endfunction

