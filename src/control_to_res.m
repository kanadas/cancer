
function res = control_to_res(g, points, x0, h, discr_fun, constants)

  function [x, dx] = approx_dyn(t, g, points, x0, h, discr_fun, constants)
    t0 = points(1);
    ghat = @(t) discr_fun(t, g, points);
    if nargout == 1
      x = rk4(@(t, V) dynamics(t, V, ghat, constants), t0, t, x0, h);
    else
      [x, dx] = rk4(@(t, V, dV) dynamics(t, V, ghat, constants, dV), t0, t, x0, h, length(g));
    endif
  endfunction
  
  res = @(t) approx_dyn(t, g, points, x0, h, discr_fun, constants);
endfunction

