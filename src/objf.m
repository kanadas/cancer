
function [y, dy] = objf(g, points, x0, h, discr_fun, constants)
  t0 = points(1);
  T = points(end);
  approx_dyn = control_to_res(g, points, x0, h, discr_fun, constants);

  if nargout == 1
    y = integrate(@(x) integral_func(x, approx_dyn, constants), t0, T, h);
  else
    [y, dy] = integrate(@(x) integral_func(x, approx_dyn, constants), t0, T, h);
  endif

  if nargout == 1
    disp(y);
  else
    disp(strcat(num2str(y), ", ", num2str(norm(dy, 1))))
  endif
endfunction
