
function [y, dy] = objf(g, points, x0, h)
  t0 = points(1);
  T = points(end);
  approx_dyn = control_to_res(g, points, x0, h);

  if nargout == 1
    y = integrate(@(x) integral_func(x, approx_dyn), t0, T, h);
  else
    [y, dy] = integrate(@(x) integral_func(x, approx_dyn), t0, T, h);

  endif
  
  disp(y);
%  if nargout > 1
%    disp(dy);
%  endif

  if nargout > 1
    %Seems that it wants - gradient
%    dy = -dy;
    dy = dy';
  endif
endfunction
