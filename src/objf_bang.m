%x[1] - first value
%x[2] - second value
%x[3] - discontinouity point
function y = objf_bang(x, x0, t0, T, h, constants, discr_fun = @const_discr)
  points = [t0; x(3); T];
  g = [x(1); x(2); x(2)];
  y = objf(g, points, x0, h, discr_fun, constants);
endfunction
