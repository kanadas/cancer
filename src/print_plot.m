
function print_plot(fname, p, points, x0, t0, T, discr_fun, constants)
  v = control_to_res(p, points, x0, (T - t0) / 2000, discr_fun, constants);
  control = @(t) discr_fun(t, p, points)
  subplot(2, 1, 1);
  x = t0:0.5:T;
  plot(x, control(x));
  axis([0, 200, 0, 3.5]);
  leg = legend('control');
  legend boxoff;
  set(leg, "fontsize", 12);
  subplot(2, 1, 2);
  plot(x, v(x));
  leg = legend({'V1', 'V2', 'K'});
  legend boxoff;
  set(leg, "fontsize", 12);

  print("-djpg", "-color", "-S800,800", fname);
endfunction
