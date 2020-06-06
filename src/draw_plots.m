function draw_plots(g, points, x0)
  t0 = points(1);
  T = points(end);
  v = control_to_res(g, points, x0);
  step = (T - t0) / 100;
  x = t0:step:(T - step);
  function y = control(x)
    pos = lookup(points, x);
    y = g(pos);
  endfunction

  subplot(2, 1, 1);
  plot(x, control(x));
  legend('control');
  subplot(2, 1, 2);
  plot(x, v(x'));
  legend('V1', 'V2', 'K');
endfunction