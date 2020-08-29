
function [y, dy] = objf(g, points, x0, h, discr_fun, constants)
  t0 = points(1);
  T = points(end);
  approx_dyn = control_to_res(g, points, x0, h, discr_fun, constants);

  if nargout == 1
    y = integrate(@(x) integral_func(x, approx_dyn, constants), t0, T, h);
  else
    [y, dy] = integrate(@(x) integral_func(x, approx_dyn, constants), t0, T, h);
  endif
  
  disp(y);

  global DEBUG;
  if DEBUG > 1 && nargout > 1
    test_grad = dfpdp(g, @(g) objf(g, points, x0, h, discr_fun, constants));
    diff_norm = norm(dy - test_grad, 1);
    disp(strcat("OBJF grad diff norm: ", num2str(diff_norm)));
    if diff_norm >= 1
      disp("g: ")
      disp(g)
      disp("Computed gradient");
      disp(dy');
      disp("Finite diff gradient");
      disp(test_grad');
      save "bug_log" g dy test_grad;
      exit;
    endif
  endif
endfunction
