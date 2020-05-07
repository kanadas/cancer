function res = control_to_res(g, points, x0)
  N = length(points) - 1;
  X = [x0];
  t0 = points(1);
  T = points(end);

  %Preprocess values at collocation points to return faster function
  for i = 1 : N - 1
    X = [X ; ode(@(t, V) dynamics(t, V, @(x) g(i)), points(i), points(i + 1), X(end, :))];
  end
  
  function x = approx_dyn(t, points, g, X)
    pos = lookup(points, t);
    ster = g(pos);
    zero_idx = points(pos)' == t;
    nz = points(pos)' != t;
    x(zero_idx, :) = X(pos(zero_idx), :);
    x(nz, :) = ode(@(t, V) dynamics(t, V, @(x) ster(nz)), points(pos(nz)), t(nz), X(pos(nz), :));
  endfunction

  res = @(t) approx_dyn(t, points, g, X);
endfunction
