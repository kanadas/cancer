#{
function res = control_to_res(g, points, x0)
  N = length(points) - 1;
  t0 = points(1);
  T = points(end);

%  X = [x0];
  %Preprocess values at collocation points to return faster function
%  for i = 1 : N - 1
%    X = [X ; ode(@(t, V) dynamics(t, V, @(x) g(i)), points(i), points(i + 1), X(end, :))];
%  end

  function x = approx_dyn(t, points, g, X)

%    persistent calls = 0;
%    printf ("called %d times\n", ++calls);
%    if(issorted(t, "ascend")) disp("sorted") else disp("not sorted") endif
        
    pos = lookup(points, t);
    ster = g(pos);
    zero_idx = points(pos)' == t;
    nz = points(pos)' != t;    
    x(zero_idx, :) = X(pos(zero_idx), :);
    x(nz, :) = ode(@(t, V) dynamics(t, V, @(x) ster(nz)), points(pos(nz)), t(nz), X(pos(nz), :));
  endfunction


  res = @(t) approx_dyn(t, points, g);
endfunction
#}

function res = control_to_res(g, points, x0)

  function x = approx_dyn(t, g, points, x0)
    t = t';
    t0 = points(1);
    ster = @(t) g(lookup(points, t));
    [st, id] = sort(t);
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

