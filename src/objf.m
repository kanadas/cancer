
function y = objf(g, points, x0)
  t0 = points(1);
  T = points(end);
%  N = length(points);

  approx_dyn = control_to_res(g, points, x0);
  
  function y = integral(t)
    epsilon = 0.01;
    omega = 1000; %???TODO
    function y = G(x)
      y = (1 + tanh(x))/2;
    endfunction

    V = approx_dyn(t);
    y = V(:, 1) + V(:, 2) + omega * G((V(:, 2) - V(:, 1)) / epsilon);
  endfunction
  
  y = quadgk(@(x) integral(x), t0, T);
endfunction
