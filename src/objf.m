
function [y, dy] = objf(g, points, x0)
  STEP = 0.1;
  t0 = points(1);
  T = points(end);
  approx_dyn = control_to_res(g, points, x0, STEP);
  
  function [y, dy] = integral(t)
    epsilon = 0.01;
    omega = 1000; %???TODO
    function y = G(x)
      y = (1 + tanh(x))/2;
    endfunction

    function y = dG(x)
      y = sech(x) .^ 2 / 2;
    endfunction

    if nargout == 1
      V = approx_dyn(t);
    else
      [V, dV] = approx_dyn(t);
    endif
    y = V(:, 1) + V(:, 2) + omega * G((V(:, 2) - V(:, 1)) / epsilon);
    if nargout > 1
      dy = dV(:, 1, :) + dV(:, 2, :) + omega * dG((V(:, 2) - V(:, 1)) / epsilon) .* ((dV(:, 1, :) + dV(:, 2, :)) /epsilon);
      dy = squeeze(dy); %remove singleton dimension
    endif
  endfunction

  if nargout == 1
    y = integrate(@(x) integral(x), t0, T, STEP);
  else
    [y, dy] = integrate(@(x) integral(x), t0, T, STEP);
  endif
  
  disp(y);
%  if nargout > 1
%    disp(dy);
%  endif

  if nargout > 1
    %Seems that it wants - gradient
    dy = -dy;
  endif
endfunction
