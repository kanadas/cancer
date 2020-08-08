
function [y, dy] = integral_func(t, dyn)
  epsilon = 0.01;
%  omega = 1000;
  omega = 2000;
  function y = G(x)
    y = (1 + tanh(x))/2;
  endfunction

  function y = dG(x)
    y = sech(x) .^ 2 / 2;
  endfunction

  if nargout == 1
    V = dyn(t);
  else
    [V, dV] = dyn(t);
  endif
  y = V(:, 1) + V(:, 2) + omega * G((V(:, 2) - V(:, 1)) / epsilon);
  if nargout > 1
    dy = dV(:, 1, :) + dV(:, 2, :) + omega * dG((V(:, 2) - V(:, 1)) / epsilon) .* ((dV(:, 1, :) + dV(:, 2, :)) /epsilon);
    dy = squeeze(dy); %remove singleton dimension
  endif
endfunction

%TEST FUNCTION
%function [y, dy] = integral_func(t, dyn)
%  if nargout == 1
%    V = dyn(t);
%  else
%    [V, dV] = dyn(t);
%  endif
%
%  y = V(:,2) .^ 2;
%  if nargout > 1
%    dy = 2*V(:, 2) .* dV(:, 2, :);
%  endif
%endfunction

