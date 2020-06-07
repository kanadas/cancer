%Piecewise constant function
%TODO: linear and Simpson collocation
function [y, grad] = ster(t, g, points)
  pos = lookup(points, t);
  y = g(pos);

  if nargout > 1 %gradient required
    grad = zeros(length(t), length(g));
    grad(:,pos) = 1;
  endif
endfunction
