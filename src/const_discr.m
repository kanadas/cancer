%Piecewise constant discretization
function [y, grad] = const_discr(t, g, points)
  pos = lookup(points, t);
  y = g(pos);

  if nargout > 1 %gradient required
    grad = zeros(length(t), length(g));
    grad(:,pos) = 1;
  endif
endfunction

