%Piecewise linear discretization
function [y, grad] = linear_discr(t, g, points)
  pos = lookup(points, t);
  
  d1 = t - points(pos);
  if d1 > 1e-10
    d2 = points(pos + 1) - t;
    y = (d2 .* g(pos) + d1 .* g(pos + 1)) ./ (d1 + d2);
  else
    y = g(pos);
  endif

  if nargout > 1 %gradient required
    grad = zeros(length(t), length(g));
    if d1 > 1e-10
      grad(:,pos) = d2 ./ (d1 + d2);
      grad(:,pos + 1) = d1 ./ (d1 + d2);
    else
      grad(:,pos) = 1;
    endif
  endif
endfunction

