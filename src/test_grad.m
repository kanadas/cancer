
function test_grad(param, fun)
  [res, grad] = fun(param);  
  grad = squeeze(grad);
  hook.diffp = .00001 * ones (size (param))
  test_grad = squeeze(dfpdp(param, fun, hook));
  disp("computed gradient")
  if length(grad) <= 100
    disp(grad')
  else disp(grad(1:50)')
  endif
  disp("finite differences")
  if length(test_grad) <= 100
    disp(test_grad')
  else
    disp(test_grad(1:50)')
  endif
  disp("result")
  disp(norm(grad - test_grad, 1))

%  row_errors = sum(abs(grad - test_grad), 2);
%  [w, iw] = max(row_errors);
%  disp("max difference")
%  disp(w)
%  disp("on position")
%  disp(iw)
%  disp("grad values:")
%  disp(grad(iw,:)')
%  disp("fin dif values:")
%  disp(test_grad(iw,:)')
  
%  disp("res")
%  disp(res)
endfunction
