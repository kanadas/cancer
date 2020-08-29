function [v, dv, dg] = dynamics(t, V, g, constants)  
  lambda1 = constants.lambda1;
  lambda2 = constants.lambda2;
  mu = constants.mu;
  b1 = constants.b1;
  b2 = constants.b2;
  d = constants.d;
  beta1 = constants.beta1;
  beta2 = constants.beta2;
  beta = constants.beta;
  alpha12 = constants.alpha12;
  alpha21 = constants.alpha21;
  function y = F(x)
    y = -log(x);
  endfunction

  function y = dF(x)
    y = -1 ./ x;
  endfunction

  if nargout == 1
    ggt = g(t);
  else
    [ggt, dgt] = g(t);
  endif
  
  V1 = V(1);
  V2 = V(2);
  K = V(3);
  frac1 = (V1 + alpha12 * V2) ./ K;
  frac2 = (V2 + alpha21 * V1) ./ K;
  pow = (V1 + V2).^(2/3);
  v1 = lambda1 * V1 .* F(frac1) - beta1 * V1 .* ggt;
  v2 = lambda2 * V2 .* F(frac2) - beta2 * V2 .* ggt;
  k = -mu * K + (V1 + V2) - d * pow .* K - beta * K .* ggt;
  v = [v1; v2; k];

  if nargout > 1 %gradient required
    dpow = 2/3 * pow ./ (V1 + V2);
    dv1dv1 = lambda1 * (F(frac1) + V1 .* dF(frac1) ./ K) - beta1 .* ggt;
    dv1dv2 = lambda1 * V1 .* dF(frac1) .* (alpha12 / K);
    dv1dk = lambda1 * V1 .* dF(frac1) .* ( -1 * frac1 ./ K);
    dv2dv1 = lambda2 * V2 .* dF(frac2) .* (alpha21 / K);
    dv2dv2 = lambda2 * (F(frac2) + V2 .* dF(frac2) ./ K) - beta2 .* ggt;
    dv2dk = lambda2 * V2 .* dF(frac2) .* ( -1 * frac2 ./ K);
    dkdv1 = 1 - d * dpow .* K;
    dkdv2 = 1 - d * dpow .* K;
    dkdk = -mu - d * pow - beta * ggt;
    dv = [dv1dv1, dv1dv2, dv1dk; dv2dv1, dv2dv2, dv2dk; dkdv1, dkdv2, dkdk];

    dv1dg = -beta1 .* V1 .* dgt;
    dv2dg = -beta2 .* V2 .* dgt;
    dkdg = -beta .* K .* dgt;
    dg = [dv1dg; dv2dg; dkdg];
  endif
endfunction

%TEST FUNCTION
%function [v, dv, dg] = dynamics(t, V, g, dV = [])
%  if nargout == 1
%    u = g(t);
%  else
%    [u, du] = g(t);
%  endif
%  px = V(2);
%  pv = u;
%  v = [px; pv];
%
%  if nargout > 1
%    dxdx = 0;
%    dxdv = 1;
%    dvdx = 0;
%    dvdv = 0;
%    dv = [dxdx, dxdv; dvdx, dvdv];
%
%    dxdg = dV(2,:);
%    dvdg = du;
%    dg = [dxdg; dvdg];
%  endif
%endfunction

