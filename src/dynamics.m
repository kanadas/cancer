function v = dynamics(t, V, g)
  lambda1 = 0.192;
  lambda2 = 0.192;
  mu = 0.0;
  b1 = 5.85;
  b2 = 5.85;
  d = 0.00873;
  beta1 = 0.15;
  beta2 = 0.1;
  beta = 0.05;
  alpha12 = 0.1; %???TODO
  alpha21 = 0.15; %???TODO
  function y = F(x)
    y = -log(x);
  endfunction

  V1 = V(1);
  V2 = V(2);
  K = V(3);
  v1 = lambda1 * V1 .* F((V1 + alpha12 .* V2) ./ K) - beta1 * V1 .* g(t);
  v2 = lambda2 * V2 .* F((V2 + alpha21 .* V1) ./ K) - beta2 * V2 .* g(t);
  k = -mu * K + (V1 + V2) - d * ((V1 + V2).^(2/3)) .* K - beta * K .* g(t);
  v = [v1; v2; k];
endfunction
