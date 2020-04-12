
%%% Constants

global lambda1 = 0.192;
global lambda2 = 0.192;
global mu = 0.0;
global b1 = 5.85;
global b2 = 5.85;
global d = 0.00873;
global beta1 = 0.15;
global beta2 = 0.1;
global beta = 0.05;
global gmax = 3;
global alpha12 = 0.1; %???TODO
global alpha21 = 0.15; %???TODO
global t0 = 0;
global T = 200;
global x0 = [20, 280, 650];
function y = F(x)
  y = -log(x);
endfunction

global epsilon = 0.01;
global omega = 1000; %???TODO
function y = G(x)
  y = (1 + tanh(x))/2;
endfunction

%%% Auxiliary definitions

function x = ode(f, t0, T, x0)
  %Constant length of step for now
  STEP = 0.1; 
  N = ceil((T - t0)/STEP);
  x = rk4(f, t0, x0, T, N);
endfunction

%%% Differential equation

function v = f(t, V, g)
  global lambda1;
  global lambda2;
  global alpha12;
  global alpha21;
  global beta1;
  global beta2;
  global beta;
  global mu;
  global d;

%  disp(V);
  
  [V1, V2, K] = num2cell(V){:};
  v1 = lambda1 * V1 .* F((V1 + alpha12 .* V2) ./ K) - beta1 * V1 .* g(t);
  v2 = lambda2 * V2 .* F((V2 + alpha21 .* V1) ./ K) - beta2 * V2 .* g(t);
  k = -mu * K + (V1 + V2) - d * ((V1 + V2).^(2/3)) .* K - beta * K .* g(t);
  v = [v1 v2 k];
endfunction

%%% Goal function

%% collocation points
step = 0.2;
global tt = t0 : step : T;
global N = (T - t0) / step;

function y = J(g)
  global t0;
  global T;
  global N;
  global tt;
  global x0;

%  disp(g);
  
  X = [x0];
  for i = 1 : N - 1
%    disp(X(end, :))   
    X = [X ; ode(@(t, V) f(t, V, @(x) g(i)), tt(i), tt(i + 1), X(end, :))];
  end
  
  gg = @(t) g(lookup(tt, t)); 
  
  function x = approxf(t)
    pos = lookup(tt, t);
    ster = g(pos);
    x = ode(@(t, V) f(t, V, @(x) ster), tt(pos), tt(pos + 1), X(pos, :))(end, :, :);
  endfunction

  function y = integral(t)
    global omega;
    global epsilon;
    V = approxf(t);
    y = V(1) + V(2) + omega * G((V(2) - V(1)) / epsilon);

    disp(V);
    disp(omega);
    disp(epsilon);
    disp(y);
  endfunction
  
  y = quadgk(@(x) integral(x), 0, T)
endfunction

[p, objf, cvg, outp] = nonlin_min(@J, zeros(N, 1), optimset("lbound", zeros(N, 1), "ubound", gmax * ones(N, 1)));

