function u = RungeKutta4Step(f, t0, u0, h)
% u = RK4Step(f, t0, u0, h)
% Executa um passo do método numérico RK4.
% INPUTS:
%   f = função du/dt = f(t, u)
%   t0 = valor atual de t
%   u0 = valor atual de u
%   h = passo da malha
% OUTPUTS:
%   u = próximo valor de u

  persistent A = [0,    0,    0,  0;
                  1/2,  0,    0,  0;
                  0,    1/2,  0,  0;
                  0,    0,    1,  0];
  persistent B = [1/6, 1/3, 1/3, 1/6];
  persistent C = [0, 1/2, 1/2, 1];
  k = zeros(1, 4);
  for i=1:4
    k(i) = f(t0 + h*C(i), u0 + h*dot(A(i,:), k));
  endfor
  u = u0 + h*dot(B, k);
endfunction