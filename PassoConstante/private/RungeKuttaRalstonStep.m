function [u, f0] = RungeKuttaRalstonStep(f, t0, u0, h)
% u = RungeKuttaRalstonStep(f, t0, u0, h)
% Executa um passo do método numérico RK de Ralston.
% INPUTS:
%   f = função du/dt = f(t, u)
%   t0 = valor atual de t
%   u0 = valor atual de u
%   h = passo da malha
% OUTPUTS:
%   u = próximo valor de u
%   f0 = valor de f(t0, u0)

  persistent A = [0,          0,            0,          0;
                  0.4,        0,            0,          0;
                  0.29697761, 0.15875964,   0,          0;
                  0.21810040, -3.05096516,  3.83286476, 0];
  persistent B = [0.17476028, -0.55148066, 1.20553560, 0.17118478];
  persistent C = [0, 0.4, 0.45573725, 1];
  k = zeros(1, 4);
  for i=1:4
    k(i) = f(t0 + h*C(i), u0 + h*dot(A(i,:), k));
  endfor
  u = u0 + h*dot(B, k);
  f0 = k(1);
endfunction