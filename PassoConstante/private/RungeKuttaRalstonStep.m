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

  k1 = f(t0, u0);
  k2 = f(t0 + 0.4*h, u0 + 0.4*h*k1);
  k3 = f(t0 + 0.45573725*h, u0 + 0.29697761*h*k1 + 0.15875964*h*k2);
  k4 = f(t0 + h, u0 + 0.21810040*h*k1 - 3.05096516*h*k2 + 3.83286476*h*k3);
  u = u0 + h*(0.17476028*k1 - 0.55148066*k2 + 1.20553560*k3 + 0.17118478*k4);
  f0 = k1;
endfunction