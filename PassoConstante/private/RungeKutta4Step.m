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

  k1 = f(t0, u0);
  k2 = f(t0 + h/2, u0 + (h/2)*k1);
  k3 = f(t0 + h/2, u0 + (h/2)*k2);
  k4 = f(t0 + h, u0 + h*k3);
  u = u0 + h*(k1/6 + k2/3 + k3/3 + k4/6);
endfunction