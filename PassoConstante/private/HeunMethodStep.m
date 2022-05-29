function u = HeunMethodStep(f, t0, u0, h)
% u = HeunMethodStep(f, t0, u0, h)
% Executa um passo do método numérico de Heun.
% INPUTS:
%   f = função du/dt = f(t, u)
%   t0 = valor atual de t
%   u0 = valor atual de u
%   h = passo da malha
% OUTPUTS:
%   u = próximo valor de u

  k1 = f(t0, u0);
  predicted = u0 + h*k1;
  u = u0 + h*((k1 + f(t0 + h, predicted))/2);
endfunction