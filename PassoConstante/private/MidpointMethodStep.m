function u = MidpointMethodStep(f, t0, u0, h)
% u = MidpointMethodStep(f, t0, u0, h)
% Executa um passo do método numérico do ponto médio.
% INPUTS:
%   f = função du/dt = f(t, u)
%   t0 = valor atual de t
%   u0 = valor atual de u
%   h = passo da malha
% OUTPUTS:
%   u = próximo valor de u

  u = u0 + h*f(t0 + h/2, u0 + (h/2)*f(t0, u0));
endfunction