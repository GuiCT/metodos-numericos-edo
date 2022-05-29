function u = AdamsBashforthStep(f, f0, u0, h, s)
% u = AdamsBashforthStep(f, t0, u0, h, s)
% Executa um passo do método numérico de Adams-Bashforth.
% INPUTS:
%   f = função du/dt = f(t, u)
%   f0 = valores anteriores de f(t, u)
%   u0 = último valor anterior de u
%   h = passo da malha
%   s = quantidade de estágios utilizados
% OUTPUTS:
%   u = próximo valor de u

  % Tabela de coeficientes de Adams-Bashforth.
  % Variável persistente que é declarada uma vez.
  persistent C = [
    [1, 0, 0, 0, 0, 0, 0, 0];
    [3/2, -1/2, 0, 0, 0, 0, 0, 0];
    [23/12, -4/3, 5/12, 0, 0, 0, 0, 0];
    [55/24, -59/24, 37/24, -3/8, 0, 0, 0, 0];
    [1901/720, -1387/360, 109/30, -637/360, 251/720, 0, 0, 0];
    [4277/1440, -2641/480, 4991/720, -3649/720, 959/480, -95/288, 0, 0];
    [198721/60480, -18637/2520, 235183/20160, -10754/945, 135713/20160, -5603/2520, 19087/60480, 0];
    [16083/4480, -1152169/120960, 242653/13440, -296053/13440, 2102243/120960, -115747/13440, 32863/13440, -5257/17280]
  ];
  
  % Último valor conhecido
  u = u0;
  % Adiciona termo a termo do método de Adams-Bashforth.
  u += h*dot(C(s,1:s), fliplr(f0));
endfunction