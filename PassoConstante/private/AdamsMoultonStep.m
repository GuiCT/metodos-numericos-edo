function u = AdamsMoultonStep(f, f0, u0, f_predicted, h, s)
% u = AdamsMoultonStep(f, t0, u0, h, s)
% Executa um passo do método numérico de Adams-Moulton.
% INPUTS:
%   f = função du/dt = f(t, u)
%   f0 = valores anteriores de f(t, u)
%   u0 = último valor de anterior u
%   f_predicted = valor de f(t, u) utilizando preditor
%   h = passo da malha
%   s = quantidade de estágios utilizados
% OUTPUTS:
%   u = próximo valor de u

  % Tabela de coeficientes de Adams-Moulton.
  % Variável persistente que é declarada uma vez.
  persistent C = [
    [1, 0, 0, 0, 0, 0, 0, 0];
    [1/2, 1/2, 0, 0, 0, 0, 0, 0];
    [5/12, 2/3, -1/12, 0, 0, 0, 0, 0];
    [3/8, 19/24, -5/24, 1/24, 0, 0, 0, 0];
    [251/720, 323/360, -11/30, 53/360, -19/720, 0, 0, 0];
    [95/288, 1427/1440, -133/240, 241/720, -173/1440, 3/160, 0, 0];
    [19087/60480, 2713/2520, -15487/20160, 586/945, -6737/20160, 263/2520, -863/60480, 0];
    [5257/17280, 139849/120960, -4511/4480, 123133/120960, -88547/120960, 1537/4480, -11351/120960, 275/24192]
  ];
  
  % Último valor conhecido.
  u = u0;
  % Somando primeiro termo (utiliza preditor).
  u += h*C(s, 1)*f_predicted;
  % Somando termos restantes.
  u += h*dot(C(s, 2:s), fliplr(f0));
endfunction