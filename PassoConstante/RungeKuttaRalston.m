function [u, info] = RungeKuttaRalston(f, tSpan, u0, n)
% [u, info] = RungeKuttaRalston(f, tSpan, u0, n)
% Executa o método numérico de Runge-Kutta de A. Ralston.
% INPUTS:
%   f = função du/dt = f(t, u)
%   tSpan = intervalo do domínio em t
%   u0 = valor inicial de u
%   n = número de divisões da malha, quanto maior, mais refinada a malha.
% OUTPUTS:
%   u = [X, n] = vetor contendo valores da função u em todo o domínio em t.
%   X indica a quantidade de EDOs solucionadas.
%     X = 1 sendo uma única EDO.
%     X = k sendo um sistema de EDOs com k equações.
%   info
%     .nEvals = número de vezes que a função f foi calculada.
%     .tElapsed = tempo de execução total do método (em milissegundos).
%     (não inclui etapas preparatórias, tais como as alocações).

  % Quantidade de EDOs no sistema, inferindo a partir de u0.
  % (considera-se que u0 terá as mesmas dimensões de um resultado
  % retornado por f).
  X = size(u0)(1);
  % Passo da malha (constante).
  h = (tSpan(2) - tSpan(1))/(n - 1);
  % Alocação do u resultante.
  u = zeros(X, n); u(:, 1) = u0;
  % Alocação do domínio discreto.
  t = linspace(tSpan(1), tSpan(2), n);
  % Inicializando número de cálculos de f.
  info.nEvals = 0;
  % Guarda clock do momento antes do cálculo.
  clock0 = clock();
 
  % Calculando os valores restantes para prencheer o vetor.
  % Nesse caso, i indica a coluna do valor a ser calculado,
  % não o valor atual conhecido.
  for i = 2 : n
    % Todos as linhas da coluna i são calculadas de uma vez.
    u(:, i) = RungeKuttaRalstonStep(f, t(i-1), u(:, i-1), h);
    % A função f é chamada quatro vezes no método RK-Ralston.
    info.nEvals += 4;
  endfor
  % Calculando tempo total de execução.
  % Diferença entre o momento atual dado por clock() e
  % clock0 registrado anteriormente.
  % O valor é dado em segundos pela função etime e
  % transformado em milissegundos ao multiplicar por 1000.
  info.tElapsed = etime(clock(), clock0)*1000;
endfunction

function u = RungeKuttaRalstonStep(f, t0, u0, h)
% u = RungeKuttaRalstonStep(f, t0, u0, h)
% Executa um passo do método numérico RK de Ralston.
% INPUTS:
%   f = função du/dt = f(t, u)
%   t0 = valor atual de t
%   u0 = valor atual de u
%   h = passo da malha
% OUTPUTS:
%   u = próximo valor de u

  k1 = f(t0, u0);
  k2 = f(t0 + 0.4*h, u0 + 0.4*h*k1);
  k3 = f(t0 + 0.45573725*h, u0 + 0.29697761*h*k1 + 0.15875964*h*k2);
  k4 = f(t0 + h, u0 + 0.21810040*h*k1 - 3.05096516*h*k2 + 3.83286476*h*k3);
  u = u0 + h*(0.17476028*k1 - 0.55148066*k2 + 1.20553560*k3 + 0.17118478*k4);
endfunction