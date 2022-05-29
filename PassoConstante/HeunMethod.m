function [u, info] = HeunMethod(f, tSpan, u0, n)
% [u, info] = HeunMethod(f, tSpan, u0, n)
% Executa o método numérico de Heun
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
    u(:, i) = HeunMethodStep(f, t(i-1), u(:, i-1), h);
    % A função f é chamada duas vezes no método de Heun.
    info.nEvals += 2;
  endfor
  % Calculando tempo total de execução.
  % Diferença entre o momento atual dado por clock() e
  % clock0 registrado anteriormente.
  % O valor é dado em segundos pela função etime e
  % transformado em milissegundos ao multiplicar por 1000.
  info.tElapsed = etime(clock(), clock0)*1000;
endfunction