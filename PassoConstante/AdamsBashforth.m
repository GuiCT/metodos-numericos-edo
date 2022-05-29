function [u, info] = AdamsBashforth(f, tSpan, u0, n, s)
% [u, info] = AdamsBashforth(f, tSpan, u0, n)
% Executa o método numérico de Adams-Bashforth.
% INPUTS:
%   f = função du/dt = f(t, u)
%   tSpan = intervalo do domínio em t
%   u0 = valor inicial de u
%   n = número de divisões da malha, quanto maior, mais refinada a malha.
%   s = número de estágios utilizados (4..8)
% OUTPUTS:
%   u = [X, n] = vetor contendo valores da função u em todo o domínio em t.
%   X indica a quantidade de EDOs solucionadas.
%     X = 1 sendo uma única EDO.
%     X = k sendo um sistema de EDOs com k equações.
%   info
%     .nEvals = número de vezes que a função f foi calculada.
%     .tElapsed = tempo de execução total do método (em milissegundos).
%     (não inclui etapas preparatórias, tais como as alocações).
  
  % Verificando se o número de estágios é válido.
  if s < 4 || s > 8
    printf("Número de estágios deve estar entre 4 e 8.");
    return;
  endif

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
  % Valores anteriores de f(t, u)
  f0 = zeros(X, s);
  % Guarda clock do momento antes do cálculo.
  clock0 = clock();
  
  % Descobrindo os primeiros s pontos através de métodos de passo simples.
  % Método de passo simples utilizado foi o Runge-Kutta de Ralston.
  for i = 2 : min(s, n)
    % Todos as linhas da coluna i são calculadas de uma vez.
    [u(:, i), f0(:, i-1)] = RungeKuttaRalstonStep(f, t(i-1), u(:, i-1), h);
    % A função f é chamada quatro vezes no método RK-Ralston.
    info.nEvals += 4;
  endfor
  
  % A partir do s+1 ponto, utiliza Adams-Bashforth
  for i = s+1 : n
    % Se NÃO for primeiro passo múltiplo
    if i != s+1
      % Deslocar cada valor uma unidade para trás.
      f0 = shift(f0, -1);
    endif
    f0(:, end) = f(t(i-1), u(:, i-1));
    % São passadas mais de uma coluna para a função AdamsBashforthStep,
    % visto que o método é de passo múltiplo e requer mais de um ponto
    % para todas as equações no sistema.
    % No caso, são passados os s valores anteriores de f(t, u).
    u(:, i) = AdamsBashforthStep(f, f0, u(:, i-1), h, s);
    % A função f é chamada 1 vez.
    info.nEvals += 1;
  endfor
  % Calculando tempo total de execução.
  % Diferença entre o momento atual dado por clock() e
  % clock0 registrado anteriormente.
  % O valor é dado em segundos pela função etime e
  % transformado em milissegundos ao multiplicar por 1000.
  info.tElapsed = etime(clock(), clock0)*1000;
endfunction