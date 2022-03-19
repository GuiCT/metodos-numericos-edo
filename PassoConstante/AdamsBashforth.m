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
  % Guarda clock do momento antes do cálculo.
  clock0 = clock();
  
  % Descobrindo os primeiros s pontos através de métodos de passo simples.
  % Método de passo simples utilizado foi o Runge-Kutta de Ralston.
  for i = 2 : min(s, n)
    % Todos as linhas da coluna i são calculadas de uma vez.
    u(:, i) = RungeKuttaRalstonStep(f, t(i-1), u(:, i-1), h);
    % A função f é chamada quatro vezes no método RK-Ralston.
    info.nEvals += 4;
  endfor
  
  % A partir do s+1 ponto, utiliza Adams-Bashforth
  for i = s+1 : n
    % São passadas mais de uma coluna para a função AdamsBashforthStep,
    % visto que o método é de passo múltiplo e requer mais de um ponto
    % para todas as equações no sistema.
    % No caso, são passados os s pontos anteriores, visto que
    % são utilizados s estágios.
    u(:, i) = AdamsBashforthStep(f, t(i-s : i-1), u(:, i-s : i-1), h, s);
    % A função f é chamada s vezes quando há s estágios presentes.
    info.nEvals += s;
  endfor
  % Calculando tempo total de execução.
  % Diferença entre o momento atual dado por clock() e
  % clock0 registrado anteriormente.
  % O valor é dado em segundos pela função etime e
  % transformado em milissegundos ao multiplicar por 1000.
  info.tElapsed = etime(clock(), clock0)*1000;
endfunction
function u = AdamsBashforthStep(f, t0, u0, h, s)
% u = AdamsBashforthStep(f, t0, u0, h, s)
% Executa um passo do método numérico de Adams-Bashforth.
% INPUTS:
%   f = função du/dt = f(t, u)
%   t0 = valores anteriores de t
%   u0 = valores anteriores de u
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
  
  % Tamanho do vetor de pontos conhecidos
  n = size(u0)(2);
  % Último valor conhecido
  u = u0(:, n);
  % Valor auxiliar
  aux = n + 1;
  
  for i = 1 : s
    % Adiciona termo a termo do método de Adams-Bashforth.
    u += h*C(s, i)*f(t0(aux-i), u0(:, aux-i));
  endfor
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