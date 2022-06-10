function sol = simple_step(f, tspan, u0, h, method)
% sol = simple_step(f, tspan, u0, h)
% Executa um método numérico de passo simples.
% INPUTS:
%   f = função du/dt = f(t, u).
%   tspan = intervalo do domínio em t [t1, t2]
%   u0 = valor inicial de u [u0_1, u0_2, ..., u0_m]
%   h = passo da malha: se não atingir o final, será adaptado para um valor que o faça.
%   method = string indicando o método de passo simples utilizado:
%     'Euler': Euler explícito
%     'Heun': Método de Heun
%     'Midpoint': Método do ponto médio
%     'RK4': Método de Runge-Kutta clássico
%     'RKR': Método de Runge-Kutta de Ralston (4a ordem)
%   Caso a string não seja nenhuma dessas, será utilizado o método do ponto médio.
% OUTPUTS:
%   sol
%     .t = vetor t contendo os valores da variável independente
%     .u = vetor/matriz u contendo os valores de u para cada t.
%     .fevals = número de vezes que a função f foi chamada.
%     .telapsed = tempo de execução total em milissegundos.

  % Quantidade de EDOs no sistema, inferindo a partir de u0.
  % (considera-se que u0 terá as mesmas dimensões de um resultado
  % retornado por f).
  m = max(size(u0));
  n = (tspan(2) - tspan(1))/h;
  % Arredondando pra caso não seja um inteiro.
  % Soma um para obter o passo correto.
  n = round(n) + 1;
  % Calculando novo passo da malha, garantindo que atinge o final.
  h = (tspan(2) - tspan(1))/n;
  % Alocação do resultado u.
  sol.u = zeros(n, m); sol.u(1,:) = u0(:)';
  % Alocação do domínio discreto.
  sol.t = linspace(tspan(1), tspan(2), n)';
  sol.fevals = 0;
  % Criando handler para a função que calcula o passo, dependente do método
  switch(method)
    case 'Euler'
      step = @ExplicitEulerStep;
      fevals_per_step = 1;
    case 'Heun'
      step = @HeunMethodStep;
      fevals_per_step = 2;
    case 'Midpoint'
      step = @MidpointMethodStep;
      fevals_per_step = 2;
    case 'RK4'
      step = @RungeKutta4Step;
      fevals_per_step = 4;
    case 'RKR'
      step = @RungeKuttaRalston;
      fevals_per_step = 4;
    otherwise
      step = @MidpointMethodStep;
      fevals_per_step = 2;
  end
  % Guarda clock do momento antes do cálculo.
  clock1 = clock();
  
  for i = 2 : n
    % Todos as linhas da coluna i são calculadas de uma vez.
    sol.u(i, :) = step(f, sol.t(i-1), sol.u(i-1,:), h);
    % Atualiza número de chamadas a função.
    sol.fevals = sol.fevals + fevals_per_step;
  end
  % Guarda clock do momento após o cálculo.
  clock2 = clock();
  % Calculando tempo total de execução.
  % Diferença entre o momento atual dado por clock() e
  % clock0 registrado anteriormente.
  % O valor é dado em segundos pela função etime e
  % transformado em milissegundos ao multiplicar por 1000.
  sol.telapsed = etime(clock2, clock1)*1000;
end