function sol = multi_step(f, tspan, u0, h, method_simple, method_multi, s)
% sol = multi_step(f, tspan, u0, h)
% Executa um método numérico de passo múltiplo.
% INPUTS:
%   f = função du/dt = f(t, u).
%   tspan = intervalo do domínio em t [t1, t2]
%   u0 = valor inicial de u [u0_1, u0_2, ..., u0_m]
%   h = passo da malha: se não atingir o final, será adaptado para um valor que o faça.
%   method_simple = string indicando o método de passo simples utilizado para
%   INICIALIZAR o método de passo múltiplo:
%     'Euler': Euler explícito
%     'Heun': Método de Heun
%     'Midpoint': Método do ponto médio
%     'RK4': Método de Runge-Kutta clássico
%     'RKR': Método de Runge-Kutta de Ralston (4a ordem)
%   Caso a string não seja nenhuma dessas, será utilizado o método do ponto médio.
%   method_multi = string indicando o método de passo múltiplo utilizado
%     'AB': Adams-Bashforth
%     'AM': Adams-Moulton
%   Caso a string não seja nenhuma dessas, será utilizado o método de Adams-Moulton.
%   s = quantidade de estágios no método do passo múltiplo (4..8).
% OUTPUTS:
%   sol
%     .t = vetor t contendo os valores da variável independente
%     .u = vetor/matriz u contendo os valores de u para cada t.
%     .fevals = número de vezes que a função f foi chamada.
%     .telapsed = tempo de execução total em milissegundos.

% Verifica número de estágios (deve estar entre 4..8)
if (s < 4 || s > 8)
    warning('Número de estágios deve estar entre 4 e 8.');
    return;
end
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
% Armazenando últimos s valores de f(t, u).
last_f = zeros(s, m); last_f(1, :) = f(sol.t(1), sol.u(1,:));
% Criando handler para a função que calcula o passo simples
switch(method_simple)
    case 'Euler'
        step_simple = @ExplicitEulerStep;
        fevals_per_step_simple = 1;
    case 'Heun'
        step_simple = @HeunMethodStep;
        fevals_per_step_simple = 2;
    case 'Midpoint'
        step_simple = @MidpointMethodStep;
        fevals_per_step_simple = 2;
    case 'RK4'
        step_simple = @RungeKutta4Step;
        fevals_per_step_simple = 4;
    case 'RKR'
        step_simple = @RungeKuttaRalstonStep;
        fevals_per_step_simple = 4;
    otherwise
        step_simple = @MidpointMethodStep;
        fevals_per_step_simple = 2;
end
% Criando handler para a função que calcula o passo múltiplo.
switch(method_multi)
    case 'AB'
        step_multi = @AdamsBashforthStep;
        fevals_per_step_multi = 1;
    otherwise
        step_multi = @AdamsMoultonStep;
        fevals_per_step_multi = 2;
end
% Guarda clock do momento antes do cálculo.
clock1 = clock();

% Descobrindo os primeiros s pontos através de métodos de passo simples.
for i = 2 : min(s, n)
    % Todos as linhas da coluna i são calculadas de uma vez.
    sol.u(i,:) = step_simple(f, sol.t(i-1), sol.u(i-1,:), h);
    % Armazenando valor de f(t,u).
    last_f(i, :) = f(sol.t(i), sol.u(i,:));
    % A função f é chamada várias vezes, a depender do método.
    % Uma a mais para calcular e armazenar f(t,u).
    sol.fevals = sol.fevals + fevals_per_step_simple + 1;
end

% A partir do s+1 ponto, utiliza método de passo múltiplo.
for i = s+1 : n
    % Caso estiver utilizando Adams-Moulton, é necessário calcular um preditor
    % e guardá-lo na última posição de last_f.
    if strcmp(method_multi, 'AM')
        last_f = circshift(last_f, -1);
        sol.u(i,:) = AdamsBashforthStep(last_f, sol.u(i-1,:), h, s);
        last_f(end,:) = f(sol.t(i), sol.u(i,:));
    end
    % São passadas mais de uma linha para a função de passo múltiplo,
    % visto que o método é de passo múltiplo e requer mais de um ponto
    % para todas as equações no sistema.
    % No caso, são passados os s valores anteriores de f(t, u).
    sol.u(i,:) = step_multi(last_f, sol.u(i-1,:), h, s);
    if strcmp(method_multi, 'AB')
        % Deslocar o vetor last_f uma unidade para trás.
        % Em Adams-Moulton isso não é necessário pois o último valor de f
        % é o preditor, que será corrigido.
        last_f = circshift(last_f, -1);
    end
    % Armazenar último valor de f(t,u)
    last_f(end,:) = f(sol.t(i), sol.u(i,:));
    % A função f é várias vezes dependente do método.
    sol.fevals = sol.fevals + fevals_per_step_multi;
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