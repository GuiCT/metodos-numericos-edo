function sol = AdamsBashforth(f, tspan, u0, h, method_simple, s)
% sol = AdamsBashforth(f, tspan, u0, h)
% Executa um método numérico de AdamsBashforth de s estágios.
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
%   s = quantidade de estágios no método do passo múltiplo (4..8).
% OUTPUTS:
%   sol
%     .t = vetor t contendo os valores da variável independente
%     .u = vetor/matriz u contendo os valores de u para cada t.
%     .fevals = número de vezes que a função f foi chamada.
%     .telapsed = tempo de execução total em milissegundos.
sol = multi_step(f, tspan, u0, h, method_simple, 'AB', s);
endfunction