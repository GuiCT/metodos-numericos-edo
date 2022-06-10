function sol = ExplicitEuler(f, tspan, u0, h)
% sol = ExplicitEuler(f, tspan, u0, h)
% Executa um método numérico de Euler explícito.
% INPUTS:
%   f = função du/dt = f(t, u).
%   tspan = intervalo do domínio em t [t1, t2]
%   u0 = valor inicial de u [u0_1, u0_2, ..., u0_m]
%   h = passo da malha: se não atingir o final, será adaptado para um valor que o faça.
% OUTPUTS:
%   sol
%     .t = vetor t contendo os valores da variável independente
%     .u = vetor/matriz u contendo os valores de u para cada t.
%     .fevals = número de vezes que a função f foi chamada.
%     .telapsed = tempo de execução total em milissegundos.
sol = simple_step(f, tspan, u0, h, 'Euler');
end