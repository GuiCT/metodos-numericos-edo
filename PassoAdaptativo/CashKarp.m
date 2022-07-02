function [t, u] = CashKarp(f, tspan, u0, first_step, tol)
% Executa o método numérico de Cash-Karp, com passo adaptativo
% Escrito com base na obra "Numerical Recipes in C"

% Quantidade de equações no sistema
m = max(size(u0));
% Pré-aloca vetores que irão conter os valores de t e u
% Para tal, utiliza o primeiro passo (first_step) como base
array_size = diff(tspan)/first_step;
t = zeros(array_size, 1);
u = zeros(array_size, m);
u(1, :) = u0;
t(1) = tspan(1);
% Tamanho do passo atual e valor de t
current_h = first_step;
current_t = tspan(1);
% Indíce da posição atual nos vetores t e u
i = 2;
while(current_t < tspan(2))
    % Se o t passa do valor final do intervalo
    if (current_t + current_h) > tspan(2)
        % Corrige h para terminar em t(2)
        current_h = tspan(2) - current_t;
    end
    [temporary_u, e, ~] = cash_karp_step(f, current_t, u(i-1, :), current_h);
    % Se alguma equação do sistema se aproximar de 0
    if (min(abs(temporary_u))) < 1e-1
        % Utilizar erro absoluto
        max_error = max(abs(e));
    else
        % Caso contrário, utilizar erro relativo
        max_error = max(abs(e./temporary_u));
    end
    if max_error > tol
        % Necessário realizar o passo novamente
        % Atualizando h
        current_h = adapt_stepsize(current_h, tol, max_error);
    else
        % Passo aceito
        t(i) = current_t + current_h;
        u(i, :) = temporary_u;
        i = i + 1;
        % Se passar do tamanho alocado
        if i > array_size
            % Alocando mais 30 entradas no vetor
            t(i:i+29, 1) = zeros(30, 1);
            u(i:i+29, :) = zeros(30, m); 
            array_size = array_size + 30;
        end
        % Atualizando h e t
        current_t = current_t + current_h;
        current_h = adapt_stepsize(current_h, tol, max_error);
    end
end
% Se os vetores t e u não ultrapassarem o tamanho pré-alocado, ignorar os
% valores adiante.
if array_size > i-1
    t = t(1:i-1, 1);
    u = u(1:i-1, :);
end
end
