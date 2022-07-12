function [z_h, f_h] = NewtonMethod(f, z0, tol, maxitr)
% Executa método de Newton para uma determinada função f.
% Nesse caso, f calcula não só o valor de f(z), mas de f'(z).
% Guarda valores de z e f(z) para cada z.
% Útil para mostrar a convergência de cada passo.
z_h = z0;
[f0, fprime0] = f(z0);
f_h = f0;
itr = 0;
while itr < maxitr
    % Calcula novo z a partir do Método de Newton
    z1 = newton_method_step(f0, fprime0, z0);
    % Calcula f(z) e f'(z) para o novo z
    [f1, fprime1] = f(z1);
    % Adiciona valores ao histórico.
    z_h = [z_h;z1];
    f_h = [f_h;f1];
    % Se f1 está suficientemente perto de zero dada a tolerância
    % especificada, encerra.
    if abs(f1) < tol
        break
    end
    % Atualizando z0, f0 e fprime0
    z0 = z1;
    f0 = f1;
    fprime0 = fprime1;
    % Próxima iteração
    itr = itr + 1;
end
end