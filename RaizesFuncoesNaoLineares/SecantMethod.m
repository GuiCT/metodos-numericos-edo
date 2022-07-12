function [z_h, f_h] = SecantMethod(f, z0, z1, tol, maxitr)
% Executa método da Secante para uma determinada função f.
% Guarda valores de z e f(z) para cada z.
% Útil para mostrar a convergência de cada passo.
f0 = f(z0);
f1 = f(z1);
z_h = [z0;z1];
f_h = [f0;f1];
itr = 0;
while itr < maxitr
    % Encontra novo z a partir do Método da Secante
    z2 = secant_method_step(f0, f1, z0, z1);
    % Se z2 for NaN (ocorre quando f(z1) = f(z0)), encerra.
    % Não garante a convergência.
    if isnan(z2)
        break
    end
    % Encontra novo f(z)
    f2 = f(z2);
    % Adiciona valores ao histórico.
    z_h = [z_h;z2];
    f_h = [f_h;f2];
    % Se f2 está suficientemente perto de zero dada a tolerância
    % especificada, encerra.
    if abs(f2) < tol
        break
    end
    % Atualizando z0, z1, f0 e f1.
    z0 = z1;
    f0 = f1;
    z1 = z2;
    f1 = f2;
    % Próxima iteração.
    itr = itr + 1;
end
end