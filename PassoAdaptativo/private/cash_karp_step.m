function [u, e, k1] = cash_karp_step(f, t0, u0, h)
%C Realiza um passo de integração do método numérico de Cash-Karp
%  O Método retorna um  valor de quinta ordem e uma estimativa de erro
%  de quarta ordem, respectivamente.
persistent A B B_DIFF C;

if isempty(A)
    A = [0,         0,      0,          0,              0           0;
        1/5,        0,      0,          0,              0           0;
        3/40,       9/40,   0,          0,              0           0;
        3/10,       -9/10,  6/5,        0,              0           0;
        -11/54,     5/2,    -70/27,     35/27,          0           0;
        1631/55296, 175/512, 575/13824,  44275/110592,   253/4096,   0];
    B = [37/378, 0, 250/621, 125/594, 0, 512/1771];
    % Segunda linha de coeficientes B
    B_prime = [2825/27648, 0, 18575/48384, 13525/55296, 277/14336, 1/4];
    % Diferença entre a segunda linha dos coeficientes B e a primeira.
    B_DIFF = B - B_prime;
    C = [0, 1/5, 3/10, 3/5, 1, 7/8];
end

k_size = 6;
k = zeros(k_size, max(size(u0)));
for i=1:k_size
    k(i, :) = f(t0 + h*C(i), u0 + h*sum(A(i,:)'.*k, 1));
end
u = u0 + h*sum(B'.*k, 1);
% e = (u0 + h*(B*k)) - (u0 + h*(B_prime*K))
% e = h*((B-B_prime)*k)
% e = h*B_DIFF*k
e = h*sum(B_DIFF'.*k, 1);
k1 = k(1, :);
end

