clc
clear
format long

[z,fz] = SecantMethod(@integrateBlasius, 0.3, 0.4, 1e-15, 100);
csv = [[1:max(size(z))]' z fz];

% Função R(x)
function residual = integrateBlasius(z)
f = @(t,u) [u(2),u(3),-0.5*u(1)*u(3)];
tspan = [0,10];
u0 = [0,0,z];
first_step = 0.1;
tol = 1e-15;
[~, u] = CashKarp(f,tspan,u0,first_step,tol);
residual = u(end,2) - 1;
end