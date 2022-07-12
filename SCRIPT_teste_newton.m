clc
clear
format long

[z,fz] = NewtonMethod(@integrateBlasius, 0.3, 1e-15, 100);
csv = [[1:max(size(z))]' z fz];

% Função R(x), que calcula R(z) e g'(b)
function [residual, derivative] = integrateBlasius(z)
% Sistema de equações estendido
f = @(t,u) [u(2),u(3),-0.5*u(1)*u(3),u(5),u(6),-0.5*(u(6)*u(1)+u(3)*u(4))];
tspan = [0,10];
u0 = [0,0,z,0,0,1];
first_step = 0.1;
tol = 1e-15;
[~, u] = CashKarp(f,tspan,u0,first_step,tol);
residual = u(end,2) - 1;
derivative = u(end,5);
end