clc
clear

% list_of_methods = {'ExplicitEuler','RungeKutta4','AdamsBashforth','AdamsMoulton'};
% is_multistep = [false, false, true, true];
% s = 6;
% auxiliar_method = 'RK4';

f = @(t,u) u;
%f_c = @(t,u) [u(2); u(3); -0.5*u(1)*u(3)];
tspan = [0, 12];
u0 = 1;
%u0_c = u0';
size_u0 = max(size(u0));
h = 0.4;
q = 4;
h_2 = h/q;
%times = zeros(10, 5);
t = tspan(1):h:tspan(2);

% for i=1:10
%     for j=1:4
%         if is_multistep(j)
%             sol = feval(list_of_methods{j}, f, tspan, u0, h, auxiliar_method, 6);
%         else
%             sol = feval(list_of_methods{j}, f, tspan, u0, h);
%         end
%         times(i,j)=sol.telapsed;
%     end
%     tic
%     sol = ode45(f_c, t, u0_c, odeset('RelTol', 1e-15, 'AbsTol', 1e-15));
%     time = toc*1000;
%     times(i, 5) = time;
% end

f_exact = @(t) exp(t);

sol = RungeKutta4(f, tspan, u0, h);
u = sol.u;
sol_2 = RungeKutta4(f, tspan, u0, h_2);
u_2 = sol_2.u;
u_2_c = u_2(1:q:end);
estimated_error = abs((u - u_2_c)./(4.^-5 - 1));
real_error = abs(u - f_exact(sol.t));
error_in_estimate = abs(estimated_error - real_error);
normal_table = [estimated_error, real_error, error_in_estimate];