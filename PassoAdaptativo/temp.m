f = @(t, u) u;
tspan = [0, 12];
u0 = 1;
h = 0.4;
f_abs = @(t) exp(t);

m = max(size(u0));
t = tspan(1):h:tspan(2);
n = max(size(t));
u = zeros(n, m); u(1, :) = u0;
e = zeros(n, m);
for i = 2 : n
    [temp_1, temp_2] = CashKarpStep(f, t(i-1), u(i-1), h);
    u(i, :) = temp_1;
    e(i, :) = temp_2;
end

u_abs = f_abs(t');
e_abs = abs(u_abs - u);