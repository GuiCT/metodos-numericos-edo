function z = secant_method_step(equation, z0, z1)
if z1 < z0
    aux = z1;
    z1 = z0;
    z0 = aux;
end
numerical_derivative = (z1-z0)/(equation(z1)-equation(z0));
z = z1 - equation(z1)*numerical_derivative;
end

