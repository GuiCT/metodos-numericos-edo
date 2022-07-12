function z = secant_method_step(f0, f1, z0, z1)
numerical_derivative = (f1-f0)/(z1-z0);
z = z1 - f1/numerical_derivative;
end