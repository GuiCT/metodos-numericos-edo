function z = newton_method_step(f0, fprime0, z0)
z = z0 - f0/fprime0;
end