function new_h = adapt_stepsize(old_h, tol, error)
persistent S p_grow p_shrink;

if isempty(S)
    S = 0.9;
    p_grow = 0.2;
    p_shrink = 0.25;
end

if tol >= error
    p = p_shrink;
else
    p = p_grow;
end

new_h = S*old_h*(tol/error).^p;

end