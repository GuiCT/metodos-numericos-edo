function u = MidpointMethodStep(f, t0, u0, h)
  persistent A = [0,    0;
                  1/2,  0];
  persistent B = [0, 1];
  persistent C = [0, 1/2];
  k_size = max(size(B));
  k = zeros(k_size, max(size(u0)));
  for i=1:k_size
    k(i, :) = f(t0 + h*C(i), u0 + h*sum(A(i,:)'.*k, 1));
  endfor
  u = u0 + h*sum(B'.*k, 1);
endfunction