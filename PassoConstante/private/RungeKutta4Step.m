function u = RungeKutta4Step(f, t0, u0, h)
  persistent A = [0,    0,    0,  0;
                  1/2,  0,    0,  0;
                  0,    1/2,  0,  0;
                  0,    0,    1,  0];
  persistent B = [1/6, 1/3, 1/3, 1/6];
  persistent C = [0, 1/2, 1/2, 1];
  k_size = max(size(B));
  k = zeros(k_size, max(size(u0)));
  for i=1:k_size
    k(i, :) = f(t0 + h*C(i), u0 + h*sum(A(i,:)'.*k, 1));
  endfor
  u = u0 + h*sum(B'.*k, 1);
endfunction