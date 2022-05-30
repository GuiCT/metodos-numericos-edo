function [u, f0] = RungeKuttaRalstonStep(f, t0, u0, h)
  persistent A = [0,          0,            0,          0;
                  0.4,        0,            0,          0;
                  0.29697761, 0.15875964,   0,          0;
                  0.21810040, -3.05096516,  3.83286476, 0];
  persistent B = [0.17476028, -0.55148066, 1.20553560, 0.17118478];
  persistent C = [0, 0.4, 0.45573725, 1];
  k_size = max(size(B));
  k = zeros(k_size, max(size(u0)));
  for i=1:k_size
    k(i, :) = f(t0 + h*C(i), u0 + h*sum(A(i,:)'.*k, 1));
  endfor
  u = u0 + h*sum(B'.*k, 1);
  f0 = k(1,:);
endfunction