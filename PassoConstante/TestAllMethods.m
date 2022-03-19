function TestAllMethods(f, tSpan, u0, n, s, file_name, sol)
% TestAllMethods(f, tSpan, u0, n, s, file_name, sol)
% Testa todos os métodos numéricos ao mesmo tempo,
% plota seus devidos gráficos, erro local de truncamento
% e erro acumulado, assim como seus tempos de execução
% e quantidade de chamadas a função.
% INPUTS:
%   f = função du/dt = f(t, u)
%   tSpan = intervalo do domínio em t
%   u0 = valor inicial de u
%   n = número de divisões da malha, quanto maior, mais refinada a malha.
%   s = número de estágios para os Métodos de Adams-Bashforth e Adams-Moulton.
%   file_name = nome do arquivo para armazenar os resultados obtidos.
%   sol = vetor contendo a solução exata do caso apresentado.
%   O parâmetro sol pode ser omitido, mas sem o mesmo não será
%   possível calcular os erros locais de truncamento e erros globais.
% OUTPUTS:
%   Não há.

  % Quantidade de EDOs no sistema.
  X = size(u0)(1);
    
  % Aplicando os métodos
  [ee, i_ee] = ExplicitEuler(f, tSpan, u0, n);
  [hm, i_hm] = HeunMethod(f, tSpan, u0, n);
  [mm, i_mm] = MidpointMethod(f, tSpan, u0, n);
  [rk4, i_rk4] = RungeKutta4(f, tSpan, u0, n);
  [rkr, i_rkr] = RungeKuttaRalston(f, tSpan, u0, n);
  [ab, i_ab] = AdamsBashforth(f, tSpan, u0, n, s);
  [am, i_am] = AdamsMoulton(f, tSpan, u0, n, s);
  t = linspace(tSpan(1), tSpan(2), n);

  % ==========================================================================
  % =========================== PLOTANDO GRÁFICOS ============================
  % ==========================================================================
  % Um gráfico para cada equação no sistema
  for i = 1:X
    figure(i);
    plot(
      t, ee(i, :), 'k-',
      t, hm(i, :), 'r-',
      t, mm(i, :), 'g-',
      t, rk4(i, :),'b-',
      t, rkr(i, :),'y-',
      t, ab(i, :), 'm-',
      t, am(i, :), 'c-'
    );
    xlabel("Valor de t");
    ylabel("Valor de u");
    legend(
      "Euler Explícito",
      "Método de Heun",
      "Método do ponto médio",
      "Runge-Kutta 4 (RK4)",
      "Runge-Kutta Ralston",
      "Adams-Bashforth",
      "Adams-Moulton"
    );
    title(["Equação " num2str(i) " do sistema."]);
  endfor
  % ==========================================================================
  % ================ CALCULANDO ELT E PLOTANDO (SE POSSÍVEL) =================
  % ==========================================================================
  % Se não existir o parâmetro que dá a solução exata,
  % os vetores de erro não serão alocados.
  if exist("sol")
    % Erro local de truncamento de Euler Explícito.
    elt_ee = ee - sol;
    % Erro local de truncamento do Método de Heun.
    elt_hm = hm - sol;
    elt_mm = mm - sol;
    elt_rk4 = rk4 - sol;
    elt_rkr = rkr - sol;
    elt_ab = ab - sol;
    elt_am = am - sol;
    % Calculando erros globais de truncamento
    egt_ee = sum(elt_ee, 2);
    egt_hm = sum(elt_hm, 2);
    egt_mm = sum(elt_mm, 2);
    egt_rk4 = sum(elt_rk4, 2);
    egt_rkr = sum(elt_rkr, 2);
    egt_ab = sum(elt_ab, 2);
    egt_am = sum(elt_am, 2);
    
    % Erro Local de Truncamento para todas as equações do sistema
    for i = 1:X
      figure(X + i);
      plot(
        t, elt_ee(1, :), 'k-',
        t, elt_hm(1, :), 'r-',
        t, elt_mm(1, :), 'g-',
        t, elt_rk4(1, :),'b-',
        t, elt_rkr(1, :),'y-',
        t, elt_ab(1, :), 'm-',
        t, elt_am(1, :), 'c-'
      );
      legend(
        "ELT de Euler Explícito",
        "ELT de Método de Heun",
        "ELT de Método do ponto médio",
        "ELT de Runge-Kutta 4 (RK4)",
        "ELT de Runge-Kutta Ralston",
        "ELT de Adams-Bashforth",
        "ELT de Adams-Moulton"
      );
      title(["Erros Locais de Truncamento, equação " num2str(i)]);
    endfor
  endif
  % ==========================================================================
  % ============================= FIM DOS PLOTS ==============================
  % ==========================================================================
  % Criando o arquivo destino.
  file = fopen(file_name, "w");
  % Função que imprime na tela e salva no arquivo ao mesmo tempo.
  function print_and_save(s)
    % Imprime na tela uma string s
    printf(s);
    % Salva no arquivo uma string s
    fprintf(file, s);
  endfunction
  % ==========================================================================
  % ===================== CABEÇALHO DA IMPRESSÃO/ARQUIVO =====================
  % ==========================================================================
  print_and_save("Metodo,Numero de chamadas a f,Tempo decorrido(ms)");
  % Caso exista EGT disponível
  if exist("sol")
    % Para todas as equações do sistema, criando cabeçalho para 
    % Erro Global de Trucamento (EGT)
    for i = 1:X
      % Caso seja a última equação
      if i == X
        print_and_save(["EGT da Equacao " num2str(i) "\n"]);
      else
        print_and_save([",EGT da Equacao " num2str(i)]);
      endif
    endfor
  else
    % Encerra o cabeçalho
    print_and_save("\n");
  endif
  % ==========================================================================
  % ============================ FIM DO CABEÇALHO ============================
  % ==========================================================================
  % Formata uma linha da saída
  % Parâmetros s e egt são opcionais
  function print_line(method_name, nEvals, tElapsed, s, useS, egt)
    if useS
      print_and_save([method_name " de " num2str(s) " estagios," num2str(nEvals) "," num2str(tElapsed)]);
    else
      print_and_save([method_name "," num2str(nEvals) "," num2str(tElapsed)]);
    endif
    
    if exist("egt")
      for i = 1:X
        if i == X
          print_and_save([num2str(egt(i)) "\n"]);
        else
          print_and_save(["," num2str(egt(i))]);
        endif
      endfor
    else
      print_and_save("\n");
    endif
  endfunction
  % ==========================================================================
  % ============================ IMPRIMINDO LINHAS ===========================
  % ==========================================================================
  if exist("sol")
    % Incluindo erro global de truncamento
    print_line("Euler Explicito", i_ee.nEvals, i_ee.tElapsed, s, false, egt_ee);
    print_line("Metodo de Heun", i_hm.nEvals, i_hm.tElapsed, s, false, egt_hm);
    print_line("Metodo do Ponto Medio", i_mm.nEvals, i_mm.tElapsed, s, false, egt_mm);
    print_line("Metodo de Runge-Kutta Quarta Ordem", i_rk4.nEvals, i_rk4.tElapsed, s, false, egt_rk4);
    print_line("Metodo de Runge-Kutta Ralston Quarta Ordem", i_rkr.nEvals, i_rkr.tElapsed, s, false, egt_rkr);
    print_line("Metodo Adams-Bashforth", i_ab.nEvals, i_ab.tElapsed, s, true, egt_ab);
    print_line("Metodo Adams-Moulton", i_am.nEvals, i_am.tElapsed, s, true, egt_am);
  else
    % Não incluindo erro global de truncamento
    print_line("Euler Explicito", i_ee.nEvals, i_ee.tElapsed, s, false);
    print_line("Metodo de Heun", i_hm.nEvals, i_hm.tElapsed, s, false);
    print_line("Metodo do Ponto Medio", i_mm.nEvals, i_mm.tElapsed, s, false);
    print_line("Metodo de Runge-Kutta Quarta Ordem", i_rk4.nEvals, i_rk4.tElapsed, s, false);
    print_line("Metodo de Runge-Kutta Ralston Quarta Ordem", i_rkr.nEvals, i_rkr.tElapsed, s, false);
    print_line("Metodo Adams-Bashforth", i_ab.nEvals, i_ab.tElapsed, s, true);
    print_line("Metodo Adams-Moulton", i_am.nEvals, i_am.tElapsed, s, true);
  endif
  % ==========================================================================
  % ============================ FIM DA IMPRESSÃO ============================
  % ==========================================================================
  fclose(file);
endfunction