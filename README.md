# Métodos Numéricos para resolução de Equações Diferenciais Ordinárias

Esse repositório visa organizar minhas implementações de métodos numéricos utilizados na resolução de EDOs de 1ª Ordem e/ou Sistema de EDOs de 1ª Ordem. Atualmente, todas as implementações foram feitas utilizando a linguagem GNU Octave/MATLAB, mas é minha intenção reescrever os mesmos códigos utilizando a linguagem C++. Abaixo está uma lista de todos os métodos implementados:

- Método de Euler (explícito)
- Método de Heun
- Método do Ponto Médio
- Método de Runge-Kutta de Quarta Ordem (versão original, RK4)
- Método de Runge-Kutta Ralston de Quarta Ordem*
- Método de Adams-Bashforth
- Método de Adams-Moulton

> *Este método foi apresentado por Anthony Ralston em seu artigo "Runge-Kutta Methods with Minimum Error Bounds", o mesmo pode ser encontrado em https://www.ams.org/journals/mcom/1962-16-080/S0025-5718-1962-0150954-0/

Além de implementar esses métodos, o repositório contém um script que permite automaticamente comparar os métodos, apresentando:
- Tempo de execução (em milissegundos)
- Número de chamadas a função ``f(u, t)``
- Erro Local de Truncamento para cada passo dos métodos (caso seja apresentada uma solução exata a ser utilizada para comparação, visto que o script ainda não é capaz de estimar tal erro)
- Erro Global de Truncamento, a soma de todos os ELT (sob as mesmas condições)

Dentre os futuros objetivos deste repositório estão:
- Implementar métodos com incremento variável (malhas irregulares) visando melhorar a eficiência do código
- Permitir a implementação de qualquer método da família de Métodos de Runge-Kutta a partir de uma Matriz de Butcher**
- Implementar o Método de Bulirsch-Stoer***

> **Uma matriz utilizada para armazenar os coeficientes de um Método Runge-Kutta, formulada por John C. Butcher

> ***Método que utiliza do Método do Ponto Médio Modificado e da Extrapolação de Richardson e se apresenta como um método mais eficiente em baixas tolerâncias. Formulado por Roland Bulirsch e Josef Stoer.
