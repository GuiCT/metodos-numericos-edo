Nesta pasta estão presentes métodos numéricos para a solução de Problemas de Valor Inicial (PVI´s). A lista abaixo apresenta todos os métodos implementados.

- Método de Euler (explícito)
- Método de Heun
- Método do Ponto Médio
- Método de Runge-Kutta de Quarta Ordem (versão original, RK4)
- Método de Runge-Kutta Ralston de Quarta Ordem*
- Método de Adams-Bashforth
- Método de Adams-Moulton

> *Este método foi apresentado por Anthony Ralston em seu artigo "Runge-Kutta Methods with Minimum Error Bounds", o mesmo pode ser encontrado em https://www.ams.org/journals/mcom/1962-16-080/S0025-5718-1962-0150954-0/

Além de implementar esses métodos, a pasta contém um script (``TestAllMethods.m``) que permite automaticamente comparar os métodos, apresentando:
- Tempo de execução (em milissegundos)
- Número de chamadas a função ``f(t, u)``
- Erro Local de Truncamento para cada passo dos métodos (caso seja apresentada uma solução exata a ser utilizada para comparação, visto que o script ainda não é capaz de estimar tal erro)
- Erro Global de Truncamento, a soma de todos os ELT (sob as mesmas condições)

Dentre os futuros objetivos desta pasta estão:
- Permitir a implementação de qualquer método da família de Métodos de Runge-Kutta a partir de uma Matriz de Butcher**

> **Uma matriz utilizada para armazenar os coeficientes de um Método Runge-Kutta, formulada por John C. Butcher

Os métodos presentes nesta pasta utilizam um passo constante, utilizando de um domínio discreto regular. Planejo implementar métodos numéricos adaptativos em outra pasta, futuramente.

---
OBS: Alguns scripts escritos aqui utilizam de uma API antiga, e não irão funcionar. É intenção minha atualizar esses scripts no futuro.
