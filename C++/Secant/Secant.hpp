/**
 * @file Secant.hpp
 * @author Guilherme Cesar Tomiasi (gtomiasi@gmail.com)
 * @brief Método da Secante
 * @date 2022-05-05
 */

#include <functional>

/*
* @brief Rotina utilizada para encontrar a raíz da função de uma variável
* Função passada deve obedecer ao protótipo double(double)
* @param[in] fun Função cuja raíz deseja-se encontrar (entrada)
* @param[in] x0 Valor utilizado para primeira iteração (entrada) 
* @param[in] x1 Valor utilizado para primeira iteração (entrada)
* @param[in] tolerance Tolerância utilizada na execução do método (entrada)
* @param[in] maxIterations Número máximo de iterações realizadas (entrada)
*/
double secant(
	std::function<double(double)>& fun,
	double x0,
	double x1,
	double tolerance,
	std::size_t maxIterations
);
