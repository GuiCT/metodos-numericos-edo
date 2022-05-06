/**
 * @file Secant.hpp
 * @author Guilherme Cesar Tomiasi (gtomiasi@gmail.com)
 * @brief M�todo da Secante
 * @date 2022-05-05
 */

#include <functional>

/*
* @brief Rotina utilizada para encontrar a ra�z da fun��o de uma vari�vel
* Fun��o passada deve obedecer ao prot�tipo double(double)
* @param[in] fun Fun��o cuja ra�z deseja-se encontrar (entrada)
* @param[in] x0 Valor utilizado para primeira itera��o (entrada) 
* @param[in] x1 Valor utilizado para primeira itera��o (entrada)
* @param[in] tolerance Toler�ncia utilizada na execu��o do m�todo (entrada)
* @param[in] maxIterations N�mero m�ximo de itera��es realizadas (entrada)
*/
double secant(
	std::function<double(double)>& fun,
	double x0,
	double x1,
	double tolerance,
	std::size_t maxIterations
);
