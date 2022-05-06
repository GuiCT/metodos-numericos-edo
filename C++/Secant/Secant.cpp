/**
 * @file Secant.cpp
 * @author Guilherme Cesar Tomiasi (gtomiasi@gmail.com)
 * @brief M�todo da Secante
 * @date 2022-05-05
 */

#include "Secant.hpp"
#include <functional>
#include <cmath>

/*
* Verifica se n0 est� pr�ximo de n1 de acordo com uma toler�ncia (relativa)
* @param[in] n0 Valor a ser comparado (entrada)
* @param[in] n1 Valor a ser comparado (toler�ncia aplicada a este valor, entrada)
* @param[in] tolerance Toler�ncia considerada (entrada)
*/
inline bool inRange(double n0, double n1, double tolerance) {
	return (std::abs(n1 - n0) <= tolerance * std::abs(n1));
}

double secant(
	std::function<double(double)>& fun,
	double x0,
	double x1,
	double tolerance,
	std::size_t maxIterations
) {
	/*
	* Verificando se os dois valores utilizados para o chute inicial
	* s�o iguais. Se sim, o processo � abortado a partir de uma instru��o
	* throw, alertando para o problema.
	*/
	if (x1 == x0) {
		throw ("x1 n�o pode ter o mesmo valor de x0.");
	}
	
	double f0 = fun(x0), f1 = fun(x1);
	double xNew;

	// Verifica se alguma das entradas j� � uma raiz
	if (inRange(f0, 0.0, tolerance))
		return x0;
	if (inRange(f1, 0.0, tolerance))
		return x1;

	/*
	* Verifica se f1 � menor que f0. Se sim, os valores
	* ser�o trocados para que f1 seja maior que f0.
	*/
	if (std::abs(f1) < std::abs(f0)) {
		xNew = x0;
		x0 = x1;
		x1 = xNew;
		xNew = f0;
		f0 = f1;
		f1 = xNew;
	}

	std::size_t iterations;
	for (iterations = 0; iterations < maxIterations; iterations++) {
		/*
		* Se os valores f(x0) e f(x1) convergirem, verifica se
		* x0 e x1 tamb�m convergiram. Se esse for o caso, retorna
		* o x resultante. Caso contr�rio, � executada uma instru��o
		* throw alertando para a n�o-converg�ncia do m�todo.
		*/
		if (f1 == f0) {
			if (x1 != x0)
				throw ("O m�todo falhou em atingir converg�ncia.");
			else
				return x0;
		}
		/*
		* Realiza uma nova itera��o do m�todo.
		* � realizada uma manipula��o alg�brica do passo de itera��o
		* (anteriormente: (x0 * f1 - x1 * f0) / (x1 - x0))
		* De forma que � necess�rio alternar entre
		* x0 / x1 <===> x1 / x0, e
		* x0 / x1 * f1 <===> x1 / x0 * f0
		*/
		if (std::abs(f1) > std::abs(f0))
			xNew = (f0 - x0 / x1 * f1) / (1 - x0 / x1);
		else
			xNew = (f1 - x1 / x0 * f0) / (1 - x1 / x0);
		/*
		* Se o valor de xNew render um resultado muito pr�ximo
		* do valor anterior (x1), retorna xNew.
		*/
		if (inRange(xNew, x1, tolerance))
			return xNew;
		// Atualizando valores para a pr�xima itera��o.
		x0 = x1;
		f0 = f1;
		x1 = xNew;
		f1 = fun(x1);
	}
	// M�todo n�o obteve converg�ncia dentro do n�mero m�ximo de itera��es.
	throw ("O m�todo falhou em atingir converg�ncia.");
}