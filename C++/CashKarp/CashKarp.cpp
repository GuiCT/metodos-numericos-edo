/**
* @file CashKarp.cpp
* @author Guilherme Cesar Tomiasi (gtomiasi@gmail.com)
* @brief Algoritmo de Cash-Karp (Runge-Kutta com passo adaptativo)
* @date 2022-04-03
*/

/*
	*  Essa implementa��o � baseada na se��o 16.2 do livro "Numerical Recipes in
	C", escrito pelos autores:
	- William H. Press
	- Saul A. Teukolsky
	- William T. Vetterling
	- Brian P. Flannery
	cujo c�digo ISBN � 0-521-43108-5

	*  As implementa��es deste arquivo, no entanto, n�o utilizam a Linguagem C,
	mas sim a Linguagem C++. N�o interprete essas implementa��es como uma
	simples c�pia dos m�todos apresentados no livro, mas uma adapta��o para
	facilitar a leitura dos m�todos � luz de uma vis�o acad�mica.
*/

#include "CashKarp.hpp"
#include <utility>
#include <vector>
#include <functional>
#include <iostream>
#include <intrin.h>
#include <cmath>

void CashKarp::CashKarpStep(
	std::vector<double>& u,
	std::vector<double>& dudt,
	double t,
	double stepSize,
	std::vector<double>& uOutput,
	std::vector<double>& uError,
	std::function<
	void(
		double,
		std::vector<double>&,
		std::vector<double>&
		)
	>& dynFun)
{
	/*
		Valor dos coeficientes � constante, logo s�o utilizadas
		vari�veis est�ticas.

		Coefficient values are constants, hence the use of
		static variables.
	*/

	/*
		Coeficientes 'c', determinam mudan�a no valor de 't' para
		cada valor intermedi�rio.

		'C' coefficients, determining the change in the value of 't'
		for each intermediate.
	*/
	static double c2 = 1.0 / 5.0,
		c3 = 3.0 / 10.0,
		c4 = 3.0 / 5.0,
		c5 = 1.0,
		c6 = 7.0 / 8.0;

	/*
		Coeficientes 'a', determinam participa��o de cada valor
		intermedi�rio na mudan�a do valor de 'u' para os intermedi�rios
		subsequentes.

		'A' coefficients, determining the weight of intermediate values
		in the change of 'u' when calculating the forthcoming intermediates.
	*/
	static double a21 = 1.0 / 5.0,
		a31 = 3.0 / 40.0, a32 = 9.0 / 40.0,
		a41 = 3.0 / 10.0, a42 = -9.0 / 10.0, a43 = 6.0 / 5.0,
		a51 = -11.0 / 54.0, a52 = 5.0 / 2.0, a53 = -70.0 / 27.0,
		a54 = 35.0 / 27.0,
		a61 = 1631.0 / 55296.0, a62 = 175.0 / 512.0,
		a63 = 575.0 / 13824.0, a64 = 44275.0 / 110592.0,
		a65 = 253.0 / 4096.0;

	/*
		Coeficientes 'b', determinam participa��o de cada valor
		intermedi�rio no c�lculo do valor final de 'u'.

		'B' coefficients, determining the weight of intermediate
		values when calculating the final value of 'u'.
	*/
	static double b1 = 37.0 / 378.0,
		b3 = 250.0 / 621.0,
		b4 = 125.0 / 594.0,
		b6 = 512.0 / 1771.0;

	/*
		Coeficientes 'd', diferen�a entre o coeficiente b do m�todo
		principal e o m�todo embarcado. � utilizado para estimar o erro.

		'D' coefficients, the difference between the 'b' coefficients
		of the main method and the embedded method. It is used to
		estimate the error.
	*/
	static double d1 = -0.0042937748015873,
		d3 = 0.0186685860938579,
		d4 = -0.0341550268308081,
		d5 = -0.0193219866071429,
		d6 = 0.0391022021456804;

	std::size_t uSize = u.size();
	std::size_t i;
	std::vector<double> k2(uSize);
	std::vector<double> k3(uSize);
	std::vector<double> k4(uSize);
	std::vector<double> k5(uSize);
	std::vector<double> k6(uSize);
	std::vector<double> uTemporary(uSize);

	/*
		Calculando valores intermedi�rios k1, k2, ..., k6
		Calculating intermediate values

		Uma itera��o do loop para cada equa��o presente no sistema
		One iteration of the for-loop for each equation in the system
	*/
	for (i = 0; i < uSize; i++)
		uTemporary[i] = u[i] + a21 * stepSize * dudt[i];
	dynFun(t + c2 * stepSize, uTemporary, k2);

	for (i = 0; i < uSize; i++)
		uTemporary[i] = u[i] + stepSize * (a31 * dudt[i] + a32 * k2[i]);
	dynFun(t + c3 * stepSize, uTemporary, k3);

	for (i = 0; i < uSize; i++)
	{
		uTemporary[i] =
			u[i] + stepSize * (a41 * dudt[i] +
				a42 * k2[i] +
				a43 * k3[i]);
	}
	dynFun(t + c4 * stepSize, uTemporary, k4);

	for (i = 0; i < uSize; i++)
	{
		uTemporary[i] =
			u[i] + stepSize * (a51 * dudt[i] +
				a52 * k2[i] +
				a53 * k3[i] +
				a54 * k4[i]);
	}
	dynFun(t + c5 * stepSize, uTemporary, k5);

	for (i = 0; i < uSize; i++)
	{
		uTemporary[i] =
			u[i] + stepSize * (a61 * dudt[i] +
				a62 * k2[i] +
				a63 * k3[i] +
				a64 * k4[i] +
				a65 * k5[i]);
	}
	dynFun(t + c6 * stepSize, uTemporary, k6);

	/*
		Calculando valor na precis�o de quarta ordem
	*/
	for (i = 0; i < uSize; i++)
	{
		uOutput[i] =
			u[i] + stepSize * (b1 * dudt[i] +
				b3 * k3[i] +
				b4 * k4[i] +
				b6 * k6[i]);
	}

	/*
		Estimando erro a partir da diferen�a entre quarta ordem e quinta ordem
		N�o � necess�rio calcular o valor de quinta ordem, visto que
		algebricamente � poss�vel prever qual ser� a diferen�a,
		tendo em vista que os coeficientes da tabela de Butcher s�o constantes.
	*/
	for (i = 0; i < uSize; i++)
	{
		uError[i] =
			stepSize * (d1 * dudt[i] +
				d3 * k3[i] +
				d4 * k4[i] +
				d5 * k5[i] +
				d6 * k6[i]);
	}

	/*
		Fim do m�todo
		End of the function
	*/
}

void CashKarp::CashKarpStepAVX2(
	std::vector<double>& u,
	std::vector<double>& dudt,
	double t,
	double stepSize,
	std::vector<double>& uOutput,
	std::vector<double>& uError,
	std::function<
	void(
		double,
		std::vector<double>&,
		std::vector<double>&
		)
	>& dynFun)
{
	/*
		Valor dos coeficientes � constante, logo s�o utilizadas
		vari�veis est�ticas.

		Coefficient values are constants, hence the use of
		static variables.
	*/

	/*
		Coeficientes 'c', determinam mudan�a no valor de 't' para
		cada valor intermedi�rio.

		'C' coefficients, determining the change in the value of 't'
		for each intermediate.
	*/
	static double c2 = 1.0 / 5.0,
		c3 = 3.0 / 10.0,
		c4 = 3.0 / 5.0,
		c5 = 1.0,
		c6 = 7.0 / 8.0;

	/*
		Coeficientes 'a', determinam participa��o de cada valor
		intermedi�rio na mudan�a do valor de 'u' para os intermedi�rios
		subsequentes.

		'A' coefficients, determining the weight of intermediate values
		in the change of 'u' when calculating the forthcoming intermediates.
	*/
	static double a21 = 1.0 / 5.0,
		a31 = 3.0 / 40.0, a32 = 9.0 / 40.0,
		a41 = 3.0 / 10.0, a42 = -9.0 / 10.0, a43 = 6.0 / 5.0,
		a51 = -11.0 / 54.0, a52 = 5.0 / 2.0, a53 = -70.0 / 27.0,
		a54 = 35.0 / 27.0,
		a61 = 1631.0 / 55296.0, a62 = 175.0 / 512.0,
		a63 = 575.0 / 13824.0, a64 = 44275.0 / 110592.0,
		a65 = 253.0 / 4096.0;

	/*
		Coeficientes 'b', determinam participa��o de cada valor
		intermedi�rio no c�lculo do valor final de 'u'.

		'B' coefficients, determining the weight of intermediate
		values when calculating the final value of 'u'.
	*/
	static double b1 = 37.0 / 378.0,
		b3 = 250.0 / 621.0,
		b4 = 125.0 / 594.0,
		b6 = 512.0 / 1771.0;

	/*
		Coeficientes 'd', diferen�a entre o coeficiente b do m�todo
		principal e o m�todo embarcado. � utilizado para estimar o erro.

		'D' coefficients, the difference between the 'b' coefficients
		of the main method and the embedded method. It is used to
		estimate the error.
	*/
	static double d1 = -0.0042937748015873,
		d3 = 0.0186685860938579,
		d4 = -0.0341550268308081,
		d5 = -0.0193219866071429,
		d6 = 0.0391022021456804;

	std::size_t uSize = u.size();
	std::size_t i;
	std::vector<double> kTemporary(uSize);
	std::vector<double> uTemporary(uSize);

	static int* mask = new int[4]{};
	static double* aux = new double[4]{};
	for (i = 0; i < uSize; i++)
		mask[i] = -1;

	__m256i _mask = _mm256_set_epi64x(mask[3], mask[2], mask[1], mask[0]);
	__m256d _u = _mm256_maskload_pd(u.data(), _mask);
	__m256d _dudt = _mm256_maskload_pd(dudt.data(), _mask);
	__m256d _aux, _aux2, _k2, _k3, _k4, _k5, _k6;

	/*
		Calculando valores intermedi�rios k1, k2, ..., k6
		Calculating intermediate values
	*/

	// k2 = u + a21 * stepSize * dudt;

	// a21 * stepSize * dudt
	aux[0] = a21 * stepSize;
	_aux = _mm256_broadcast_sd(aux);
	_aux = _mm256_mul_pd(_aux, _dudt);
	// u + a21 * stepSize * dudt
	_aux = _mm256_add_pd(_aux, _u);
	_mm256_store_pd(aux, _aux);
	uTemporary.assign(aux, aux + uSize);

	dynFun(t + c2 * stepSize, uTemporary, kTemporary);
	_k2 = _mm256_maskload_pd(kTemporary.data(), _mask);

	// k3 = u + stepSize * (a31 * dudt + a32 * k2)

	// a31 * dudt
	_aux = _mm256_broadcast_sd(&a31);
	_aux = _mm256_mul_pd(_aux, _dudt);

	// a32 * k2
	_aux2 = _mm256_broadcast_sd(&a32);
	_aux2 = _mm256_mul_pd(_aux2, _k2);
	
	// a31 * dudt + a32 * k2
	_aux = _mm256_add_pd(_aux, _aux2);

	// stepSize * (a31 * dudt + a32 * k2)
	_aux2 = _mm256_broadcast_sd(&stepSize);
	_aux = _mm256_mul_pd(_aux, _aux2);

	// u + stepSize * (a31 * dudt + a32 * k2)
	_aux = _mm256_add_pd(_aux, _u);
	_mm256_store_pd(aux, _aux);
	uTemporary.assign(aux, aux + uSize);

	dynFun(t + c3 * stepSize, uTemporary, kTemporary);
	_k3 = _mm256_maskload_pd(kTemporary.data(), _mask);

	/*
		k4 =
			u + stepSize * (a41 * dudt +
								a42 * k2 +
								a43 * k3);
	*/

	// a41 * dudt
	_aux = _mm256_broadcast_sd(&a41);
	_aux = _mm256_mul_pd(_aux, _dudt);

	// a42 * k2
	_aux2 = _mm256_broadcast_sd(&a42);
	_aux2 = _mm256_mul_pd(_aux2, _k2);

	// a41 * dudt + a42 * k2
	_aux = _mm256_add_pd(_aux, _aux2);

	// a43 * k3
	_aux2 = _mm256_broadcast_sd(&a43);
	_aux2 = _mm256_mul_pd(_aux2, _k3);

	// a41 * dudt + a42 * k2 + a43 * k3
	_aux = _mm256_add_pd(_aux, _aux2);

	// stepSize * (a41 * dudt + a42 * k2 + a43 * k3)
	_aux2 = _mm256_broadcast_sd(&stepSize);
	_aux = _mm256_mul_pd(_aux, _aux2);

	// u + stepSize * (a41 * dudt + a42* k2 + a43 * k3);
	_aux = _mm256_add_pd(_aux, _u);
	_mm256_store_pd(aux, _aux);
	uTemporary.assign(aux, aux + uSize);

	dynFun(t + c4 * stepSize, uTemporary, kTemporary);
	_k4 = _mm256_maskload_pd(kTemporary.data(), _mask);

	/*
		k5 =
			u + stepSize * (a51 * dudt +
								a52 * k2 +
								a53 * k3 +
								a54 * k4);
	*/

	// a51 * dudt
	_aux = _mm256_broadcast_sd(&a51);
	_aux = _mm256_mul_pd(_aux, _dudt);

	// a52 * k2
	_aux2 = _mm256_broadcast_sd(&a52);
	_aux2 = _mm256_mul_pd(_aux2, _k2);

	// a51 * dudt + a52 * k2
	_aux = _mm256_add_pd(_aux, _aux2);

	// a53 * k3
	_aux2 = _mm256_broadcast_sd(&a53);
	_aux2 = _mm256_mul_pd(_aux2, _k3);

	// a51 * dudt + a52 * k2 + a53 * k3
	_aux = _mm256_add_pd(_aux, _aux2);

	// a54 * k4
	_aux2 = _mm256_broadcast_sd(&a54);
	_aux2 = _mm256_mul_pd(_aux2, _k4);

	// a51 * dudt + a52 * k2 + a53 * k3 + a54 * k4
	_aux = _mm256_add_pd(_aux, _aux2);

	// stepSize * (a51 * dudt + a52 * k2 + a53 * k3 + a54 * k4)
	_aux2 = _mm256_broadcast_sd(&stepSize);
	_aux = _mm256_mul_pd(_aux, _aux2);

	// u + stepSize * (a51 * dudt + a52 * k2 + a53 * k3 + a54 * k4)
	_aux = _mm256_add_pd(_u, _aux);
	_mm256_store_pd(aux, _aux);
	uTemporary.assign(aux, aux + uSize);

	dynFun(t + c5 * stepSize, uTemporary, kTemporary);
	_k5 = _mm256_maskload_pd(kTemporary.data(), _mask);

	/*
		k6 =
			u + stepSize * (a61 * dudt +
							   a62 * k2 +
							   a63 * k3 +
							   a64 * k4 +
							   a65 * k5);
	*/

	// a61 * dudt
	_aux = _mm256_broadcast_sd(&a61);
	_aux = _mm256_mul_pd(_aux, _dudt);

	// a62 * k2
	_aux2 = _mm256_broadcast_sd(&a62);
	_aux2 = _mm256_mul_pd(_aux2, _k2);

	// a61 * dudt + a62 * k2
	_aux = _mm256_add_pd(_aux, _aux2);

	// a63 * k3
	_aux2 = _mm256_broadcast_sd(&a63);
	_aux2 = _mm256_mul_pd(_aux2, _k3);

	// a61 * dudt + a62 * k2 + a63 * k3
	_aux = _mm256_add_pd(_aux, _aux2);

	// a64 * k4
	_aux2 = _mm256_broadcast_sd(&a64);
	_aux2 = _mm256_mul_pd(_aux2, _k4);

	// a61 * dudt + a62 * k2 + a63 * k3 + a64 * k4
	_aux = _mm256_add_pd(_aux, _aux2);

	// a65 * k5
	_aux2 = _mm256_broadcast_sd(&a65);
	_aux2 = _mm256_mul_pd(_aux2, _k5);

	// a61 * dudt + a62 * k2 + a63 * k3 + a64 * k4 + a65 * k5
	_aux = _mm256_add_pd(_aux, _aux2);

	// stepSize * (a61 * dudt + a62 * k2 + a63 * k3 + a64 * k4 + a65 * k5)
	_aux2 = _mm256_broadcast_sd(&stepSize);
	_aux = _mm256_mul_pd(_aux, _aux2);

	// u + stepSize * (a61 * dudt + a62 * k2 + a63 * k3 + a64 * k4 + a65 * k5)
	_aux = _mm256_add_pd(_aux, _u);
	_mm256_store_pd(aux, _aux);
	uTemporary.assign(aux, aux + uSize);

	dynFun(t + c6 * stepSize, uTemporary, kTemporary);
	_k6 = _mm256_maskload_pd(kTemporary.data(), _mask);

	/*
		Calculando valor na precis�o de quarta ordem
		uOutput =
			u + stepSize * (b1 * dudt +
							   b3 * k3 +
							   b4 * k4 +
							   b6 * k6);
	*/

	// b1 * dudt
	_aux = _mm256_broadcast_sd(&b1);
	_aux = _mm256_mul_pd(_aux, _dudt);

	// b3 * k3
	_aux2 = _mm256_broadcast_sd(&b3);
	_aux2 = _mm256_mul_pd(_aux2, _k3);

	// b1 * dudt + b3 * k3
	_aux = _mm256_add_pd(_aux, _aux2);

	// b4 * k4
	_aux2 = _mm256_broadcast_sd(&b4);
	_aux2 = _mm256_mul_pd(_aux2, _k4);

	// b1 * dudt + b3 * k3 + b4 * k4
	_aux = _mm256_add_pd(_aux, _aux2);

	// b6 * k6
	_aux2 = _mm256_broadcast_sd(&b6);
	_aux2 = _mm256_mul_pd(_aux2, _k6);

	// b1 * dudt + b3 * k3 + b4 * k4 + b6 * k6
	_aux = _mm256_add_pd(_aux, _aux2);

	// stepSize * (b1 * dudt + b3 * k3 + b4 * k4 + b6 * k6)
	_aux2 = _mm256_broadcast_sd(&stepSize);
	_aux = _mm256_mul_pd(_aux, _aux2);

	// u + stepSize * (b1 * dudt + b3 * k3 + b4 * k4 + b6 * k6)
	_aux = _mm256_add_pd(_aux, _u);
	_mm256_store_pd(aux, _aux);
	uOutput.assign(aux, aux + uSize);

	/*
		Estimando erro a partir da diferen�a entre quarta ordem e quinta ordem
		N�o � necess�rio calcular o valor de quinta ordem, visto que
		algebricamente � poss�vel prever qual ser� a diferen�a,
		tendo em vista que os coeficientes da tabela de Butcher s�o constantes.
	*/

	// d1 * dudt
	_aux = _mm256_broadcast_sd(&d1);
	_aux = _mm256_mul_pd(_aux, _dudt);

	// d3 * k3
	_aux2 = _mm256_broadcast_sd(&d3);
	_aux2 = _mm256_mul_pd(_aux2, _k3);

	// d1 * dudt + d3 * k3
	_aux = _mm256_add_pd(_aux, _aux2);

	// d4 * k4
	_aux2 = _mm256_broadcast_sd(&d4);
	_aux2 = _mm256_mul_pd(_aux2, _k4);

	// d1 * dudt + d3 * k3 + d4 * k4
	_aux = _mm256_add_pd(_aux, _aux2);

	// d5 * k5
	_aux2 = _mm256_broadcast_sd(&d5);
	_aux2 = _mm256_mul_pd(_aux2, _k5);

	// d1 * dudt + d3 * k3 + d4 * k4 + d5 * k5
	_aux = _mm256_add_pd(_aux, _aux2);

	// d6 * k6
	_aux2 = _mm256_broadcast_sd(&d6);
	_aux2 = _mm256_mul_pd(_aux2, _k6);

	// d1 * dudt + d3 * k3 + d4 * k4 + d5 * k5 + d6 * k6
	_aux = _mm256_add_pd(_aux, _aux2);

	// stepSize * (d1 * dudt + d3 * k3 + d4 * k4 + d5 * k5 + d6 * k6)
	_aux2 = _mm256_broadcast_sd(&stepSize);
	_aux = _mm256_mul_pd(_aux, _aux2);
	_mm256_store_pd(aux, _aux);
	uError.assign(aux, aux + uSize);

	/*
		Fim do m�todo
		End of the function
	*/
}

void CashKarp::CashKarpQualityStep(
	std::vector<double>& u,
	std::vector<double>& dudt,
	std::vector<double>& uScaled,
	double& t,
	double stepSizeTry,
	double tolerance,
	double& previousStepSize,
	double& nextStepSize,
	std::function<
	void(
		double,
		std::vector<double>&,
		std::vector<double>&)>
	& dynFun)
{
	std::size_t i;
	std::size_t uSize = u.size();
	std::vector<double> uTemporary(uSize);
	std::vector<double> uError(uSize);
	double maximumError, stepSize, temporaryStepSize, tNew;
	static double errorComparingValue = std::pow(5.0 / 0.9, 1.0 / -0.2);

	/*
		* Verificando se � possivel utilizar fun��o que faz uso dos
		intr�nsecos AVX2.
		* Ser� poss�vel se:
		-> O sistema de equa��es tiver 4 equa��es ou menos.
		-> O sistema oferecer suporte � instru��es AVX2.
	*/
	bool useAVX = (__AVX2__ == 1) && (uSize <= 4);

	/*
		Primeira tentativa ser� feita utilizando o par�metro stepSizeTry.
	*/
	stepSize = stepSizeTry;
	while (true)
	{
		if (useAVX)
			CashKarpStepAVX2(u, dudt, t, stepSize, uTemporary, uError, dynFun);
		else
			CashKarpStep(u, dudt, t, stepSize, uTemporary, uError, dynFun);

		/*
			Identificando maior erro no sistema de equa��es.
			Aqui o erro � definido como o m�dulo do erro estimado na
			solu��o de uma equa��o do sistema dividido pela toler�ncia da mesma.

			Equa��es possuem toler�ncias diferentes pois fun��es que
			apresentam valores muito maiores tendem a apresentar
			erro proporcionalmente maior tamb�m. Para contrapor tal efeito
			� utilizado o vetor de valores uScaled, que leva a ordem de
			grandeza destes valores em conta para apresentar suas toler�ncias.
		*/
		maximumError = 0.0;
		for (i = 0; i < uSize; i++) {
			double newError = std::abs(uError[i] / uScaled[i]);
			if (newError > 1.0e16) {
				newError = std::abs(uError[i] / uTemporary[i]);
			}
			maximumError = std::max<double>(maximumError, newError);
		}

		/*
			Comparando esse erro com a toler�ncia especificada.
			Se for menor, o loop � finalizado pois foi encontrada
			uma solu��o dentro da toler�ncia exigida.

			Aqui n�o � necess�rio utilizar uScaled, pois o erro
			j� foi normalizado na etapa anterior.
		*/
		maximumError /= tolerance;

		if (maximumError <= 1.0)
			break;

		/*
			Caso contr�rio, � necess�rio calcular um novo stepSize.
		*/
		temporaryStepSize = 0.9 * stepSize * std::pow(maximumError, -0.25);
		stepSize =
			(stepSize >= 0.0)
			? std::max(temporaryStepSize, 0.1 * stepSize)
			: std::min(temporaryStepSize, 0.1 * stepSize);
		/*
			Avaliando qual ser� o pr�ximo valor de t com base no
			novo stepSize. Se esse novo valor for igual ao antigo,
			alerta-se para um erro matem�tico.
		*/
		tNew = t + stepSize;

		if (tNew == t)
			throw "Mathematical error: step size is equal to zero.";
	}

	/*
		Calculando valor de stepSize para o pr�ximo passo adaptativo.
		Se o valor do erro no passo adaptativo atual for pequeno,
		tentaremos um stepSize maior, assumindo que ser� suficiente
		para o pr�ximo passo adaptativo. Caso contr�rio, ser�
		novamente diminuído.
	*/
	nextStepSize =
		(maximumError > errorComparingValue)
		? 0.9 * stepSize * std::pow(maximumError, -0.2)
		: 5.0 * stepSize;

	/*
		Armazenando valor do stepSize utilizado e
		atualizando valor da vari�vel independente t.
	*/
	previousStepSize = stepSize;
	t += previousStepSize;

	/*
		Salvando valores de u e encerrando o m�todo.
	*/
	u = uTemporary;
}

void CashKarp::CashKarpRange(
	std::vector<double>& uInitial,
	std::pair<double, double>& tSpan,
	double tolerance,
	double initialStep,
	double minimumStep,
	std::size_t maximumNumberOfSteps,
	std::function<
	void(double,
		std::vector<double>&,
		std::vector<double>&)
	>& dynFun,
	std::vector<double>& tValues,
	std::vector<
	std::vector<double>
	>& uValues)
{
	std::size_t i, numberOfSteps, uSize = uInitial.size();
	std::vector<double> uScaled(uSize);
	std::vector<double> u(uSize);
	std::vector<double> dudt(uSize);

	double t, previousStepSize, stepSize, nextStepSize;

	tValues.clear();
	uValues.clear();

	t = tSpan.first;
	tValues.push_back(t);

	u = uInitial;
	uValues.push_back(uInitial);

	stepSize =
		(tSpan.second - tSpan.first >= 0.0)
		? fabs(initialStep)
		: -fabs(initialStep);

	for (
		numberOfSteps = 0;
		numberOfSteps <= maximumNumberOfSteps;
		numberOfSteps++)
	{
		dynFun(t, u, dudt);

		for (i = 0; i < uSize; i++)
		{
			uScaled[i] =
				std::abs(u[i]) +
				std::abs(dudt[i] * stepSize) +
				1.0e-30;
		}

		double tNext = t + stepSize;
		if ((tNext - tSpan.second) * (tNext - tSpan.first) > 0.0)
			stepSize = tSpan.second - t;

		CashKarpQualityStep(
			u, dudt, uScaled, t, stepSize,
			tolerance, previousStepSize,
			nextStepSize, dynFun);

		tValues.push_back(t);
		uValues.push_back(u);

		if ((t - tSpan.second) * (tSpan.second - tSpan.first) >= 0.0)
			return; 

		stepSize = nextStepSize;
	}
}