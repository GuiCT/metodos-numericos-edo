/**
 * @file CashKarp.hpp
 * @author Guilherme Cesar Tomiasi (gtomiasi@gmail.com)
 * @brief Algoritmo de Cash-Karp (Runge-Kutta com passo adaptativo)
 * @date 2022-04-03
 */

 /*
  *  Essa implementa��o � baseada na se��o 16.2 do livro "Numerical Recipes in
	 C",
	 escrito pelos autores:
	 - William H. Press
	 - Saul A. Teukolsky
	 - William T. Vetterling
	 - Brian P. Flannery
	 cujo c�digo ISBN � 0-521-43108-5

  *  Arquivo de cabe�alho, n�o cont�m implementa��es.
 */

#include <vector>
#include <functional>

#ifdef __AVX2__
#define CASH_KARP_STEP CashKarpStepAVX2
#else
#define CASH_KARP_STEP CashKarpStep
#endif

namespace CashKarp {
	/**
	   * @brief Rotina utilizada para calcular um passo utilizando o Runge-Kutta de
	   * Cash-Karp.
	   * Retorna por refer�ncia o pr�ximo valor de u e o erro estimado a partir do
	   * m�todo embarcado.
	   * @param[in] u Vetor contendo atuais valores de u (entrada)
	   * @param[in] dudt Vetor contendo valores de du/dt (entrada)
	   * @param[in] t Valor de t (entrada)
	   * @param[in] stepSize Tamanho do passo (entrada)
	   * @param[out] uOutput Vetor contendo novos valores de u (sa�da)
	   * @param[out] uError Vetor contendo erros estimados de u (sa�da)
	   * @param[in] dynFun Fun��o que calcula as derivadas de primeira ordem (entrada)
	   */
	void CashKarpStep(
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
		>& dynFun);

	/**
	  * @brief Rotina utilizada para calcular um passo utilizando o Runge-Kutta de
	  * Cash-Karp.
	  * Utiliza fun��es intr�nsecas do conjunto de instru��es AVX2.
	  * Retorna por refer�ncia o pr�ximo valor de u e o erro estimado a partir do
	  * m�todo embarcado.
	  * @param[in] u Vetor contendo atuais valores de u (entrada)
	  * @param[in] dudt Vetor contendo valores de du/dt (entrada)
	  * @param[in] t Valor de t (entrada)
	  * @param[in] stepSize Tamanho do passo (entrada)
	  * @param[out] uOutput Vetor contendo novos valores de u (sa�da)
	  * @param[out] uError Vetor contendo erros estimados de u (sa�da)
	  * @param[in] dynFun Fun��o que calcula as derivadas de primeira ordem (entrada)
	  */
	void CashKarpStepAVX2(
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
		>& dynFun);

	/**
	 * @brief Rotina utilizada para calcular um passo adaptativo via Runge-Kutta de
	 * Cash-Karp.
	 * Tenta realizar tal passo dentro da toler�ncia especificada utilizando
	 * stepSizeTry.
	 * Refina a malha at� alcan�ar tal objetivo ou atinge o limite de refinamento
	 * da malha.
	 * @param[in, out] u Vetor contendo valores de u (entrada e sa�da)
	 * @param[in] dudt Vetor contendo valores de du/dt (entrada)
	 * @param[in] uScaled Vetor armazenando a toler�ncia para equa��o do sistema
	 * (entrada)
	 * @param[in, out] t Valor da vari�vel independente t (entrada e sa�da)
	 * @param[in] stepSizeTry Primeiro valor assumido pelo passo (entrada)
	 * @param[in] tolerance Toler�ncia de erro.
	 * Dado um erro superior a esse valor, a malha ser� refinada (entrada)
	 * @param[in] previousStepSize Valor do passo na itera��o anterior (entrada)
	 * @param[out] nextStepSize Valor do passo na pr�xima itera��o (sa�da)
	 * @param[in] dynFun Fun��o que calcula as derivadas de primeira ordem (entrada)
	 */
	void CashKarpQualityStep(
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
		& dynFun);

	/**
		* @brief Rotina que aplica o m�todo de Cash-Karp para realizar a integra��o
		* de um determinado sistema de EDO`s em um intervalo espec�fico.
		* @param[in] uInitial Valores iniciais do sistema de EDO`s (entrada)
		* @param[in] tSpan Intervalo de integra��o, com in�cio e fim (entrada)
		* @param[in] tolerance Toler�ncia aceita pelo algoritmo (entrada)
		* @param[in] initialStep Passo inicial (entrada)
		* @param[in] minimumStep Passo m�nimo, atualmente n�o implementado (entrada)
		* @param[in] maximumNumberOfSteps Quantidade m�xima de itera��es (entrada)
		* @param[in] dynFun Fun��o que computa os valores do sistema de EDO`s (entrada)
		* @param[in, out] tValues Valores de t (vari�vel independente) (entrada e sa�da)
		* @param[in, out] uValues Valores de u (vari�vel dependente) (entrada e sa�da)
	*/
	void CashKarpRange(
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
		>& uValues);
}