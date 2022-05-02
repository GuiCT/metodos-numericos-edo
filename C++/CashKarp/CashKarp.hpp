/**
 * @file CashKarp.hpp
 * @author Guilherme Cesar Tomiasi (gtomiasi@gmail.com)
 * @brief Algoritmo de Cash-Karp (Runge-Kutta com passo adaptativo)
 * @date 2022-04-03
 */

 /*
  *  Essa implementação é baseada na seção 16.2 do livro "Numerical Recipes in
	 C",
	 escrito pelos autores:
	 - William H. Press
	 - Saul A. Teukolsky
	 - William T. Vetterling
	 - Brian P. Flannery
	 cujo código ISBN é 0-521-43108-5

  *  Arquivo de cabeçalho, não contém implementações.
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
	   * Retorna por referência o próximo valor de u e o erro estimado a partir do
	   * método embarcado.
	   * @param[in] u Vetor contendo atuais valores de u (entrada)
	   * @param[in] dudt Vetor contendo valores de du/dt (entrada)
	   * @param[in] t Valor de t (entrada)
	   * @param[in] stepSize Tamanho do passo (entrada)
	   * @param[out] uOutput Vetor contendo novos valores de u (saída)
	   * @param[out] uError Vetor contendo erros estimados de u (saída)
	   * @param[in] dynFun Função que calcula as derivadas de primeira ordem (entrada)
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
	  * Utiliza funções intrínsecas do conjunto de instruções AVX2.
	  * Retorna por referência o próximo valor de u e o erro estimado a partir do
	  * método embarcado.
	  * @param[in] u Vetor contendo atuais valores de u (entrada)
	  * @param[in] dudt Vetor contendo valores de du/dt (entrada)
	  * @param[in] t Valor de t (entrada)
	  * @param[in] stepSize Tamanho do passo (entrada)
	  * @param[out] uOutput Vetor contendo novos valores de u (saída)
	  * @param[out] uError Vetor contendo erros estimados de u (saída)
	  * @param[in] dynFun Função que calcula as derivadas de primeira ordem (entrada)
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
	 * Tenta realizar tal passo dentro da tolerância especificada utilizando
	 * stepSizeTry.
	 * Refina a malha até alcançar tal objetivo ou atinge o limite de refinamento
	 * da malha.
	 * @param[in, out] u Vetor contendo valores de u (entrada e saída)
	 * @param[in] dudt Vetor contendo valores de du/dt (entrada)
	 * @param[in] uScaled Vetor armazenando a tolerância para equação do sistema
	 * (entrada)
	 * @param[in, out] t Valor da variável independente t (entrada e saída)
	 * @param[in] stepSizeTry Primeiro valor assumido pelo passo (entrada)
	 * @param[in] tolerance Tolerância de erro.
	 * Dado um erro superior a esse valor, a malha será refinada (entrada)
	 * @param[in] previousStepSize Valor do passo na iteração anterior (entrada)
	 * @param[out] nextStepSize Valor do passo na próxima iteração (saída)
	 * @param[in] dynFun Função que calcula as derivadas de primeira ordem (entrada)
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
		* @brief Rotina que aplica o método de Cash-Karp para realizar a integração
		* de um determinado sistema de EDO`s em um intervalo específico.
		* @param[in] uInitial Valores iniciais do sistema de EDO`s (entrada)
		* @param[in] tSpan Intervalo de integração, com início e fim (entrada)
		* @param[in] tolerance Tolerância aceita pelo algoritmo (entrada)
		* @param[in] initialStep Passo inicial (entrada)
		* @param[in] minimumStep Passo mínimo, atualmente não implementado (entrada)
		* @param[in] maximumNumberOfSteps Quantidade máxima de iterações (entrada)
		* @param[in] dynFun Função que computa os valores do sistema de EDO`s (entrada)
		* @param[in, out] tValues Valores de t (variável independente) (entrada e saída)
		* @param[in, out] uValues Valores de u (variável dependente) (entrada e saída)
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