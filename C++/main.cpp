#include "CashKarp.hpp"
#include <vector>
#include <iostream>

int main(void)
{
	std::vector<double> uInitial = { 0.0, 0.0, 0.33206 };
	std::pair<double, double> tSpan = { 0.0, 50000.0 };
	double tolerance = 1e-5;
	double initialStep = 1e-1;
	double minimumStep = 1e-10;
	std::size_t maximumNumberOfSteps = 100000;
	std::function<
		void(
			double,
			std::vector<double>&,
			std::vector<double>&)>
		dynFun = [](
			double t,
			std::vector<double>& u,
			std::vector<double>& dudt)
	{
		dudt[0] = u[1];
		dudt[1] = u[2];
		dudt[2] = (-1.0 / 2.0) * u[0] * u[2];
	};
	std::vector<double> tValues;
	tValues.reserve(1000);
	std::vector<
		std::vector<double>>
		uValues;
	uValues.reserve(1000);

	CashKarp::CashKarpRange(
		uInitial,
		tSpan,
		tolerance,
		initialStep,
		minimumStep,
		maximumNumberOfSteps,
		dynFun,
		tValues,
		uValues);
	std::cout << "u'[" << tValues[tValues.size() - 1] << "]: " << uValues[uValues.size() - 1][1];

	exit(EXIT_SUCCESS);
}