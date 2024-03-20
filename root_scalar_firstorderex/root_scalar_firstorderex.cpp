/* CMSC 117 FINAL REQUIREMENT
* Name: Gisselle Derije and Kianne Lee
* Date: 12 December 2020
*/

#include "rootscalar.hpp"

double f(double x)						{return 0.25*std::pow(std::cos(2.*x), 2) - std::pow(x, 2);}
double f_derivative(double x)	{return -0.5*std::sin(4.*x) - 2*x;}
double forward(double x)			{return (f(x + std::sqrt(std::numeric_limits<double>::epsilon())) - f(x)) / std::sqrt(std::numeric_limits<double>::epsilon());}
double backward(double x)			{return (f(x) - (f(x - std::sqrt(std::numeric_limits<double>::epsilon()))))/std::sqrt(std::numeric_limits<double>::epsilon());}
double center(double x)				{return (f(x + std::sqrt(std::numeric_limits<double>::epsilon())) - f(x - std::sqrt(std::numeric_limits<double>::epsilon()))) / (2.*std::sqrt(std::numeric_limits<double>::epsilon()));}

int main()
{
	double x = 0.5;

	// print filename
	std::cout << "File: " << __FILE__ << std::endl;

	// parameter object
	cs117::root::scalar::param parameter;
	parameter.tol = 10e-15;
	parameter.maxit = 100;
	parameter.wa = 0.9;
	parameter.wf = 0.1;


	// approximate roots
	cs117::root::scalar::RootScalarResult
	result = Secant(f, x, 0, parameter);
	result.print();
	result = Newton(f, f_derivative, 0.5, parameter);
	result.print();
	result = Steffensen(f, 0.5, parameter);
	result.print();
	result = Center(f, center, x, parameter);
	result.print();
	result = Forward(f, forward, x, parameter);
	result.print();
	result = Backward(f, backward, x, parameter);
	result.print();
	return 0;
}
