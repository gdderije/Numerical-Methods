/* CMSC 117 FINAL REQUIREMENT
 * Name: Gisselle Derije and Kianne Lee
 * Date: 12 December 2020
 */

#include "rootscalar.hpp"

double f(double x) {
	return 0.5 * std::cos(2 * x);
}

int main() {
	double x = f(x);

	// print filename
	std::cout << "File: " << __FILE__ << std::endl;

	// parameter object
	cs117::root::scalar::param parameter;
	parameter.tol = 10e-16;
	parameter.maxit = 100;

	// approximate roots
	cs117::root::scalar::RootScalarResult result = FixPoint(f, x, parameter);
	result.print();
	result = Aitken(f, x, parameter);
	result.print();
	return 0;
}