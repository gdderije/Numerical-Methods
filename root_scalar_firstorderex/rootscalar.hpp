/* CMSC 117 FINAL REQUIREMENT
* Name: Gisselle Derije and Kianne Lee
* Date: 12 December 2020
*/

/* C++ header file for root-finding algorithms. */
#ifndef ROOTSCALAR_HPP_INCLUDE
#define ROOTSCALAR_HPP_INCLUDE

// standard library includes
#include <iostream>
#include <iomanip>
#include <limits>
#include <cmath>

// local include
#include "timer.hpp"

namespace cs117 { namespace root { namespace scalar
{
// type for a funciton with input and output having the type double
using UniVarFunction = double(double);

// struct for parameters in scalar-root finding algorithms
struct param
{
	double 		tol;
	int 			maxit;
	double 		wa;
	double 		wf;
	double 		h;
};

// struct for solution to scalar-root finding problems
struct RootScalarResult
{
	int 				numit;
	int 				maxit;
	double 			x;
	double 			funval;
	double 			error;
	double			tol;
	double 			elapsed_time;
	std::string method_name;
	std::string	termination_flag;

	// default constructor
	RootScalarResult(){}

	// user-defined constructor
	RootScalarResult(const int &numit, const int &maxit, const double &x,
		const double &funval, const double &error, const double &tol,
		const double &elapsed_time, const std::string &method_name,
		const std::string &termination_flag)
	{
		this->numit 							= numit;
		this->maxit 							= maxit;
		this->x 									= x;
		this->funval 							= funval;
		this->error 							= error;
		this->tol 								= tol;
		this->elapsed_time 				= elapsed_time;
		this->method_name					= method_name;
		this->termination_flag  	= termination_flag;
	}

	void print()
	{
		std::cout << "ROOT FINDER:                     " << method_name << std::endl;
		std::cout << std::setprecision(17);
		std::cout << std::fixed;
		std::cout << "APPROXIMATE ROOT / LAST ITERATE: " << x << std::endl;
		std::cout << "TERMINATION:                     " << termination_flag
			<< std::endl;
		std::cout << std::scientific;
		std::cout << "FUNCTION VALUE:                  " << funval << std::endl;
		std::cout << "ERROR:                           " << error << std::endl;
		std::cout << "TOLERANCE:                       " << tol << std::endl;
		std::cout << "NUM ITERATIONS:                  " << numit << std::endl;
		std::cout << "MAX ITERATIONS:                  " << maxit << std::endl;
		std::cout << "ELAPSED TIME:                    " << elapsed_time
							<< " seconds												 " << std::endl << std::endl;
		std::cout << std::defaultfloat;
	}
};

RootScalarResult Secant(UniVarFunction &f, double x0, double x1, param &parameter){
	timer stopwatch;
	stopwatch.start();

	std::string term_flag = "Success";

	double err = parameter.tol + 1.0;
	double f0 = f(x0);
	double f1 = f(x1);
	int k = 1;

	while((err > parameter.tol) && (k < parameter.maxit)){
		double q = (f1 - f0)/(x1 - x0);
		double xtemp = x1;
		x1 -= f1/q;
		x0 = xtemp;
		f0 = f1;
		f1 = f(x1);
		err = (parameter.wa*std::abs(x1-x0)) + (parameter.wf*std::abs(f1));
		k++;
	}
	if((err > parameter.tol) && (k == parameter.maxit)){
		std::cerr << "Method Fails!" << std::endl;
		term_flag = "Fail";
		return RootScalarResult();
	}
	stopwatch.stop();
	return RootScalarResult(k, parameter.maxit, x1, f1, err, parameter.tol,
		stopwatch.get_elapsed_time(), "SECANT METHOD", term_flag);
}

RootScalarResult Steffensen(UniVarFunction &f, double x, param &parameter){
	timer stopwatch;
	stopwatch.start();

	std::string term_flag = "Success";

	double err = parameter.tol + 1.0;
	double fx = f(x);
	int k = 0;

	while((err > parameter.tol) && (k < parameter.maxit)){
		double x_old = x;
		double q = (f(x + fx) - fx)/fx;
		x -= fx/q;
		fx = f(x);
		err = parameter.wa*std::abs(x - x_old) + parameter.wf*std::abs(fx);
		k++;
	}
	if((err > parameter.tol) && (k == parameter.maxit)){
		std::cerr << "Method Fails!" << std::endl;
		term_flag = "Fail";
		return RootScalarResult();
	}
	stopwatch.stop();
	return RootScalarResult(k, parameter.maxit, x, fx, err, parameter.tol,
		stopwatch.get_elapsed_time(), "STEFFENSEN METHOD", term_flag);
}

RootScalarResult Newton(UniVarFunction &f, UniVarFunction &f_derivative, double x, param &parameter){
	timer stopwatch;
	stopwatch.start();

	std::string term_flag = "Success";

	double err = parameter.tol + 1.0;
	double fx = f(x);
	int k = 0;

	while((err > parameter.tol) && (k < parameter.maxit)){
		double x_old = x;
		x -= fx/f_derivative(x);
		fx = f(x);
		err = parameter.wa*std::abs(x - x_old) + parameter.wf*std::abs(fx);
		k++;
	}
	if((err > parameter.tol) && (k == parameter.maxit)){
		std::cerr << "Method Fails!" << std::endl;
		term_flag = "Fail";
		return RootScalarResult();
	}
	stopwatch.stop();
	return RootScalarResult(k, parameter.maxit, x, fx, err, parameter.tol,
		stopwatch.get_elapsed_time(), "NEWTON METHOD", term_flag);
}

RootScalarResult Center(UniVarFunction &f, UniVarFunction &center, double x, param &parameter){
	timer stopwatch;
	stopwatch.start();

	std::string term_flag = "Success";

	double err = parameter.tol + 1.0;
	double fx = f(x);
	int k = 0;

	while((err > parameter.tol) && (k < parameter.maxit)){
		double x_old = x;
		x -= fx/center(x);
		fx = f(x);
		err = parameter.wa*std::abs(x - x_old) + parameter.wf*std::abs(fx);
		k++;
	}
	if((err > parameter.tol) && (k == parameter.maxit)){
		std::cerr << "Method Fails!" << std::endl;
		term_flag = "Fail";
		return RootScalarResult();
	}
	stopwatch.stop();
	return RootScalarResult(k, parameter.maxit, x, fx, err, parameter.tol,
		stopwatch.get_elapsed_time(), "CENTER", term_flag);
}

RootScalarResult Forward(UniVarFunction &f, UniVarFunction &forward, double x, param &parameter){
	timer stopwatch;
	stopwatch.start();

	std::string term_flag = "Success";

	double err = parameter.tol + 1.0;
	double fx = f(x);
	int k = 0;

	while((err > parameter.tol) && (k < parameter.maxit)){
		double x_old = x;
		x -= fx/forward(x);
		fx = f(x);
		err = parameter.wa*std::abs(x - x_old) + parameter.wf*std::abs(fx);
		k++;
	}
	if((err > parameter.tol) && (k == parameter.maxit)){
		std::cerr << "Method Fails!" << std::endl;
		term_flag = "Fail";
		return RootScalarResult();
	}
	stopwatch.stop();
	return RootScalarResult(k, parameter.maxit, x, fx, err, parameter.tol,
		stopwatch.get_elapsed_time(), "FORWARD", term_flag);
}

RootScalarResult Backward(UniVarFunction &f, UniVarFunction &backward, double x, param &parameter){
	timer stopwatch;
	stopwatch.start();

	std::string term_flag = "Success";

	double err = parameter.tol + 1.0;
	double fx = f(x);
	int k = 0;

	while((err > parameter.tol) && (k < parameter.maxit)){
		double x_old = x;
		x -= fx/backward(x);
		fx = f(x);
		err = parameter.wa*std::abs(x - x_old) + parameter.wf*std::abs(fx);
		k++;
	}
	if((err > parameter.tol) && (k == parameter.maxit)){
		std::cerr << "Method Fails!" << std::endl;
		term_flag = "Fail";
		return RootScalarResult();
	}
	stopwatch.stop();
	return RootScalarResult(k, parameter.maxit, x, fx, err, parameter.tol,
		stopwatch.get_elapsed_time(), "BACKWARD", term_flag);
}
}}}	// end of namespace cs117::root::scalar


#endif
