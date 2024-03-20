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
		std::cout << std::setprecision(16);
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

RootScalarResult FixPoint(UniVarFunction &f, double x, param &parameter){
  timer stopwatch;
	stopwatch.start();

	std::string term_flag = "Success";

  double err = parameter.tol + 1.0;
  int k = 0;

  while ((err > parameter.tol) && (k < parameter.maxit)){
    double x_old = x;
    x = f(x);
    err = std::abs(x - x_old);
    k++;
  }
  if ((err > parameter.tol) && (k == parameter.maxit)){
    std::cerr << "Method Fails!" << std::endl;
		term_flag = "Fail";
		return RootScalarResult();
  }

  stopwatch.stop();
	return RootScalarResult(k, parameter.maxit, x, f(x), err, parameter.tol,
		stopwatch.get_elapsed_time(), "FIXED POINT METHOD", term_flag);
}

RootScalarResult Aitken(UniVarFunction &f, double x, param &parameter){
  timer stopwatch;
	stopwatch.start();

  std::string term_flag = "Success";

  double err = parameter.tol + 1.0;
  x = 0;
  int k = 0;

  while((err > parameter.tol) && (k < parameter.maxit)){
    double x_old = x;
    double x_onethird = f(x);
    double x_twothird = f(x_onethird);
    x = x_twothird - (std::pow(x_twothird - x_onethird, 2)/(x_twothird - (2 * x_onethird) + x));
    err = std::abs(x - x_old);
    k++;
  }
  if((err > parameter.tol) && (k == parameter.maxit)){
    std::cerr << "Method Fails!" << std::endl;
		term_flag = "Fail";
		return RootScalarResult();
  }

  stopwatch.stop();
	return RootScalarResult(k, parameter.maxit, x, f(x), err, parameter.tol,
		stopwatch.get_elapsed_time(), "AITKEN ACCELERATION METHOD", term_flag);
}
}}}
#endif
