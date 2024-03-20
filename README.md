# Numerical-Methods
Methods Based on First-Order Approximation and Fixed-Point Methods

I. Methods Based on First-Order Approximation

Consider the function f : R → R defined by
f(x) = 0.25 cos(2x)^2 − x^2.

Implement the methods of computing a root of the equation f(x) = 0 based
on the first-order approximation of the function. In each of the methods below, take the
parameter values ωa = 0.9, ωf = 0.1, tol = 10−15 and maxit = 100 in the stopping criterion
of the corresponding algorithms.

Put the C++ implementations of each algorithm in the rootscalar.hpp file. Use the struct

cs117::root::scalar::RootScalarResult

as the return value of each function associated to the algorithm. Also, use the struct

cs117::root::scalar::param

for the input parameters. You may add additional members to this struct as needed. For the
implementation to the specific function given by (fun), use the root_scalar_firstorderex.cpp
as the filename, that is, the file where the main function is located. Print the results of each
method by invoking the member function print of RootScalarResult.

In the following list, use the initial point x_0 = 0.5 and the second point x_1 = 0 in the case
of the secant method.
1 Secant Method

![image](https://github.com/gdderije/Numerical-Methods/assets/71222985/b6464116-a725-432b-ad96-7c3516dea27d)

2 Newton Method

![image](https://github.com/gdderije/Numerical-Methods/assets/71222985/0b28b37c-5e0e-4e17-acf4-8ae7ebc4d831)

3 Steffensen Method

![image](https://github.com/gdderije/Numerical-Methods/assets/71222985/4c833bb7-4e81-40bc-8eca-befc26e857ad)

4 Inexact Newton Method (Apply the backward, forward and centered finite difference
approximations of the derivative and use the step size h = √eps.)

II. Fixed-Point Methods

A point x ∈ R is called a fixed point of g : R → R if g(x) = x. Observe that a fixed point
of g(x) = 0.5 cos(2x) is a root of the equation f(x) = 0 in Item 1. Repeat the tasks in the
previous item to the following methods. Use the filename fixpointex.cpp.

1 Fixed-Point Method

![image](https://github.com/gdderije/Numerical-Methods/assets/71222985/9bdea949-cb60-4c24-b881-7ea8a8ae4254)

2 Aitken Acceleration Method

![image](https://github.com/gdderije/Numerical-Methods/assets/71222985/d676ee52-1705-4b6f-93c4-a283517840c1)
