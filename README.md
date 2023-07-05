This directory contains the Global BFGS algorithm, the  BFGS-Wolfe algorithm and the Cautious BFGS-Armijo algorithm for solving multiobjective optimization problems described in the paper:

L. F. Prudente and D. R. Souza, Global convergence of a BFGS-type algorithm for nonconvex multiobjective optimization problems, technical report, 2023.

- MOPsolverBFGS.f90: routine containing the Global BFGS and BFGS-Wolfe algorithms
- MOPsolverCautiousBFGS.f90: routine containing the Standard Cautious BFGS-Armijo algorithm

This folder also contains the third-party free codes:

software Algencan 3.1.1

E. G. Birgin and J. M. Martı́nez, Practical augmented Lagrangian methods for constrained optimization, SIAM, 2014.
https://www.ime.usp.br/~egbirgin/tango/
Algencan is used to compute the search directions; see innersolver.f90.

subroutines dcsrch and dcstep of Moré and Thuente

J. J. Moré and D. J. Thuente, Line Search Algorithms with Guaranteed Sufficient Decrease, ACM Trans. Math. Softw., 20 (1994), pp. 286–307.
http://ftp.mcs.anl.gov/pub/MINPACK-2/csrch/
These subroutines are used as the inner solver of lsvecopt.f90 which computes a step size satisfying the (vector) Wolfe conditions.

Instructions:

File main.f90 contains the main program where you can choose the algorithm to be used. Modify myproblem.f90 routine to solve your own problem. Alternatively, set a test problem in main.f90 routine; see myproblem.f90.

The codes are written in Fortran 90. Users need to install gfortran.

In the terminal, go to the folder and type:

make

Run typing:

./MOPsolver

and see the output in the screen.

out: outer iteration number
|theta|: optimality measure
LS: flag of the line search routine to compute the step size (0 means success)
IS: flag of the inner solver routine to compute the search direction (0 means success)
#evalf: number of function evaluations
#evalg: number of gradient evaluations
