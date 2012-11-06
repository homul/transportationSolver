transportationSolver
====================

Transportation solver implementation to obtain a primal solution of the LP relaxation of the MRF MAP-inference problem
-----------------------------------------------------------------------------------
The software implements a simplex solver for transportation problems: http://en.wikibooks.org/wiki/Operations_Research/Transportation_and_Assignment_Problem
used to compute a primal solution of the LP-relaxed MAP-inference problem for graphical models. Versions of this code were used to produce resuls for papers

[1] B. Savchynskyy, J. H. Kappes, S. Schmidt, C. Schnörr
A Study of Nesterov's Scheme for Lagrangian Decomposition and MAP Labeling
CVPR 2011 

[2] S. Schmidt, B. Savchynskyy, J. H. Kappes, C. Schnörr
Evaluation of a First-Order Primal-Dual Algorithm for MRF Energy Minimization
In EMMCVPR, 2011. Springer, pp.89-103. 

[3] B. Savchynskyy, S. Schmidt, J. H. Kappes, C. Schnörr
Efficient MRF Energy Minimization via Adaptive Diminishing Smoothing
In UAI, 2012. Note: In press. 

[4] B. Savchynskyy, S. Schmidt
Getting Feasible Variable Estimates From Infeasible Ones: MRF Local Polytope Study
arXiv:1210.4081 Submitted Oct. 2012 [Bib] [PDF]

available at http://hci.iwr.uni-heidelberg.de/Staff/bsavchyn/

Please site [4] or [1] if you use it. The code was tested under Ubuntu 11 operating system with gcc v.4.4 compiler.

The solver itself is in the /src folder. It contains only headers and hence does not require any precompilation and is dependent on the standard library only. Look into /example/example.cpp to see how it can be used. To compile /example.cpp BOOST::UBLAS library is required. Further details - in src/primalSolver.h .

The folder /tests contains unit and regression tests. Use ccmake and make commants from the package root to build them. You will need a CPPUNIT and BOOST::UBLAS libraries (both are included into Ubuntu 10-12 distributions) to do build the tests. Run ./tests from /tests directory if you want to check tests. You should have right permissions for /tests to do this, since tests will create .tst files, which content will be compared to the content of existing .chk files.
---------------------------------------------------
Versions history:
1.0 - initial version
1.1 - FIXED: Numerical bug when dealing with degenerate problems 
