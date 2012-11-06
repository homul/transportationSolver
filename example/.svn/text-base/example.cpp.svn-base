#include <iostream>
#include "primalSolver.h"
#include "testcommon.h"

using namespace std;
using namespace DD;

int main()
{
    size_t a_size=3,b_size=4;
	// The solver can be initialized during construction like here:
    Matrix bin(a_size, b_size, 0.0);
	PrimalSolver<false, Matrix, SparseArray> solver(a_size,b_size,bin);//<~ minimization problem

	//creating 2 empty sparse arrays
	SparseArray a(a_size, 0.0);
	SparseArray b(b_size, 0.0);

	//assigning values
	a(0) = 0.7;
	a(1) = 0.3;
	b(0) = 1.0;

	std::pair<bool, floatType> res = solver.Solve(a, b);

	if (res.first)
	{
	 cout << "Transportation problem was successfully solved! Optimal value="<<res.second<< endl;
	}
	else
	 cout << "Transportation problem was not solved. Please contact the program author(s)" <<endl;

	//============================================================================
	//alternatively an empty constructor  can be used. The object can be initialized afterwards
	PrimalSolver<true, Matrix, SparseArray> solver1();//maximization problem
	solver.Init(a_size,b_size,bin);

	bin(0,0)=10; //solver keeps only a pointer to the matrix bin, hence keeping it in memory is in your responsibility
	res = solver.Solve(a, b);

		if (res.first)
		{
		 cout << "Transportation problem was successfully solved! Optimal value="<<res.second<< endl;
		 Matrix solution(a_size,b_size);//allocating memory for a solution matrix
		 floatType optimalValue=solver.GetSolution(&solution);
		 cout << "Solution = "<<solution << endl;
		 cout << "Its quality="<<optimalValue <<endl;//optimalValue =  res.second
		}
		else
		 cout << "Transportation problem was not solved. Please contact the program author(s)" <<endl;

	//============================================================================
	return 0;
}
