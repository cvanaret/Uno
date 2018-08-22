#ifndef MA57SOLVER_H
#define MA57SOLVER_H

#include <vector>
#include "Matrix.hpp"
#include "LocalSolution.hpp"

struct MA57Data {
	std::vector<double> fact;
	int lfact;
	std::vector<int> ifact;
	int lifact;
	std::vector<int> iwork;	
};

/*! \class MA57Solver
* \brief Interface for MA57
* see https://github.com/YimingYAN/linSolve
*
*  Interface to the sparse symmetric linear solver MA57
*/
class MA57Solver {
	public:
		MA57Solver();
		
		int use_fortran;
		
		MA57Data factorize(COOMatrix& matrix);
		std::vector<double> solve(COOMatrix& matrix, std::vector<double>& rhs, MA57Data& data);
		int number_negative_eigenvalues();
                bool matrix_is_singular();
		
	private:
		/* for ma57id_ */
		std::vector<double> cntl_;
		std::vector<int> icntl_;
		std::vector<int> info_;
		std::vector<double> rinfo_;
};

#endif // MA57SOLVER_H
