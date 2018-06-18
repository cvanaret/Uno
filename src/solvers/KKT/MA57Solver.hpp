#ifndef MA57SOLVER_H
#define MA57SOLVER_H

#include <vector>
#include "Matrix.hpp"
#include "LocalSolution.hpp"

/*! \class MA57Solver
* \brief Interface for MA57
*
*  Interface to the sparse symmetric linear solver MA57
*/
class MA57Solver {
	public:
		MA57Solver();
		
		int use_fortran;
		LocalSolution solve(COOMatrix& matrix, std::vector<double> rhs);
		
	private:
		/* for ma57id_ */
		std::vector<double> cntl_;
		std::vector<int> icntl_;
		std::vector<int> info_;
		std::vector<double> rinfo_;
};

#endif // MA57SOLVER_H
