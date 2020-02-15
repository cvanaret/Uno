#ifndef QP_H
#define QP_H

#include <ostream>
#include <vector>
#include <map>
#include "Matrix.hpp"
#include "Constraint.hpp"

/*! \struct QP
* \brief QP
*
*  Quadratic problem
*/
struct QP {
	/*!
	 *  Constructor
	 * 
	 * \param hessian_column_start: column description of sparse Hessian
	 * \param hessian_row_number: row description of sparse Hessian
	 */
	QP(int number_variables, int number_constraints, const CSCMatrix& hessian);

	//int number_variables;
	//int number_constraints;
	
    std::vector<Range> variables_bounds;
    std::vector<Range> constraints_bounds;
	
	std::map<int,double> linear_objective;
	std::vector<std::map<int,double> > constraints;
	const CSCMatrix& hessian; /*!< Sparse Hessian */
	
	friend std::ostream& operator<< (std::ostream &stream, QP& qp);
};

#endif // QP_H
