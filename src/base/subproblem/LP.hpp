#ifndef LP_H
#define LP_H

#include <ostream>
#include <vector>
#include <map>

/*! \struct LP
* \brief LP
*
*  Linear problem
*/
struct LP {
	LP(int number_variables, int number_constraints);
	
	int number_variables;
	int number_constraints;
	
	std::vector<double> variable_lb;
	std::vector<double> variable_ub;
	std::vector<double> constraint_lb;
	std::vector<double> constraint_ub;
	
	std::map<int,double> objective;
	std::vector<std::map<int,double> > constraints;
	
	friend std::ostream& operator<< (std::ostream &stream, LP& lp);
};

#endif // LP_H
