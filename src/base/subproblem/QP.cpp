#include "QP.hpp"

QP::QP(int number_variables, int number_constraints, const CSCMatrix& hessian) :
variables_bounds(number_variables),
constraints_bounds(number_constraints),
constraints(number_constraints),
hessian(hessian) {
    //this->number_variables = number_variables;
    //this->number_constraints = number_constraints;
}

std::ostream& operator<<(std::ostream &stream, QP& qp) {
    stream << "min 1/2 d^T W d + g^T d\n";
    stream << " d\n";

    /* variables */
    stream << qp.variables_bounds.size() << " variables\n";
    int i = 0;
    for (Range& range: qp.variables_bounds) {
        stream << range.lb << " <= " << "d" << i << " <= " << range.ub << "\n";
        i++;
    }

    /* objective */
    stream << "g = ";
    int number_terms = 0;
    for (std::pair<int, double> term : qp.linear_objective) {
        int index = term.first;
        double value = term.second;
        if (0 < number_terms) {
            stream << " + ";
        }
        stream << value << "*d" << index;
        number_terms++;
    }
    stream << "\n";

    /* constraints */
    stream << qp.constraints.size() << " constraints\n";
    int j = 0;
    for (Range& range: qp.constraints_bounds) {
        stream << range.lb << " <= ";
        int number_terms = 0;
        for (std::pair<int, double> term : qp.constraints[j]) {
            int index = term.first;
            double value = term.second;
            if (0 < number_terms) {
                stream << " + ";
            }
            stream << value << "*d" << index;
            number_terms++;
        }
        stream << " <= " << range.ub << "\n";
        j++;
    }

    /* Hessian */
    stream << qp.hessian;

    return stream;
}
