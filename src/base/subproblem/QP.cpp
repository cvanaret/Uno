#include "QP.hpp"

QP::QP(int number_variables, int number_constraints, const CSCMatrix& hessian) :
variable_lb(number_variables), variable_ub(number_variables),
constraint_lb(number_constraints), constraint_ub(number_constraints),
constraints(number_constraints),
hessian(hessian) {
    this->number_variables = number_variables;
    this->number_constraints = number_constraints;
}

std::ostream& operator<<(std::ostream &stream, QP& qp) {
    stream << "min 1/2 d^T W d + g^T d\n";
    stream << " d\n";

    /* variables */
    stream << qp.number_variables << " variables\n";
    for (int i = 0; i < qp.number_variables; i++) {
        stream << qp.variable_lb[i] << " <= " << "d" << i << " <= " << qp.variable_ub[i] << "\n";
    }

    /* objective */
    stream << "g = ";
    int number_terms = 0;
    for (std::pair<int, double> term : qp.objective) {
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
    stream << qp.number_constraints << " constraints\n";
    for (int j = 0; j < qp.number_constraints; j++) {
        stream << qp.constraint_lb[j] << " <= ";
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
        stream << " <= " << qp.constraint_ub[j] << "\n";
    }

    /* Hessian */
    stream << qp.hessian;

    return stream;
}
