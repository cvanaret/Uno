#ifndef QPSOLVER_H
#define QPSOLVER_H

#include <vector>
#include <type_traits>
#include "LPSolver.hpp"

/*! \class QPSolver
 * \brief QP solver
 *
 */
template<class SparseSymmetricMatrix>
class QPSolver : public LPSolver {
public:
   // make this type accessible in templates
   using matrix_type = SparseSymmetricMatrix;

   QPSolver() = default;
   ~QPSolver() override = default;
   virtual Direction solve_QP(const std::vector <Range>& variables_bounds, const std::vector <Range>& constraints_bounds,
         const SparseVector2<double>& linear_objective, const std::vector <SparseVector>& constraints_jacobian, const SparseSymmetricMatrix& hessian,
         const std::vector<double>& initial_point) = 0;
   virtual Direction solve_LP(const std::vector <Range>& variables_bounds, const std::vector <Range>& constraints_bounds, const SparseVector2<double>&
         linear_objective, const std::vector <SparseVector>& constraints_jacobian, const std::vector<double>& initial_point) = 0;
};

#endif // QPSOLVER_H
