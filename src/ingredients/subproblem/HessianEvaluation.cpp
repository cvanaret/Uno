#include <exception>
#include <cassert>
#include "HessianEvaluation.hpp"
#include "LinearSolverFactory.hpp"
#include "Vector.hpp"

//template <class MatrixType>
//HessianEvaluation<MatrixType>::HessianEvaluation(size_t dimension, size_t hessian_maximum_number_nonzeros) : dimension(dimension),
//      hessian(dimension, hessian_maximum_number_nonzeros, 1) {
//}

/* Exact Hessian */
//template <class MatrixType>
//ExactHessianEvaluation<MatrixType>::ExactHessianEvaluation(size_t dimension, size_t hessian_maximum_number_nonzeros) : HessianEvaluation<MatrixType>(dimension,
//      hessian_maximum_number_nonzeros) {
//}

//template <class MatrixType>
//void ExactHessianEvaluation<MatrixType>::compute(const Problem& problem, const std::vector<double>& primal_variables, double objective_multiplier,
//      const std::vector<double>& constraint_multipliers) {
//   /* compute Hessian */
//   problem.lagrangian_hessian(primal_variables, objective_multiplier, constraint_multipliers, this->hessian);
//   this->evaluation_count++;
//   //this->objective_multiplier_ = objective_multiplier;
//}

/* Exact Hessian with inertia control */

//template <class MatrixType>
//ExactHessianInertiaControlEvaluation<MatrixType>::ExactHessianInertiaControlEvaluation(size_t dimension, size_t hessian_maximum_number_nonzeros,
//      const std::string& linear_solver_name) : HessianEvaluation<MatrixType>(dimension, hessian_maximum_number_nonzeros),
//      linear_solver_(LinearSolverFactory::create(linear_solver_name)) {
//}

//template <class MatrixType>
//void ExactHessianInertiaControlEvaluation<MatrixType>::compute(const Problem& problem, const std::vector<double>& primal_variables,
//      double objective_multiplier, const std::vector<double>& constraint_multipliers) {
//   /* compute Hessian */
//   problem.lagrangian_hessian(primal_variables, objective_multiplier, constraint_multipliers, this->hessian);
//   this->evaluation_count++;
//   DEBUG << "hessian before convexification: " << this->hessian;
//   //assert(false && "ExactHessianInertiaControlEvaluation: inertia correction not implemented");
//   /* modify the inertia to make the problem strictly convex */
//   this->hessian = this->modify_inertia(this->hessian, *this->linear_solver_);
//}

//template <class MatrixType>
//CSCSymmetricMatrix HessianEvaluation<MatrixType>::modify_inertia(CSCSymmetricMatrix& matrix, LinearSolver& linear_solver) {
//   double beta = 1e-4;
//
//   // Nocedal and Wright, p51
//   double smallest_diagonal_entry = matrix.smallest_diagonal_entry();
//   DEBUG << "The minimal diagonal entry of the Hessian is " << matrix.smallest_diagonal_entry() << "\n";
//
//   double inertia = 0.;
//   double previous_inertia = 0.;
//   if (smallest_diagonal_entry <= 0.) {
//      inertia = beta - smallest_diagonal_entry;
//   }
//
//   if (0. < inertia) {
//      matrix = matrix.add_identity_multiple(inertia - previous_inertia);
//   }
//   COOSymmetricMatrix coo_hessian = matrix.to_COO();
//
//   DEBUG << "Testing factorization with inertia term " << inertia << "\n";
//   linear_solver.do_symbolic_factorization(coo_hessian);
//   linear_solver.do_numerical_factorization(coo_hessian);
//
//   bool good_inertia = false;
//   while (!good_inertia) {
//      DEBUG << linear_solver.number_negative_eigenvalues() << " negative eigenvalues\n";
//      if (!linear_solver.matrix_is_singular() && linear_solver.number_negative_eigenvalues() == 0) {
//         good_inertia = true;
//         DEBUG << "Factorization was a success with inertia " << inertia << "\n";
//      }
//      else {
//         previous_inertia = inertia;
//         inertia = (inertia == 0.) ? beta : 2*inertia;
//         matrix = matrix.add_identity_multiple(inertia - previous_inertia);
//         coo_hessian = matrix.to_COO();
//         DEBUG << "Testing factorization with inertia term " << inertia << "\n";
//         linear_solver.do_numerical_factorization(coo_hessian);
//      }
//   }
//   return matrix;
//}

/* BFGS Hessian */

//template <class MatrixType>
//BFGSHessianEvaluation<MatrixType>::BFGSHessianEvaluation(size_t dimension, size_t hessian_maximum_number_nonzeros) :
//   HessianEvaluation<MatrixType>(dimension, hessian_maximum_number_nonzeros),
//   previous_hessian_(dimension, hessian_maximum_number_nonzeros, 1), previous_x_(dimension) {
//}

//template <class MatrixType>
//void BFGSHessianEvaluation<MatrixType>::compute(const Problem& problem, const std::vector<double>& primal_variables, double objective_multiplier,
//      const std::vector<double>& constraint_multipliers) {
//   // the BFGS Hessian is already positive definite, do not convexify
//   problem.lagrangian_hessian(primal_variables, objective_multiplier, constraint_multipliers, this->hessian);
//   this->evaluation_count++;
//}

/* Factory */

//template <class MatrixType>
//std::unique_ptr<HessianEvaluation<MatrixType> > HessianEvaluationFactory<MatrixType>::create(const std::string& hessian_evaluation_method, size_t
//dimension, size_t hessian_maximum_number_nonzeros, bool convexify) {
//   if (hessian_evaluation_method == "exact") {
//      if (convexify) {
//         return std::make_unique<ExactHessianInertiaControlEvaluation<MatrixType> >(dimension, hessian_maximum_number_nonzeros, "MA57");
//      }
//      else {
//         return std::make_unique<ExactHessianEvaluation<MatrixType> >(dimension, hessian_maximum_number_nonzeros);
//      }
//   }
//      //    else if (hessian_evaluation_method == "BFGS") {
//      //        return std::make_unique<BFGSHessianEvaluation>(dimension, hessian_maximum_number_nonzeros);
//      //    }
//   else {
//      throw std::invalid_argument("Hessian evaluation method " + hessian_evaluation_method + " does not exist");
//   }
//}
