#ifndef HESSIANEVALUATION_H
#define HESSIANEVALUATION_H

#include <memory>
#include <vector>
#include "optimization_problem/Problem.hpp"
#include "solvers/linear/LinearSolver.hpp"
#include "solvers/linear/LinearSolverFactory.hpp"

// virtual (abstract) class
template<class SparseSymmetricMatrix>
class HessianEvaluation {
public:
   explicit HessianEvaluation(size_t dimension, size_t hessian_maximum_number_nonzeros);
   virtual ~HessianEvaluation() = default;

   SparseSymmetricMatrix hessian;
   int evaluation_count{0};

   virtual void compute(const Problem& problem, const std::vector<double>& primal_variables, double objective_multiplier,
         const std::vector<double>& constraint_multipliers) = 0;

   SparseSymmetricMatrix modify_inertia(SparseSymmetricMatrix& matrix, LinearSolver <SparseSymmetricMatrix>& linear_solver);
};

template<class SparseSymmetricMatrix>
inline HessianEvaluation<SparseSymmetricMatrix>::HessianEvaluation(size_t dimension, size_t hessian_maximum_number_nonzeros):
      hessian(dimension, hessian_maximum_number_nonzeros) {
}

template<class SparseSymmetricMatrix>
inline SparseSymmetricMatrix HessianEvaluation<SparseSymmetricMatrix>::modify_inertia(SparseSymmetricMatrix& matrix,
      LinearSolver <SparseSymmetricMatrix>& linear_solver) {
   const double beta = 1e-4;

   // Nocedal and Wright, p51
   double smallest_diagonal_entry = matrix.smallest_diagonal_entry();
   DEBUG << "The minimal diagonal entry of the Hessian is " << matrix.smallest_diagonal_entry() << "\n";

   double inertia = 0.;
   double previous_inertia = 0.;
   if (smallest_diagonal_entry <= 0.) {
      inertia = beta - smallest_diagonal_entry;
   }

   if (0. < inertia) {
      matrix.add_identity_multiple(inertia - previous_inertia);
   }

   DEBUG << "Testing factorization with inertia term " << inertia << "\n";
   linear_solver.do_symbolic_factorization(matrix);
   linear_solver.do_numerical_factorization(matrix);

   bool good_inertia = false;
   while (!good_inertia) {
      DEBUG << linear_solver.number_negative_eigenvalues() << " negative eigenvalues\n";
      if (!linear_solver.matrix_is_singular() && linear_solver.number_negative_eigenvalues() == 0) {
         good_inertia = true;
         DEBUG << "Factorization was a success with inertia " << inertia << "\n";
      }
      else {
         previous_inertia = inertia;
         inertia = (inertia == 0.) ? beta : 2 * inertia;
         matrix.add_identity_multiple(inertia - previous_inertia);
         DEBUG << "Testing factorization with inertia term " << inertia << "\n";
         linear_solver.do_numerical_factorization(matrix);
      }
   }
   return matrix;
}

// Exact Hessian

template<class SparseSymmetricMatrix>
class ExactHessianEvaluation : public HessianEvaluation<SparseSymmetricMatrix> {
public:
   explicit ExactHessianEvaluation(size_t dimension, size_t hessian_maximum_number_nonzeros) : HessianEvaluation<SparseSymmetricMatrix>(dimension,
         hessian_maximum_number_nonzeros) {
   }

   ~ExactHessianEvaluation() override = default;

   void compute(const Problem& problem, const std::vector<double>& primal_variables, double objective_multiplier,
         const std::vector<double>& constraint_multipliers) override {
      /* compute Hessian */
      problem.lagrangian_hessian(primal_variables, objective_multiplier, constraint_multipliers, this->hessian);
      this->evaluation_count++;
   }
};

//template <class SparseSymmetricMatrix>
//class ConvexifiedExactHessianEvaluation : public HessianEvaluation<SparseSymmetricMatrix> {
//public:
//   ConvexifiedExactHessianEvaluation(size_t dimension, size_t hessian_maximum_number_nonzeros, const std::string& linear_solver_name):
//      HessianEvaluation<SparseSymmetricMatrix>(dimension, hessian_maximum_number_nonzeros), linear_solver(LinearSolverFactory<MA57Solver>::create
//      (linear_solver_name)) {
//   }
//   ~ConvexifiedExactHessianEvaluation() override = default;
//
//   void compute(const Problem& problem, const std::vector<double>& primal_variables, double objective_multiplier,
//         const std::vector<double>& constraint_multipliers) override {
//      /* compute Hessian */
//      problem.lagrangian_hessian(primal_variables, objective_multiplier, constraint_multipliers, this->hessian);
//      this->evaluation_count++;
//      /* modify the inertia to make the problem strictly convex */
//      DEBUG << "hessian before convexification: " << this->hessian;
//      this->hessian = this->modify_inertia(this->hessian, *this->linear_solver);
//   }
//
//protected:
//   std::unique_ptr<LinearSolver<COOSymmetricMatrix> > linear_solver; /*!< Solver that computes the inertia */
//};

template<class SparseSymmetricMatrix>
class HessianEvaluationFactory {
public:
   static std::unique_ptr <HessianEvaluation<SparseSymmetricMatrix>> create(const std::string& hessian_evaluation_method, size_t dimension,
         size_t hessian_maximum_number_nonzeros, bool convexify) {
      if (hessian_evaluation_method == "exact") {
         if (convexify) {
            //   return std::make_unique<ConvexifiedExactHessianEvaluation<SparseSymmetricMatrix> >(dimension, hessian_maximum_number_nonzeros, "MA57");
         }
         //else {
         return std::make_unique<ExactHessianEvaluation<SparseSymmetricMatrix> >(dimension, hessian_maximum_number_nonzeros);
         //}
      }
      else {
         throw std::invalid_argument("Hessian evaluation method " + hessian_evaluation_method + " does not exist");
      }
   }
};

#endif // HESSIANEVALUATION_H
