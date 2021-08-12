#ifndef HESSIANEVALUATION_H
#define HESSIANEVALUATION_H

#include <ostream>
#include <map>
#include <memory>
#include <vector>
#include "Problem.hpp"
#include "Iterate.hpp"
#include "LinearSolver.hpp"
#include "LinearSolverFactory.hpp"

// virtual (abstract) class
template <class MatrixType>
class HessianEvaluation {
public:
   explicit HessianEvaluation(size_t dimension, size_t hessian_maximum_number_nonzeros): dimension(dimension), hessian(dimension,
         hessian_maximum_number_nonzeros) {
   }

   virtual ~HessianEvaluation() = default;

   size_t dimension;
   MatrixType hessian;
   int evaluation_count{0};

   virtual void compute(const Problem& problem, const std::vector<double>& primal_variables, double objective_multiplier,
         const std::vector<double>& constraint_multipliers) = 0;

   MatrixType modify_inertia(MatrixType& matrix, LinearSolver<MatrixType>& linear_solver) {
      double beta = 1e-4;

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
            inertia = (inertia == 0.) ? beta : 2*inertia;
            matrix.add_identity_multiple(inertia - previous_inertia);
            DEBUG << "Testing factorization with inertia term " << inertia << "\n";
            linear_solver.do_numerical_factorization(matrix);
         }
      }
      return matrix;
   }
};

template <class MatrixType>
class ExactHessianEvaluation : public HessianEvaluation<MatrixType> {
public:
   explicit ExactHessianEvaluation(size_t dimension, size_t hessian_maximum_number_nonzeros): HessianEvaluation<MatrixType>(dimension,
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

template <class MatrixType>
class ConvexifiedExactHessianEvaluation : public HessianEvaluation<MatrixType> {
public:
   ConvexifiedExactHessianEvaluation(size_t dimension, size_t hessian_maximum_number_nonzeros, const std::string& linear_solver_name):
      HessianEvaluation<MatrixType>(dimension, hessian_maximum_number_nonzeros), linear_solver_(LinearSolverFactory<MatrixType>::create
      (linear_solver_name)) {
   }
   ~ConvexifiedExactHessianEvaluation() override = default;

   void compute(const Problem& problem, const std::vector<double>& primal_variables, double objective_multiplier,
         const std::vector<double>& constraint_multipliers) override {
      /* compute Hessian */
      problem.lagrangian_hessian(primal_variables, objective_multiplier, constraint_multipliers, this->hessian);
      this->evaluation_count++;
      DEBUG << "hessian before convexification: " << this->hessian;
      /* modify the inertia to make the problem strictly convex */
      this->hessian = this->modify_inertia(this->hessian, *this->linear_solver_);
   }

protected:
   std::unique_ptr<LinearSolver<MatrixType> > linear_solver_; /*!< Solver that computes the inertia */
};

template <class MatrixType>
class BFGSHessianEvaluation : public HessianEvaluation<MatrixType> {
public:
   explicit BFGSHessianEvaluation(size_t dimension, size_t hessian_maximum_number_nonzeros):
         HessianEvaluation<MatrixType>(dimension, hessian_maximum_number_nonzeros),
         previous_hessian_(dimension, hessian_maximum_number_nonzeros), previous_x_(dimension) {
   }

   ~BFGSHessianEvaluation() override = default;

   void compute(const Problem& problem, const std::vector<double>& primal_variables, double objective_multiplier,
         const std::vector<double>& constraint_multipliers) override {
      // the BFGS Hessian is already positive definite, do not convexify
      problem.lagrangian_hessian(primal_variables, objective_multiplier, constraint_multipliers, this->hessian);
      this->evaluation_count++;
   }

private:
   MatrixType previous_hessian_;
   std::vector<double> previous_x_;
};

template <class MatrixType>
class HessianEvaluationFactory {
public:
   static std::unique_ptr<HessianEvaluation<MatrixType> > create(const std::string& hessian_evaluation_method, size_t dimension,
         size_t hessian_maximum_number_nonzeros, bool convexify) {
      if (hessian_evaluation_method == "exact") {
         if (convexify) {
            return std::make_unique<ConvexifiedExactHessianEvaluation<MatrixType> >(dimension, hessian_maximum_number_nonzeros, "MA57");
         }
         else {
            return std::make_unique<ExactHessianEvaluation<MatrixType> >(dimension, hessian_maximum_number_nonzeros);
         }
      }
      else {
         throw std::invalid_argument("Hessian evaluation method " + hessian_evaluation_method + " does not exist");
      }
   }
};

#endif // HESSIANEVALUATION_H
