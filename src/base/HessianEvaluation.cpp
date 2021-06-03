#include <exception>
#include "HessianEvaluation.hpp"
#include "LinearSolverFactory.hpp"
#include "Vector.hpp"

HessianEvaluation::HessianEvaluation(int dimension, int hessian_maximum_number_nonzeros) : dimension(dimension),
      hessian(dimension, hessian_maximum_number_nonzeros, 1) {
}

/* Exact Hessian */

ExactHessianEvaluation::ExactHessianEvaluation(int dimension, int hessian_maximum_number_nonzeros) : HessianEvaluation(dimension,
      hessian_maximum_number_nonzeros) {
}

void ExactHessianEvaluation::compute(Problem& problem, Iterate& iterate, double objective_multiplier,
      std::vector<double>& constraint_multipliers) {
   /* compute Hessian */
   problem.lagrangian_hessian(iterate.x, objective_multiplier, constraint_multipliers, this->hessian);
   return;
}

/* Exact Hessian with inertia control */

ExactHessianInertiaControlEvaluation::ExactHessianInertiaControlEvaluation(int dimension, int hessian_maximum_number_nonzeros,
      std::string linear_solver_name) : HessianEvaluation(dimension, hessian_maximum_number_nonzeros),
      linear_solver_(LinearSolverFactory::create(linear_solver_name)) {
}

void ExactHessianInertiaControlEvaluation::compute(Problem& problem, Iterate& iterate, double objective_multiplier,
      std::vector<double>& constraint_multipliers) {
   /* compute Hessian */
   problem.lagrangian_hessian(iterate.x, objective_multiplier, constraint_multipliers, this->hessian);
   DEBUG << "hessian before convexification: " << this->hessian;
   /* modify the inertia to make the problem strictly convex */
   //this->hessian = this->modify_inertia(this->hessian, *this->linear_solver_);
   return;
}

CSCMatrix HessianEvaluation::modify_inertia(CSCMatrix& hessian, LinearSolver& linear_solver) {
   double beta = 1e-4;

   // Nocedal and Wright, p51
   double smallest_diagonal_entry = hessian.smallest_diagonal_entry();
   DEBUG << "The minimal diagonal entry of the Hessian is " << hessian.smallest_diagonal_entry() << "\n";

   double inertia = 0.;
   double previous_inertia = 0.;
   if (smallest_diagonal_entry <= 0.) {
      inertia = beta - smallest_diagonal_entry;
   }

   if (0. < inertia) {
      hessian = hessian.add_identity_multiple(inertia - previous_inertia);
   }
   COOMatrix coo_hessian = hessian.to_COO();
   DEBUG << "Testing factorization with inertia term " << inertia << "\n";
   linear_solver.do_symbolic_factorization(coo_hessian);

   bool good_inertia = false;
   while (!good_inertia) {
      DEBUG << linear_solver.number_negative_eigenvalues() << " negative eigenvalues\n";
      if (!linear_solver.matrix_is_singular() && linear_solver.number_negative_eigenvalues() == 0) {
         good_inertia = true;
         DEBUG << "Factorization was a success with inertia " << inertia << "\n";
      }
      else {
         previous_inertia = inertia;
         if (inertia == 0.) {
            inertia = beta;
         }
         else {
            inertia *= 2.;
         }
         hessian = hessian.add_identity_multiple(inertia - previous_inertia);
         coo_hessian = hessian.to_COO();
         DEBUG << "Testing factorization with inertia term " << inertia << "\n";
         linear_solver.do_symbolic_factorization(coo_hessian);
      }
   }
   return hessian;
}

/* BFGS Hessian */

BFGSHessianEvaluation::BFGSHessianEvaluation(int dimension, int hessian_maximum_number_nonzeros) :
   HessianEvaluation(dimension, hessian_maximum_number_nonzeros),
   previous_hessian_(dimension, hessian_maximum_number_nonzeros, 1), previous_x_(dimension) {
}

void BFGSHessianEvaluation::compute(Problem& problem, Iterate& iterate, double objective_multiplier,
      std::vector<double>& constraint_multipliers) {
   // the BFGS Hessian is already positive definite, do not convexify
   problem.lagrangian_hessian(iterate.x, objective_multiplier, constraint_multipliers, this->hessian);
   return;
}

/* Factory */

std::unique_ptr<HessianEvaluation>
HessianEvaluationFactory::create(std::string hessian_evaluation_method, int dimension, int hessian_maximum_number_nonzeros,
      bool convexify) {
   if (hessian_evaluation_method == "exact") {
      if (convexify) {
         return std::make_unique<ExactHessianInertiaControlEvaluation>(dimension, hessian_maximum_number_nonzeros, "MA57");
      }
      else {
         return std::make_unique<ExactHessianEvaluation>(dimension, hessian_maximum_number_nonzeros);
      }
   }
      //    else if (hessian_evaluation_method == "BFGS") {
      //        return std::make_unique<BFGSHessianEvaluation>(dimension, hessian_maximum_number_nonzeros);
      //    }
   else {
      throw std::invalid_argument("Hessian evaluation method " + hessian_evaluation_method + " does not exist");
   }
}
