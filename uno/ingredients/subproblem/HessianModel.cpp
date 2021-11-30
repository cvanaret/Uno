#include "HessianModel.hpp"
#include "linear_algebra/SymmetricMatrixFactory.hpp"
#include "solvers/linear/LinearSolverFactory.hpp"

HessianModel::HessianModel(size_t dimension, size_t hessian_maximum_number_nonzeros, const std::string& sparse_format) :
      hessian(SymmetricMatrixFactory::create(sparse_format, dimension, hessian_maximum_number_nonzeros)) {
}

void HessianModel::adjust_number_variables(size_t number_variables) {
   // if the subproblem has more variables (slacks, elastic, ...) than the original problem, rectify the sparse representation
   // this assumes that all additional variables appear linearly in the functions
   for (size_t j = this->hessian->dimension; j < number_variables; j++) {
      this->hessian->finalize(j);
   }
   this->hessian->dimension = number_variables;
}

// Exact Hessian
ExactHessian::ExactHessian(size_t dimension, size_t hessian_maximum_number_nonzeros, const Options& options) :
   HessianModel(dimension, hessian_maximum_number_nonzeros, options.at("sparse_format")) {
}

void ExactHessian::evaluate(const Problem& problem, const Scaling& scaling, const std::vector<double>& primal_variables, double objective_multiplier,
      const std::vector<double>& constraint_multipliers) {
   // evaluate Lagrangian Hessian
   const double scaled_objective_multiplier = objective_multiplier*scaling.get_objective_scaling();
   std::vector<double> scaled_multipliers(problem.number_constraints);
   for (size_t j = 0; j < problem.number_constraints; j++) {
      scaled_multipliers[j] = scaling.get_constraint_scaling(j)*constraint_multipliers[j];
   }
   problem.evaluate_lagrangian_hessian(primal_variables, scaled_objective_multiplier, scaled_multipliers, *this->hessian);
   this->evaluation_count++;
}

// Convexified Hessian
ConvexifiedHessian::ConvexifiedHessian(size_t dimension, size_t hessian_maximum_number_nonzeros, const Options& options):
      HessianModel(dimension, hessian_maximum_number_nonzeros, options.at("sparse_format")),
      linear_solver(LinearSolverFactory::create(options.at("linear_solver"), dimension, hessian_maximum_number_nonzeros)),
      regularization_initial_value(stod(options.at("regularization_initial_value"))) {
}

void ConvexifiedHessian::evaluate(const Problem& problem, const Scaling& scaling, const std::vector<double>& primal_variables,
      double objective_multiplier, const std::vector<double>& constraint_multipliers) {
   // evaluate Lagrangian Hessian
   const double scaled_objective_multiplier = objective_multiplier*scaling.get_objective_scaling();
   std::vector<double> scaled_multipliers(problem.number_constraints);
   for (size_t j = 0; j < problem.number_constraints; j++) {
      scaled_multipliers[j] = scaling.get_constraint_scaling(j)*constraint_multipliers[j];
   }
   problem.evaluate_lagrangian_hessian(primal_variables, scaled_objective_multiplier, scaled_multipliers, *this->hessian);
   this->evaluation_count++;
   // regularize (only on the original variables) to make the problem strictly convex
   DEBUG << "hessian before convexification: " << *this->hessian;
   this->regularize(*this->hessian);
}

// Nocedal and Wright, p51
void ConvexifiedHessian::regularize(SymmetricMatrix& matrix) {
   //assert(size_block_to_regularize <= matrix.dimension && "The block to regularize is larger than the matrix");

   const double smallest_diagonal_entry = matrix.smallest_diagonal_entry();
   DEBUG << "The minimal diagonal entry of the matrix is " << matrix.smallest_diagonal_entry() << "\n";

   double regularization = (smallest_diagonal_entry <= 0.) ? this->regularization_initial_value - smallest_diagonal_entry : 0.;
   bool good_inertia = false;
   bool regularized = false;
   while (!good_inertia) {
      DEBUG << "Testing factorization with regularization factor " << regularization << "\n";
      if (0. < regularization) {
         // if regularization was already perform, delete the previous entries
         if (regularized) {
            for (size_t i = 0; i < matrix.dimension; i++) {
               matrix.pop();
            }
         }
         matrix.add_identity_multiple(regularization);
         regularized = true;
      }
      this->linear_solver->do_symbolic_factorization(matrix);
      this->linear_solver->do_numerical_factorization(matrix);

      if (this->linear_solver->matrix_is_positive_definite()) {
         good_inertia = true;
         DEBUG << "Factorization was a success with regularization factor " << regularization << "\n";
      }
      else {
         DEBUG << "rank: " << this->linear_solver->rank() << ", negative eigenvalues: " << this->linear_solver->number_negative_eigenvalues() << "\n";
         regularization = (regularization == 0.) ? this->regularization_initial_value : 2 * regularization;
         assert(regularization < std::numeric_limits<double>::infinity() && "The regularization coefficient diverged");
      }
   }
}

// Factory
std::unique_ptr<HessianModel> HessianModelFactory::create(const std::string& hessian_model, size_t dimension, size_t hessian_maximum_number_nonzeros,
      bool convexify, const Options& options) {
   if (hessian_model == "exact") {
      if (convexify) {
         return std::make_unique<ConvexifiedHessian>(dimension, hessian_maximum_number_nonzeros, options);
      }
      else {
         return std::make_unique<ExactHessian>(dimension, hessian_maximum_number_nonzeros, options);
      }
   }
   throw std::invalid_argument("Hessian model " + hessian_model + " does not exist");
}