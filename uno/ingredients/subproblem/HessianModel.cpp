#include "HessianModel.hpp"
#include "linear_algebra/SymmetricMatrixFactory.hpp"
#include "solvers/linear/LinearSolverFactory.hpp"

HessianModel::HessianModel(size_t dimension, size_t hessian_maximum_number_nonzeros, const std::string& sparse_format) :
      hessian(SymmetricMatrixFactory::create(sparse_format, dimension, hessian_maximum_number_nonzeros)) {
}

void HessianModel::adjust_dimension(size_t number_variables) const {
   // if the subproblem has more variables (slacks, elastic, ...) than the original problem, rectify the sparse representation
   // this assumes that all additional variables appear linearly in the functions
   for (size_t j = this->hessian->dimension; j < number_variables; j++) {
      this->hessian->finalize(j);
   }
   this->hessian->dimension = number_variables;
}

// Nocedal and Wright, p51
void HessianModel::regularize(SymmetricMatrix& matrix, LinearSolver& linear_solver) {
   const double beta = 1e-4;
   const double smallest_diagonal_entry = matrix.smallest_diagonal_entry();
   DEBUG << "The minimal diagonal entry of the matrix is " << matrix.smallest_diagonal_entry() << "\n";

   double regularization = (smallest_diagonal_entry <= 0.) ? beta - smallest_diagonal_entry : 0.;
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
      linear_solver.do_symbolic_factorization(matrix.dimension, matrix);
      linear_solver.do_numerical_factorization(matrix.dimension, matrix);

      if (linear_solver.matrix_is_positive_definite()) {
         good_inertia = true;
         DEBUG << "Factorization was a success with regularization factor " << regularization << "\n";
      }
      else {
         DEBUG << "rank: " << linear_solver.rank() << ", negative eigenvalues: " << linear_solver.number_negative_eigenvalues() << "\n";
         regularization = (regularization == 0.) ? beta : 2 * regularization;
      }
   }
}

// Exact Hessian
ExactHessian::ExactHessian(size_t dimension, size_t hessian_maximum_number_nonzeros, const std::string& sparse_format) :
   HessianModel(dimension, hessian_maximum_number_nonzeros, sparse_format) {
}

void ExactHessian::evaluate(const Problem& problem, const std::vector<double>& primal_variables, double objective_multiplier,
      const std::vector<double>& constraint_multipliers, size_t number_variables) {
   // compute Hessian
   problem.evaluate_lagrangian_hessian(primal_variables, objective_multiplier, constraint_multipliers, *this->hessian);
   this->adjust_dimension(number_variables);
   this->evaluation_count++;
}

// Convexified Hessian
ConvexifiedHessian::ConvexifiedHessian(size_t dimension, size_t hessian_maximum_number_nonzeros, const std::string& sparse_format,
      const std::string& linear_solver_name) : HessianModel(dimension, hessian_maximum_number_nonzeros, sparse_format),
      linear_solver(LinearSolverFactory::create(linear_solver_name, dimension, hessian_maximum_number_nonzeros)) {
}

void ConvexifiedHessian::evaluate(const Problem& problem, const std::vector<double>& primal_variables, double objective_multiplier,
      const std::vector<double>& constraint_multipliers, size_t number_variables) {
   // compute Hessian
   problem.evaluate_lagrangian_hessian(primal_variables, objective_multiplier, constraint_multipliers, *this->hessian);
   this->adjust_dimension(number_variables);
   this->evaluation_count++;
   // make the problem strictly convex
   DEBUG << "hessian before convexification: " << *this->hessian;
   HessianModel::regularize(*this->hessian, *this->linear_solver);
}

// Factory
std::unique_ptr<HessianModel> HessianModelFactory::create(const std::string& hessian_model, size_t dimension, size_t hessian_maximum_number_nonzeros,
      const std::string& sparse_format, bool convexify) {
   if (hessian_model == "exact") {
      if (convexify) {
         //const std::string& linear_solver_name = options.at("linear_solver");
         return std::make_unique<ConvexifiedHessian>(dimension, hessian_maximum_number_nonzeros, sparse_format, "MA57");
      }
      else {
         return std::make_unique<ExactHessian>(dimension, hessian_maximum_number_nonzeros, sparse_format);
      }
   }
   throw std::invalid_argument("Hessian model " + hessian_model + " does not exist");
}