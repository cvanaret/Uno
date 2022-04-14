#ifndef UNO_AUGMENTEDSYSTEM_H
#define UNO_AUGMENTEDSYSTEM_H

#include <memory>
#include "linear_algebra/SymmetricMatrix.hpp"
#include "solvers/linear/LinearSolver.hpp"
#include "optimization/Model.hpp"

struct UnstableRegularization : public std::exception {

   [[nodiscard]] const char* what() const throw() override {
      return "The inertia correction got unstable (delta_w > threshold)";
   }
};

class AugmentedSystem {
public:
   std::unique_ptr<SymmetricMatrix> matrix;
   std::vector<double> rhs;
   std::vector<double> solution;

   AugmentedSystem(const std::string& sparse_format, size_t max_dimension, size_t max_number_non_zeros, double regularization_failure_threshold);
   void solve(LinearSolver& linear_solver);
   void factorize_matrix(const Model& model, LinearSolver& linear_solver);
   void regularize_matrix(const Model& model, LinearSolver& linear_solver, size_t size_first_block, size_t size_second_block,
         double constraint_regularization_parameter);

protected:
   size_t number_factorizations{0};
   double previous_regularization_first_block{0.};
   const double regularization_failure_threshold;
};

#endif // UNO_AUGMENTEDSYSTEM_H
