#include "Evaluations.hpp"

namespace uno {
   Evaluations::Evaluations(size_t number_variables, size_t number_constraints):
      constraints(number_constraints),
      objective_gradient(number_variables) {
   }

   void evaluate_objective(const Model& model, const Vector<double>& primals) {

   }

   void evaluate_constraints(const Model& model, const Vector<double>& primals) {

   }

   void evaluate_objective_gradient(const Model& model, const Vector<double>& primals) {

   }
} // namespace