// Copyright (c) 2026 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <algorithm>
#include <cassert>
#include <cmath>
#include <stdexcept>
#include "NLModel.hpp"
#include "optimization/EvaluationErrors.hpp"
#include "optimization/Result.hpp"
#include "symbolic/Range.hpp"
#include "Uno.hpp"

namespace uno {
   namespace {
      cppasl::NlModel read_nl_model(const std::string& file_name) {
         cppasl::NlModel nl_model = cppasl::NlModel::read_file(file_name);
         const cppasl::NlHeader& header = nl_model.header();
         const int number_discrete_variables = header.number_binary_variables + header.number_integer_variables +
            header.number_nonlinear_integer_variables_in_both + header.number_nonlinear_integer_variables_in_constraints +
            header.number_nonlinear_integer_variables_in_objectives;
         if (0 < number_discrete_variables) {
            throw std::runtime_error("Uno does not support discrete variables");
         }
         return nl_model;
      }

      bool all_finite(const double* values, size_t size) {
         return std::all_of(values, values + size, [](double value) { return std::isfinite(value); });
      }
   } // namespace

   NLModel::NLModel(const std::string& file_name) : NLModel(file_name, read_nl_model(file_name)) {
   }

   NLModel::NLModel(const std::string& file_name, cppasl::NlModel&& nl_model) :
         Model(file_name, static_cast<size_t>(nl_model.number_variables()), static_cast<size_t>(nl_model.number_constraints()),
            (0 < nl_model.number_objectives() && nl_model.is_maximization(0)) ? -1. : 1.,
            NLModel::lagrangian_sign_convention, 0),
         nl_model(std::move(nl_model)),
         workspace(this->nl_model),
         has_objective(0 < this->nl_model.number_objectives()),
         jacobian_row_indices(this->nl_model.number_jacobian_nonzeros()),
         jacobian_column_indices(this->nl_model.number_jacobian_nonzeros()),
         jacobian_values(this->nl_model.number_jacobian_nonzeros()),
         cppasl_multipliers(this->number_constraints),
         linear_constraints_collection(this->linear_constraints),
         nonlinear_constraints_collection(this->nonlinear_constraints),
         equality_constraints_collection(this->equality_constraints),
         inequality_constraints_collection(this->inequality_constraints) {
      // Jacobian sparsity (CSR order, as the values)
      std::copy(this->nl_model.jacobian_row_indices().begin(), this->nl_model.jacobian_row_indices().end(),
         this->jacobian_row_indices.begin());
      std::copy(this->nl_model.jacobian_column_indices().begin(), this->nl_model.jacobian_column_indices().end(),
         this->jacobian_column_indices.begin());

      // linear/nonlinear constraints
      for (size_t constraint_index: Range(this->number_constraints)) {
         if (this->nl_model.constraint_structure(static_cast<int>(constraint_index)) == cppasl::FunctionStructure::Linear) {
            this->linear_constraints.emplace_back(constraint_index);
         }
         else {
            this->nonlinear_constraints.emplace_back(constraint_index);
         }
      }
      Model::find_fixed_variables(this->fixed_variables);
      Model::partition_constraints(this->equality_constraints, this->inequality_constraints);
      this->problem_type = this->determine_problem_type();
   }

   ProblemType NLModel::get_problem_type() const {
      return this->problem_type;
   }

   bool NLModel::has_jacobian_operator() const {
      return false; // available, but each product re-evaluates the Jacobian
   }

   bool NLModel::has_jacobian_transposed_operator() const {
      return false; // available, but each product re-evaluates the Jacobian
   }

   bool NLModel::has_hessian_operator() const {
      return true;
   }

   bool NLModel::has_hessian_matrix() const {
      return true;
   }

   double NLModel::evaluate_objective(const Vector<double>& x) const {
      const double objective = this->has_objective ?
         this->optimization_sense * this->nl_model.evaluate_objective(this->workspace, x.data()) : 0.;
      if (!std::isfinite(objective)) {
         throw FunctionEvaluationError();
      }
      ++this->number_model_evaluations.objective;
      return objective;
   }

   void NLModel::evaluate_constraints(const Vector<double>& x, Vector<double>& constraints) const {
      this->nl_model.evaluate_constraints(this->workspace, x.data(), constraints.data());
      if (!all_finite(constraints.data(), this->number_constraints)) {
         throw FunctionEvaluationError();
      }
      ++this->number_model_evaluations.constraints;
   }

   void NLModel::evaluate_objective_gradient(const Vector<double>& x, Vector<double>& gradient) const {
      if (!this->has_objective) {
         std::fill_n(gradient.begin(), this->number_variables, 0.);
         return;
      }
      this->nl_model.evaluate_objective_gradient(this->workspace, x.data(), gradient.data());
      if (!all_finite(gradient.data(), this->number_variables)) {
         throw GradientEvaluationError();
      }
      if (this->optimization_sense != 1.) {
         gradient.scale(this->optimization_sense);
      }
      ++this->number_model_evaluations.objective_gradient;
   }

   View<const uno_int> NLModel::get_jacobian_row_indices() const {
      return this->jacobian_row_indices.view();
   }

   View<const uno_int> NLModel::get_jacobian_column_indices() const {
      return this->jacobian_column_indices.view();
   }

   void NLModel::compute_hessian_sparsity(View<uno_int> row_indices, View<uno_int> column_indices, uno_int solver_indexing) const {
      // lower triangle, CSC order
      const std::vector<int>& nl_row_indices = this->nl_model.hessian_row_indices();
      const std::vector<int>& nl_column_indices = this->nl_model.hessian_column_indices();
      for (size_t nonzero_index: Range(this->nl_model.number_hessian_nonzeros())) {
         row_indices[nonzero_index] = static_cast<uno_int>(nl_row_indices[nonzero_index]) + solver_indexing;
         column_indices[nonzero_index] = static_cast<uno_int>(nl_column_indices[nonzero_index]) + solver_indexing;
      }
   }

   void NLModel::evaluate_jacobian(const Vector<double>& x, double* jacobian_values) const {
      this->nl_model.evaluate_jacobian(this->workspace, x.data(), jacobian_values);
      if (!all_finite(jacobian_values, this->number_jacobian_nonzeros())) {
         throw GradientEvaluationError();
      }
      ++this->number_model_evaluations.jacobian;
   }

   void NLModel::evaluate_lagrangian_hessian(const Vector<double>& x, double objective_multiplier, const Vector<double>& multipliers,
         View<double> hessian_values) const {
      objective_multiplier = this->has_objective ? this->optimization_sense * objective_multiplier : 0.;
      this->nl_model.evaluate_lagrangian_hessian(this->workspace, x.data(), objective_multiplier,
         this->negated_multipliers(multipliers), hessian_values.data());
      if (!all_finite(hessian_values.data(), this->number_hessian_nonzeros())) {
         throw HessianEvaluationError();
      }
      ++this->number_model_evaluations.hessian;
   }

   void NLModel::compute_jacobian_vector_product(const double* x, const double* vector, double* result) const {
      this->nl_model.evaluate_jacobian(this->workspace, x, this->jacobian_values.data());
      this->nl_model.multiply_jacobian(this->jacobian_values.data(), vector, result);
   }

   void NLModel::compute_jacobian_transposed_vector_product(const double* x, const double* vector, double* result) const {
      this->nl_model.evaluate_jacobian(this->workspace, x, this->jacobian_values.data());
      this->nl_model.multiply_jacobian_transpose(this->jacobian_values.data(), vector, result);
   }

   void NLModel::compute_hessian_vector_product(View<const double> x, View<const double> vector, double objective_multiplier,
         const Vector<double>& multipliers, View<double> result) const {
      objective_multiplier = this->has_objective ? this->optimization_sense * objective_multiplier : 0.;
      this->nl_model.evaluate_lagrangian_hessian_vector_product(this->workspace, x.data(), objective_multiplier,
         this->negated_multipliers(multipliers), vector.data(), result.data());
      if (!all_finite(result.data(), this->number_variables)) {
         throw HessianEvaluationError();
      }
   }

   const std::vector<double>& NLModel::get_variables_lower_bounds() const {
      return this->nl_model.variable_lower_bounds();
   }

   const std::vector<double>& NLModel::get_variables_upper_bounds() const {
      return this->nl_model.variable_upper_bounds();
   }

   const Vector<size_t>& NLModel::get_fixed_variables() const {
      return this->fixed_variables;
   }

   const std::vector<double>& NLModel::get_constraints_lower_bounds() const {
      return this->nl_model.constraint_lower_bounds();
   }

   const std::vector<double>& NLModel::get_constraints_upper_bounds() const {
      return this->nl_model.constraint_upper_bounds();
   }

   const Collection<size_t>& NLModel::get_equality_constraints() const {
      return this->equality_constraints_collection;
   }

   const Collection<size_t>& NLModel::get_inequality_constraints() const {
      return this->inequality_constraints_collection;
   }

   const Collection<size_t>& NLModel::get_linear_constraints() const {
      return this->linear_constraints_collection;
   }

   const Collection<size_t>& NLModel::get_nonlinear_constraints() const {
      return this->nonlinear_constraints_collection;
   }

   void NLModel::initial_primal_point(Vector<double>& x) const {
      assert(x.size() >= this->number_variables);
      std::copy_n(this->nl_model.initial_primal_point().begin(), this->number_variables, x.begin());
   }

   void NLModel::initial_dual_point(Vector<double>& multipliers) const {
      assert(multipliers.size() >= this->number_constraints);
      std::copy_n(this->nl_model.initial_dual_point().begin(), this->number_constraints, multipliers.begin());
   }

   void NLModel::postprocess_solution(Iterate& /*iterate*/, Evaluations& /*evaluations*/) const {
   }

   void NLModel::write_solution_to_file(Result& result) const {
      int solve_code = 400; // limit
      switch (result.solution_status) {
         case SolutionStatus::FEASIBLE_KKT_POINT: solve_code = 0; break;
         case SolutionStatus::FEASIBLE_SMALL_STEP: solve_code = 100; break;
         case SolutionStatus::INFEASIBLE_STATIONARY_POINT: solve_code = 200; break;
         case SolutionStatus::DIVERGING_ITERATE: case SolutionStatus::UNBOUNDED_OBJECTIVE: solve_code = 300; break;
         case SolutionStatus::INFEASIBLE_SMALL_STEP: solve_code = 500; break;
         default: break;
      }
      std::string message = "Uno ";
      message.append(Uno::current_version()).append(": ").append(solution_status_to_message(result.solution_status));

      // bound duals as suffixes
      const auto to_suffix = [&](const char* name, const Vector<double>& values) {
         return cppasl::Suffix{name, cppasl::SuffixTarget::Variables, true,
            std::vector<double>(values.begin(), values.begin() + static_cast<std::ptrdiff_t>(this->number_variables))};
      };
      const std::vector<cppasl::Suffix> suffixes{to_suffix("lower_bound_duals", result.lower_bound_dual_solution),
         to_suffix("upper_bound_duals", result.upper_bound_dual_solution)};

      // flip the signs of the constraint multipliers (different Lagrangian convention between AMPL and Uno)
      const Vector<double> ampl_dual_solution = -result.constraint_dual_solution;
      // as ASL: file stub (without .nl) + .sol
      std::string stub = this->name;
      if (4 <= stub.size() && stub.compare(stub.size() - 3, 3, ".nl") == 0) {
         stub.resize(stub.size() - 3);
      }
      this->nl_model.write_solution(stub + ".sol", message, result.primal_solution.data(), ampl_dual_solution.data(),
         solve_code, suffixes);
   }

   size_t NLModel::number_jacobian_nonzeros() const {
      return this->nl_model.number_jacobian_nonzeros();
   }

   size_t NLModel::number_hessian_nonzeros() const {
      return this->nl_model.number_hessian_nonzeros();
   }

   size_t NLModel::number_model_objective_evaluations() const {
      return this->number_model_evaluations.objective;
   }

   size_t NLModel::number_model_constraints_evaluations() const {
      return this->number_model_evaluations.constraints;
   }

   size_t NLModel::number_model_objective_gradient_evaluations() const {
      return this->number_model_evaluations.objective_gradient;
   }

   size_t NLModel::number_model_jacobian_evaluations() const {
      return this->number_model_evaluations.jacobian;
   }

   size_t NLModel::number_model_hessian_evaluations() const {
      return this->number_model_evaluations.hessian;
   }

   void NLModel::reset_number_evaluations() const {
      this->number_model_evaluations.reset();
   }

   // private member functions

   // exact structure detection (no convexity check, unlike cppasl::NlModel::classify())
   ProblemType NLModel::determine_problem_type() const {
      if (!this->nonlinear_constraints.empty()) {
         return ProblemType::NONLINEAR;
      }
      const cppasl::FunctionStructure objective_structure = this->has_objective ? this->nl_model.objective_structure(0) :
         cppasl::FunctionStructure::Linear;
      switch (objective_structure) {
         case cppasl::FunctionStructure::Linear: return ProblemType::LINEAR;
         case cppasl::FunctionStructure::Quadratic: return ProblemType::QUADRATIC;
         default: return ProblemType::NONLINEAR;
      }
   }

   const double* NLModel::negated_multipliers(const Vector<double>& multipliers) const {
      for (size_t constraint_index: Range(this->number_constraints)) {
         this->cppasl_multipliers[constraint_index] = -multipliers[constraint_index];
      }
      return this->cppasl_multipliers.data();
   }
} // namespace
