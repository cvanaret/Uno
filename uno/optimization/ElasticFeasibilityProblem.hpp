#ifndef UNO_ELASTICREFORMULATION_H
#define UNO_ELASTICREFORMULATION_H

#include <vector>
#include "Problem.hpp"
#include "ingredients/constraint_relaxation/ElasticVariables.hpp"

class ElasticFeasibilityProblem: public Problem {
public:
   ElasticFeasibilityProblem(const Problem& original_problem, double objective_multiplier, double elastic_objective_coefficient,
         double proximal_coefficient);

   [[nodiscard]] double get_variable_lower_bound(size_t i) const override;
   [[nodiscard]] double get_variable_upper_bound(size_t i) const override;
   [[nodiscard]] double get_constraint_lower_bound(size_t j) const override;
   [[nodiscard]] double get_constraint_upper_bound(size_t j) const override;

   [[nodiscard]] double evaluate_objective(const std::vector<double>& x) const override;
   void evaluate_objective_gradient(const std::vector<double>& x, SparseVector<double>& gradient) const override;
   void evaluate_constraints(const std::vector<double>& x, std::vector<double>& constraints) const override;
   void evaluate_constraint_jacobian(const std::vector<double>& x, std::vector<SparseVector<double>>& constraint_jacobian) const override;
   void evaluate_lagrangian_hessian(const std::vector<double>& x, double objective_multiplier, const std::vector<double>& multipliers,
         SymmetricMatrix& hessian) const override;

   [[nodiscard]] ConstraintType get_variable_status(size_t i) const override;
   [[nodiscard]] FunctionType get_constraint_type(size_t j) const override;
   [[nodiscard]] ConstraintType get_constraint_status(size_t j) const override;
   [[nodiscard]] size_t get_hessian_maximum_number_nonzeros() const override;

   void get_initial_primal_point(std::vector<double>& x) const override;
   void get_initial_dual_point(std::vector<double>& multipliers) const override;

   void set_objective_multiplier(double new_objective_multiplier);
   void set_proximal_coefficient(double new_proximal_coefficient);
   void set_proximal_reference_point(const std::vector<double>& new_proximal_reference_point);

protected:
   const Problem& original_problem;
   double objective_multiplier;
   // elastic variables
   ElasticVariables elastic_variables;
   double elastic_objective_coefficient;
   // proximal term
   double proximal_coefficient;
   std::vector<double> proximal_reference_point;

   [[nodiscard]] double compute_elastic_residual(const std::vector<double>& x) const;
   [[nodiscard]] double get_proximal_weight(size_t i) const;
};

#endif // UNO_ELASTICREFORMULATION_H