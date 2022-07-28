// Copyright (c) 2022 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_HESSIANMODEL_H
#define UNO_HESSIANMODEL_H

#include <memory>
#include <vector>
#include "ingredients/constraint_relaxation_strategy/NonlinearProblem.hpp"
#include "solvers/linear/LinearSolver.hpp"
#include "tools/Options.hpp"

class HessianModel {
public:
   HessianModel(size_t dimension, size_t maximum_number_nonzeros, const std::string& sparse_format, bool use_regularization);
   virtual ~HessianModel() = default;

   std::unique_ptr<SymmetricMatrix> hessian;
   size_t evaluation_count{0};

   virtual void evaluate(const NonlinearProblem& problem, const std::vector<double>& primal_variables, const std::vector<double>& constraint_multipliers) = 0;
};

// Exact Hessian
class ExactHessian : public HessianModel {
public:
   explicit ExactHessian(size_t dimension, size_t maximum_number_nonzeros, const Options& options);

   void evaluate(const NonlinearProblem& problem, const std::vector<double>& primal_variables, const std::vector<double>& constraint_multipliers) override;
};

// Hessian with convexification (inertia correction)
class ConvexifiedHessian : public HessianModel {
public:
   ConvexifiedHessian(size_t dimension, size_t maximum_number_nonzeros, const Options& options);

   void evaluate(const NonlinearProblem& problem, const std::vector<double>& primal_variables, const std::vector<double>& constraint_multipliers) override;

protected:
   std::unique_ptr<LinearSolver> linear_solver; /*!< Solver that computes the inertia */
   const double regularization_initial_value{};
   const double regularization_increase_factor{};

   void regularize(SymmetricMatrix& hessian, size_t number_original_variables);
};

// HessianModel factory
class HessianModelFactory {
public:
   static std::unique_ptr<HessianModel> create(const std::string& hessian_model, size_t dimension, size_t maximum_number_nonzeros,
         bool convexify, const Options& options);
};

#endif // UNO_HESSIANMODEL_H
