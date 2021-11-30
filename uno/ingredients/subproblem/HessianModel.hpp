#ifndef HESSIANMODEL_H
#define HESSIANMODEL_H

#include <memory>
#include <vector>
#include "optimization/Problem.hpp"
#include "optimization/Scaling.hpp"
#include "solvers/linear/LinearSolver.hpp"
#include "tools/Options.hpp"

class HessianModel {
public:
   // TODO handle padding
   HessianModel(size_t dimension, size_t hessian_maximum_number_nonzeros, const std::string& sparse_format);
   virtual ~HessianModel() = default;

   std::unique_ptr<SymmetricMatrix> hessian;
   size_t evaluation_count{0};

   virtual void evaluate(const Problem& problem, const Scaling& scaling, const std::vector<double>& primal_variables, double objective_multiplier,
         const std::vector<double>& constraint_multipliers) = 0;
   void adjust_number_variables(size_t number_variables);
};

// Exact Hessian
class ExactHessian : public HessianModel {
public:
   explicit ExactHessian(size_t dimension, size_t hessian_maximum_number_nonzeros, const Options& options);

   void evaluate(const Problem& problem, const Scaling& scaling, const std::vector<double>& primal_variables, double objective_multiplier,
         const std::vector<double>& constraint_multipliers) override;
};

// Hessian with convexification (inertia correction)
class ConvexifiedHessian : public HessianModel {
public:
   ConvexifiedHessian(size_t dimension, size_t hessian_maximum_number_nonzeros, const Options& options);

   void evaluate(const Problem& problem, const Scaling& scaling, const std::vector<double>& primal_variables, double objective_multiplier,
         const std::vector<double>& constraint_multipliers) override;

protected:
   std::unique_ptr<LinearSolver> linear_solver; /*!< Solver that computes the inertia */
   const double regularization_initial_value;

   void regularize(SymmetricMatrix& matrix);
};

// HessianModel factory
class HessianModelFactory {
public:
   static std::unique_ptr<HessianModel> create(const std::string& hessian_model, size_t dimension, size_t hessian_maximum_number_nonzeros,
         bool convexify, const Options& options);
};

#endif // HESSIANMODEL_H
