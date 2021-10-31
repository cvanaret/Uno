#ifndef HESSIANMODEL_H
#define HESSIANMODEL_H

#include <memory>
#include <vector>
#include "optimization/Problem.hpp"
#include "solvers/linear/LinearSolver.hpp"

class HessianModel {
public:
   // TODO handle padding
   HessianModel(size_t dimension, size_t hessian_maximum_number_nonzeros, const std::string& sparse_format);
   virtual ~HessianModel() = default;

   std::unique_ptr<SymmetricMatrix> hessian;
   size_t evaluation_count{0};

   virtual void evaluate(const Problem& problem, const std::vector<double>& primal_variables, double objective_multiplier,
         const std::vector<double>& constraint_multipliers, size_t number_variables) = 0;

   void adjust_dimension(size_t number_variables) const;
   static void regularize(SymmetricMatrix& matrix, LinearSolver& linear_solver);
};

// Exact Hessian
class ExactHessian : public HessianModel {
public:
   explicit ExactHessian(size_t dimension, size_t hessian_maximum_number_nonzeros, const std::string& sparse_format);

   void evaluate(const Problem& problem, const std::vector<double>& primal_variables, double objective_multiplier,
         const std::vector<double>& constraint_multipliers, size_t number_variables) override;
};

// Hessian with convexification (inertia correction)
class ConvexifiedHessian : public HessianModel {
public:
   ConvexifiedHessian(size_t dimension, size_t hessian_maximum_number_nonzeros, const std::string& sparse_format,
         const std::string& linear_solver_name);

   void evaluate(const Problem& problem, const std::vector<double>& primal_variables, double objective_multiplier,
         const std::vector<double>& constraint_multipliers, size_t number_variables) override;

protected:
   std::unique_ptr<LinearSolver> linear_solver; /*!< Solver that computes the inertia */
};

// HessianModel factory
class HessianModelFactory {
public:
   static std::unique_ptr<HessianModel> create(const std::string& hessian_model, size_t dimension, size_t hessian_maximum_number_nonzeros,
         const std::string& sparse_format, bool convexify);
};

#endif // HESSIANMODEL_H
