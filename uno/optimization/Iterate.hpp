#ifndef ITERATE_H
#define ITERATE_H

#include <ostream>
#include <vector>
#include "Problem.hpp"
#include "Constraint.hpp"

enum TerminationStatus {
   NOT_OPTIMAL = 0, KKT_POINT, /* feasible stationary point */
   FJ_POINT, /* infeasible stationary point */
   FEASIBLE_SMALL_STEP, INFEASIBLE_SMALL_STEP
};

struct Errors {
   double constraints;
   double KKT;
   double FJ;
   double complementarity;
};

struct ProgressMeasures {
   double infeasibility;
   double objective;
};

/*! \class Iterate
 * \brief Optimization iterate
 *
 *  Point and its evaluations during an optimization process
 */
class Iterate {
public:
   /*!
    *  Constructor
    */
   Iterate(const std::vector<double>& x, const Multipliers& multipliers);
   Iterate(size_t number_variables, size_t number_constraints);

   std::vector<double> x; /*!< \f$\mathbb{R}^n\f$ primal variables */
   Multipliers multipliers; /*!< \f$\mathbb{R}^n\f$ Lagrange multipliers/dual variables */
   static int number_eval_objective;
   static int number_eval_constraints;
   static int number_eval_jacobian;

   // functions
   double objective{std::numeric_limits<double>::infinity()}; /*!< Objective value */
   bool is_objective_computed{false};

   std::vector<double> constraints; /*!< Constraint values (size \f$m)\f$ */
   bool are_constraints_computed{false};

   SparseVector<double> objective_gradient; /*!< Sparse Jacobian of the evaluate_objective */
   bool is_objective_gradient_computed{false}; /*!< Flag that indicates if the evaluate_objective gradient has already been computed */

   std::vector<SparseVector<double>> constraints_jacobian; /*!< Sparse Jacobian of the constraints */
   bool is_constraints_jacobian_computed{false}; /*!< Flag that indicates if the constraint Jacobian has already been computed */

   // residuals
   Errors errors{0., 0., 0., 0.};
   ProgressMeasures progress{0., 0.};

   void compute_objective(const Problem& problem);
   void compute_constraints(const Problem& problem);
   void compute_objective_gradient(const Problem& problem);
   void compute_constraints_jacobian(const Problem& problem);
   std::vector<double> lagrangian_gradient(const Problem& problem, double objective_multiplier, const Multipliers& multipliers);

   friend std::ostream& operator<<(std::ostream& stream, const Iterate& iterate);
};

void reset(Iterate& iterate);

#endif // ITERATE_H
