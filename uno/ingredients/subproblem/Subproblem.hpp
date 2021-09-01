#ifndef SUBPROBLEM_H
#define SUBPROBLEM_H

#include <cmath>
#include <vector>
#include <array>
#include <memory>
#include "optimization/Problem.hpp"
#include "optimization/Iterate.hpp"
#include "optimization/Constraint.hpp"
#include "Direction.hpp"
#include "linear_algebra/Vector.hpp"
#include "solvers/linear/LinearSolver.hpp"
#include "tools/Statistics.hpp"

template <size_t NUMBER_DIRECTIONS>
class SubproblemTest {
public:
   static constexpr size_t number_directions = NUMBER_DIRECTIONS;

   virtual std::array<Direction, NUMBER_DIRECTIONS> compute_directions() = 0;
};

class SQPTest: public SubproblemTest<1> {
public:
   std::array<Direction, 1> compute_directions() override;
};

class SLQPTest: public SubproblemTest<2> {
public:
   std::array<Direction, 2> compute_directions() override;
};

/*! \class Subproblem
 * \brief Subproblem
 *
 *  Local approximation of a nonlinear optimization problem (virtual class) 
 */
class Subproblem {
public:
   Subproblem(size_t number_variables, size_t number_constraints);
   virtual ~Subproblem() = default;

   // virtual methods implemented by subclasses
   virtual Iterate generate_initial_iterate(Statistics& statistics, const Problem& problem, std::vector<double>& x, Multipliers& multipliers);
   virtual void generate(const Problem& problem, Iterate& current_iterate, double objective_multiplier, double trust_region_radius) = 0;
   virtual void update_objective_multiplier(const Problem& problem, const Iterate& current_iterate, double objective_multiplier) = 0;

   // direction computation
   virtual Direction solve(Statistics& statistics, const Problem& problem, Iterate& current_iterate) = 0;
   virtual Direction compute_second_order_correction(const Problem& problem, Iterate& trial_iterate);

   // globalization metrics
   virtual double compute_predicted_reduction(const Direction& direction, double step_length) const = 0;
   virtual void compute_progress_measures(const Problem& problem, Iterate& iterate);
   virtual void register_accepted_iterate(Iterate& iterate);

   virtual int get_hessian_evaluation_count() const = 0;
   virtual void set_initial_point(const std::vector<double>& initial_point) = 0;

   // available methods
   // feasibility subproblem
   void compute_feasibility_linear_objective(const Iterate& current_iterate, const ConstraintPartition& constraint_partition);
   void generate_feasibility_bounds(const Problem& problem, const std::vector<double>& current_constraints, const ConstraintPartition&
   constraint_partition);
   static double push_variable_to_interior(double variable_value, const Range& variable_bounds);
   void set_constraints_bounds(const Problem& problem, const std::vector<double>& current_constraints);

   template<class SparseSymmetricMatrix>
   static void compute_least_square_multipliers(const Problem& problem, SparseSymmetricMatrix& matrix, std::vector<double>& rhs, LinearSolver
   <SparseSymmetricMatrix>& solver, Iterate& current_iterate, std::vector<double>& multipliers, double multipliers_max_size = 1e3);

   static double compute_first_order_error(const Problem& problem, Iterate& iterate, double objective_multiplier);
   void compute_errors(const Problem& problem, Iterate& iterate, double objective_multiplier) const;

   double compute_complementarity_error(const Problem& problem, Iterate& iterate, const Multipliers& multipliers) const;

   const size_t number_variables;
   const size_t number_constraints;
   // when the subproblem is reformulated (e.g. when slacks are introduced), the bounds may be altered
   std::vector <Range> variables_bounds;
   std::vector<double> constraints_multipliers;
   SparseVector<double> objective_gradient;
   std::vector <SparseVector<double>> constraints_jacobian;
   std::vector <Range> constraints_bounds;
   // Hessian is optional and depends on the subproblem

   int number_subproblems_solved;
   // when the parameterization of the subproblem (e.g. penalty or barrier parameter) is updated, signal it
   bool subproblem_definition_changed;

protected:
   virtual void set_variables_bounds(const Problem& problem, const Iterate& current_iterate, double trust_region_radius);
};

// compute a least-square approximation of the multipliers by solving a linear system (uses existing linear system)
template<class SparseSymmetricMatrix>
inline void Subproblem::compute_least_square_multipliers(const Problem& problem, SparseSymmetricMatrix& matrix, std::vector<double>& rhs, LinearSolver
<SparseSymmetricMatrix>& solver, Iterate& current_iterate, std::vector<double>& multipliers, double multipliers_max_size) {
   current_iterate.compute_objective_gradient(problem);
   current_iterate.compute_constraints_jacobian(problem);

   /******************************/
   /* build the symmetric matrix */
   /******************************/
   matrix.reset();

   /* identity block */
   for (size_t i = 0; i < current_iterate.x.size(); i++) {
      matrix.insert(1., i, i);
   }
   /* Jacobian of general constraints */
   const size_t n = current_iterate.x.size();
   for (size_t j = 0; j < problem.number_constraints; j++) {
      current_iterate.constraints_jacobian[j].for_each([&](size_t i, double derivative) {
         matrix.insert(derivative, i, n + j);
      });
   }
   DEBUG << "KKT matrix for least-square multipliers:\n" << matrix << "\n";

   /********************************/
   /* generate the right-hand side */
   /********************************/
   clear(rhs);

   /* objective gradient */
   current_iterate.objective_gradient.for_each([&](size_t i, double derivative) {
      rhs[i] += problem.objective_sign * derivative;
   });

   /* variable bound constraints */
   for (size_t i = 0; i < current_iterate.x.size(); i++) {
      rhs[i] -= current_iterate.multipliers.lower_bounds[i];
      rhs[i] -= current_iterate.multipliers.upper_bounds[i];
   }

   // solve the system
   solver.factorize(matrix);
   std::vector<double> solution = solver.solve(matrix, rhs);
   DEBUG << "Solution: ";
   print_vector(DEBUG, solution);

   // if least-square multipliers too big, discard them. Otherwise, store them
   const size_t number_variables = current_iterate.x.size();
   if (norm_inf(solution, current_iterate.x.size(), problem.number_constraints) <= multipliers_max_size) {
      for (size_t j = 0; j < problem.number_constraints; j++) {
         multipliers[j] = solution[number_variables + j];
      }
   }
}

#endif // SUBPROBLEM_H
