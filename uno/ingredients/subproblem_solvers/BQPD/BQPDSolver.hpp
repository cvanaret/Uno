// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_BQPDSOLVER_H
#define UNO_BQPDSOLVER_H

#include <array>
#include <vector>
#include "ingredients/subproblem_solvers/QPSolver.hpp"
#include "BQPDEvaluationSpace.hpp"
#include "ingredients/subproblem_solvers/SubproblemStatus.hpp"

namespace uno {
   // forward declarations
   class Multipliers;
   class Options;
   class Subproblem;

   // see bqpd.f
   enum class BQPDStatus {
      OPTIMAL = 0,
      UNBOUNDED_PROBLEM = 1,
      BOUND_INCONSISTENCY = 2,
      INFEASIBLE = 3,
      INCORRECT_PARAMETER = 4,
      LP_INSUFFICIENT_SPACE = 5,
      REDUCED_HESSIAN_INSUFFICIENT_SPACE = 6,
      SPARSE_INSUFFICIENT_SPACE = 7,
      MAX_RESTARTS_REACHED = 8
   };

   enum BQPDMode {
      COLD_START = 0,
      ACTIVE_SET_EQUALITIES = 1, // cold start
      USER_DEFINED_ACTIVE_SET = 2, // hot start
      UNCHANGED_ACTIVE_SET = 3,
      UNCHANGED_ACTIVE_SET_AND_JACOBIAN = 4,
      UNCHANGED_ACTIVE_SET_AND_REDUCED_HESSIAN = 5,
      UNCHANGED_ACTIVE_SET_AND_JACOBIAN_AND_REDUCED_HESSIAN = 6, // warm start
   };

   class BQPDSolver : public QPSolver {
   public:
      explicit BQPDSolver(const Options& options);

      void initialize_memory(const Subproblem& subproblem) override;

      void solve(Statistics& statistics, Subproblem& subproblem, double trust_region_radius, const Vector<double>& initial_point,
         Direction& direction, const WarmstartInformation& warmstart_information) override;

      [[nodiscard]] EvaluationSpace& get_evaluation_space() override;

   private:
      BQPDEvaluationSpace evaluation_space;
      std::vector<double> lower_bounds{}, upper_bounds{}; // lower and upper bounds of variables and constraints

      int kmax{0};
      int (*pick_kmax)(size_t number_variables, size_t number_constraints);
      int mlp{1000};
      const size_t nprof{2000000};
      std::array<int, 100> info{};
      std::vector<double> alp{};
      std::vector<int> lp{}, active_set{};
      std::vector<double> w{}, gradient_solution{}, residuals{}, e{};
      size_t mxws{};
      size_t mxlws{};
      std::vector<double> ws{};
      std::vector<int> lws{};
      int k{0};
      int iprint{0}, nout{6};
      double fmin{-1e20};
      int peq_solution{0}, ifail{0};

      const bool print_subproblem;

      void set_up_subproblem(Statistics& statistics, const Subproblem& subproblem, double trust_region_radius,
         const WarmstartInformation& warmstart_information);
      void display_subproblem(const Subproblem& subproblem, const Vector<double>& initial_point) const;
      void solve_subproblem(const Subproblem& subproblem, const Vector<double>& initial_point, Direction& direction,
         const WarmstartInformation& warmstart_information);
      [[nodiscard]] static BQPDMode determine_mode(const WarmstartInformation& warmstart_information);
      void hide_pointers_in_workspace(Statistics& statistics, const Subproblem& subproblem);
      void set_multipliers(size_t number_variables, Multipliers& direction_multipliers) const;
      [[nodiscard]] static BQPDStatus bqpd_status_from_int(int ifail);
      [[nodiscard]] bool check_sufficient_workspace_size(BQPDStatus bqpd_status);
      [[nodiscard]] static SubproblemStatus status_from_bqpd_status(BQPDStatus bqpd_status);
   };

   // hide a pointer to an arbitrary object at a given position of lws (BQPD integer workspace)
   template <typename T>
   void hide_pointer(size_t position_in_lws, int lws[], const T& object) {
      intptr_t pointer_to_object = reinterpret_cast<intptr_t>(&object);
      std::copy(reinterpret_cast<const char *>(&pointer_to_object),
         reinterpret_cast<const char *>(&pointer_to_object) + sizeof(intptr_t),
         reinterpret_cast<char *>(lws) + position_in_lws*sizeof(intptr_t));
   }

   // retrieve a pointer to an arbitrary object at a given position of lws (BQPD integer workspace)
   template <typename T>
   T* retrieve_pointer(size_t position_in_lws, const int lws[]) {
      intptr_t pointer_to_object;
      std::copy(reinterpret_cast<const char *>(lws) + position_in_lws*sizeof(intptr_t),
         reinterpret_cast<const char *>(lws) + (position_in_lws+1)*sizeof(intptr_t),
         reinterpret_cast<char *>(&pointer_to_object));
      T* object = reinterpret_cast<T*>(pointer_to_object);
      return object;
   }
} // namespace

#endif // UNO_BQPDSOLVER_H