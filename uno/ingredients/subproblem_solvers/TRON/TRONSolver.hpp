// Copyright (c) 2026 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_TRONSOLVER_H
#define UNO_TRONSOLVER_H

#include <cmath>
#include <functional>
#include <vector>
#include "../SubproblemSolver.hpp"
#include "../SolverWorkspace.hpp"
#include "linear_algebra/Vector.hpp"

namespace uno {
   class TRONSolverWorkspace: public SolverWorkspace {
   public:
      TRONSolverWorkspace() = default;

      [[nodiscard]] double compute_hessian_quadratic_form(const Subproblem& /*subproblem*/, const Vector<double>& /*vector*/) const override {
         return 0.;
      }

      Vector<double> objective_gradient;
   };

   enum class SolveStatus {
      Unknown,
      Optimal,
      Unbounded,
      MaxIter,
      SmallStep,
      NegPred,
      Error
   };

   struct TRONStats {
      SolveStatus status = SolveStatus::Unknown;
      int iter = 0;
      double objective = 0.0;
      double dual_residual = 0.0; ///< ‖P(x - g) - x‖ (projected gradient)
      std::vector<double> solution;
   };

   struct BreakPoints {
      size_t number; // number of break points
      double min; // minimal break-point
      double max; // maximal break-point
   };

   class TRONSolver: public SubproblemSolver {
   public:
      using HessianOperator = std::function<void(const Vector<double>& x, Vector<double>& Hv)>;

      TRONSolver() = default;
      ~TRONSolver() override = default;

      void initialize_memory(const Subproblem& subproblem) override;

      void solve(Statistics& statistics, const Subproblem& subproblem, double trust_region_radius, const Vector<double>& initial_point,
         Direction& direction, Evaluations& current_evaluations, const WarmstartInformation& warmstart_information) override;
      /**
       * @param x0   Initial guess (length n); need not be feasible.
       * @param opts Solver options.
       * @return     Execution statistics including the solution vector.
       */
      TRONStats solve(const HessianOperator& hessian_operator, const std::vector<double>& x0);

      [[nodiscard]] SolverWorkspace& get_workspace() override;

   protected:
      std::vector<double> lower_bounds;
      std::vector<double> upper_bounds;
      TRONSolverWorkspace workspace{};

      // Workspace vectors (all length n_)
      Vector<double> x_, xc_, gx_, gpx_, s_, Hs_, temp_, w_;

      // parameters
      double mu0 = 1.0 / 100.0; ///< sufficient-decrease parameter  ∈ (0, 0.5)
      double mu1 = 1.0; ///< trust-region scaling            ∈ (0, ∞)
      double sigma = 10.0; ///< step-size update factor         ∈ (1, ∞)
      // trust-region constants
      double eta1_ = 0.1; // acceptance threshold (ratio)
      double eta2_ = 0.75; // "very successful" threshold
      double min_radius_ = 1e-10;
      double max_radius_ = std::min(1.0 / std::sqrt(2.0 * std::numeric_limits<double>::epsilon()), 100.0);
      // options
      int max_iter = 100000;
      int max_cgiter = 50;
      double max_time = 30.0; ///< wall-clock seconds
      double atol = 0.0; ///< absolute gradient tolerance (set in solve)
      double rtol = 0.0; ///< relative gradient tolerance (set in solve)
      double cgtol = 0.1; ///< CG sub-problem tolerance
      int verbose = 0;

      std::function<std::pair<double, Vector<double>> (const Vector<double>& x)> fg_;

      void project_bounds(Vector<double> &v) const;
      void active_set(std::vector<bool>& ifix, const Vector<double>& x) const;
      void project_step(Vector<double>& s, const Vector<double>& x, const Vector<double>& d, double alpha) const;
      [[nodiscard]] std::pair<double, double> compute_Hs_slope_qs(const HessianOperator& hessian_operator,
         const Vector<double>& s, const Vector<double> &g);
      [[nodiscard]] BreakPoints compute_break_points(const Vector<double>& x, const Vector<double>& w) const;
      void projected_line_search(const HessianOperator& hessian_operator, Vector<double> &x, const Vector<double> &d,
         const Vector<double> &g) const;
      [[nodiscard]] SolveStatus cauchy_step(const HessianOperator& hessian_operator, double& alpha, double radius);
      [[nodiscard]] static double compute_distance_to_trust_region(const Vector<double>& d, const Vector<double>& p,
         double radius) ;
      std::string projected_newton(const HessianOperator& hessian_operator, double Delta);
   };
} // namespace

#endif // UNO_TRONSOLVER_H