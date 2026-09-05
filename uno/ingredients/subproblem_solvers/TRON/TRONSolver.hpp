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
      MaxIter,
      SmallStep,
      NegPred,
      Error
   };

   enum class CGStatus {
      SUCCESS,
      ON_TR_BOUNDARY,
      MAX_ITERATIONS
   };

   struct TRONStats {
      SolveStatus status = SolveStatus::Unknown;
      int iter = 0;
      double objective = 0.0;
      double dual_residual = 0.0; ///< ‖P(x - g) - x‖ (projected gradient)
      std::vector<double> solution;
   };

   struct BreakPoints {
      double min; // minimal break-point
      double max; // maximal break-point
   };

   class TRONSolver: public SubproblemSolver {
   public:
      using ObjectiveOperator = std::function<double(const Vector<double>& x)>;
      using GradientOperator = std::function<void(const Vector<double>& x, Vector<double>& gradient)>;
      using MatrixOperator = std::function<void(const Vector<double>& x, Vector<double>& Hv)>;

      TRONSolver() = default;
      ~TRONSolver() override = default;

      void initialize_memory(const Subproblem& subproblem) override;

      void solve(Statistics& statistics, const Subproblem& subproblem, double trust_region_radius, const Vector<double>& initial_point,
         Direction& direction, Evaluations& current_evaluations, const WarmstartInformation& warmstart_information) override;
      /**
       * @param x   Initial guess (length n); need not be feasible.
       * @param opts Solver options.
       * @return     Execution statistics including the solution vector.
       */
      void solve(const ObjectiveOperator& objective_operator, const GradientOperator& gradient_operator,
         const MatrixOperator& hessian_operator, const Vector<double>& initial_point);

      [[nodiscard]] SolverWorkspace& get_workspace() override;

   protected:
      std::vector<double> lower_bounds;
      std::vector<double> upper_bounds;
      TRONSolverWorkspace workspace{};

      // Workspace vectors (all length n_)
      Vector<double> x, xc, x_copy, gpx, gx, s, Hs, temp_, w;
      Vector<double> d_masked;
      Vector<double> quadratic_gradient;
      Vector<double> d_cg; // CG point
      double fc;

      // parameters
      double mu0 = 1.0 / 100.0; ///< sufficient-decrease parameter  ∈ (0, 0.5)
      double mu1 = 1.0; ///< trust-region scaling            ∈ (0, ∞)
      double sigma = 10.0; ///< step-size update factor         ∈ (1, ∞)
      // trust-region constants
      double eta1_ = 0.1; // acceptance threshold (ratio)
      double eta2_ = 0.75; // "very successful" threshold
      double min_radius_ = 1e-10;
      double max_radius = std::min(1.0 / std::sqrt(2.0 * std::numeric_limits<double>::epsilon()), 100.0);
      // options
      size_t max_iter = 100000;
      size_t max_cgiter = 50;
      double cgtol = 0.1; ///< CG sub-problem tolerance
      const double atol = std::sqrt(std::numeric_limits<double>::epsilon());
      const double rtol = std::sqrt(std::numeric_limits<double>::epsilon());

      std::function<std::pair<double, Vector<double>> (const Vector<double>& x)> fg_;

      void project_onto_bounds(Vector<double>& v) const;
      void compute_active_set(std::vector<bool>& active, const Vector<double>& x) const;
      void project_step(const Vector<double>& x, double alpha, const Vector<double>& d, Vector<double>& s) const;
      [[nodiscard]] std::pair<double, double> compute_Hs_slope_qs(const MatrixOperator& hessian_operator,
         const Vector<double>& s, const Vector<double>& g);
      [[nodiscard]] BreakPoints compute_break_points(const Vector<double>& x, const Vector<double>& d) const;
      void projected_line_search(const MatrixOperator& hessian_operator, Vector<double>& x, const Vector<double>& d,
         const Vector<double>& g, Vector<double>& s);
      [[nodiscard]] bool compute_cauchy_step(const MatrixOperator& hessian_operator, const Vector<double>& g, double& alpha,
         double radius);
      [[nodiscard]] static double compute_distance_to_trust_region(const Vector<double>& d, const Vector<double>& p,
         double radius);
      [[nodiscard]] CGStatus CG(Vector<double>& d, const MatrixOperator& matrix_operator, const Vector<double>& rhs,
         double radius, double gfnorm_sqrt);
      std::string projected_newton(const MatrixOperator& hessian_operator, const Vector<double>& g, double radius);
   };
} // namespace

#endif // UNO_TRONSOLVER_H