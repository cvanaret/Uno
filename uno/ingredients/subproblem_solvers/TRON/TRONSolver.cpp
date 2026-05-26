// Copyright (c) 2026 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <cmath>
#include <iostream>
#include <limits>
#include <optional>
#include <stdexcept>
#include "TRONSolver.hpp"
#include "ingredients/subproblem/Subproblem.hpp"
#include "linear_algebra/BLAS.hpp"
#include "linear_algebra/Vector.hpp"
#include "symbolic/UnaryNegation.hpp"

namespace uno {
   void TRONSolver::initialize_memory(const Subproblem& subproblem) {
      this->lower_bounds.resize(subproblem.number_variables);
      this->upper_bounds.resize(subproblem.number_variables);
      this->workspace.objective_gradient.resize(subproblem.number_variables);
   }

   void TRONSolver::solve(Statistics& /*statistics*/, const Subproblem& subproblem, double trust_region_radius,
         const Vector<double>& /*initial_point*/, Direction& direction, Evaluations& current_evaluations,
         const WarmstartInformation& /*warmstart_information*/) {
      if (0 < subproblem.number_constraints) {
         throw std::runtime_error("TRONSolver cannot solve problems with general constraints");
      }
      std::cout << "TRON solver created\n";
      throw std::runtime_error("TRONSolver created");
   }

   SolverWorkspace& TRONSolver::get_workspace() {
      return this->workspace;
   }

   void TRONSolver::solve(const ObjectiveOperator& objective_operator, const GradientOperator& gradient_operator,
         const MatrixOperator& hessian_operator, const Vector<double>& initial_point) {
      SolveStatus status = SolveStatus::Unknown;

      // initialize x within bounds
      this->x = initial_point;
      project_onto_bounds(x);
      gradient_operator(this->x, gx);

      // projected-gradient norm at x0
      project_step(x, -1., gx, this->gpx); // projected_gradient = P(x - g) - x
      double pix = norm_2(this->gpx);
      double epsilon = atol + rtol * pix;

      // trust-region radius
      double radius = std::min(std::max(1., pix / 10.), this->max_radius);

      if (pix <= epsilon) {
         status = SolveStatus::Optimal;
         return;
      }

      // TODO change point at which Hessian operator is evaluated

      // compute Cauchy step and store it in s
      double alpha_c = 1.;
      const bool cauchy_success = compute_cauchy_step(hessian_operator, this->gx, alpha_c, radius);
      if (!cauchy_success) {
         status = SolveStatus::Error;
         return;
      }

      // projected Newton refinement (CG)
      std::string cg_info = projected_newton(hessian_operator, gx, radius);
   }

   // protected member functions

   // project v component-wise onto [ℓ, u]
   void TRONSolver::project_onto_bounds(Vector<double>& v) const {
      for (size_t i = 0; i < v.size(); ++i)
         v[i] = std::max(this->lower_bounds[i], std::min(v[i], this->upper_bounds[i]));
   }

   // active-set indicator: active[i] = true if x[i] at a bound
   void TRONSolver::compute_active_set(std::vector<bool>& active, const Vector<double>& x) const {
      for (size_t i = 0; i < x.size(); ++i) {
         const double delta = (-INF<double> < this->lower_bounds[i] && this->lower_bounds[i] < this->upper_bounds[i] &&
            this->upper_bounds[i] < INF<double>) ? std::min(atol, rtol * (this->upper_bounds[i] - this->lower_bounds[i])) : atol;
         if (x[i] == this->lower_bounds[i] && x[i] == this->upper_bounds[i]) {
            active[i] = true;
         }
         else if (x[i] <= this->lower_bounds[i] + delta) {
            active[i] = true;
         }
         else if (x[i] >= this->upper_bounds[i] - delta) {
            active[i] = true;
         }
         else {
            active[i] = false;
         }
      }
   }

   // s = P(x + alpha*d) - x
   void TRONSolver::project_step(const Vector<double>& x, double alpha, const Vector<double>& d, Vector<double>& s) const {
      for (size_t i = 0; i < x.size(); ++i) {
         if (x[i] + alpha * d[i] < this->lower_bounds[i]) {
            s[i] = this->lower_bounds[i] - x[i];
         }
         else if (this->upper_bounds[i] < x[i] + alpha * d[i]) {
            s[i] = this->upper_bounds[i] - x[i];
         }
         else {
            s[i] = alpha * d[i];
         }
      }
   }

   // Hs = H*s, slope = gᵀs, qs = ½sᵀHs + gᵀs
   std::pair<double, double> TRONSolver::compute_Hs_slope_qs(const MatrixOperator& hessian_operator, const Vector<double>& s,
         const Vector<double>& g) {
      hessian_operator(s, Hs); // at xc_
      double slope = dot(g, s);
      double qs = dot(s, Hs) / 2. + slope;
      return {slope, qs};
   }

   // compute the minimal and maximal break-points of the projection of x + alpha*d on the n-dimensional interval [xl,xu].
   BreakPoints TRONSolver::compute_break_points(const Vector<double>& x, const Vector<double>& d) const {
      BreakPoints break_points{INF<double>, 0.};

      for (size_t i = 0; i < x.size(); ++i) {
         std::optional<double> break_point = std::nullopt;
         if (x[i] < this->upper_bounds[i] && d[i] > 0.) {
            break_point = (this->upper_bounds[i] - x[i]) / d[i];
         }
         else if (x[i] > this->lower_bounds[i] && d[i] < 0.) {
            break_point = (this->lower_bounds[i] - x[i]) / d[i];
         }
         if (break_point.has_value()) {
            break_points.min = std::min(*break_point, break_points.min);
            break_points.max = std::max(*break_point, break_points.max);
         }
      }
      return break_points;
   }

   // backtracking projected line search: find smallest t = 2^{-k} s.t.  q(s) ≤ μ₀ gᵀs, where s = P(x + t*d) - x.
   void TRONSolver::projected_line_search(const MatrixOperator& hessian_operator, Vector<double>& x, const Vector<double>& d,
         const Vector<double>& g, Vector<double>& s) {
      double alpha = 1.;
      const BreakPoints break_points = compute_break_points(x, d);

      bool search = true;
      while (search && alpha > break_points.min) {
         project_step(x, alpha, d, s);
         const auto [slope, qs] = compute_Hs_slope_qs(hessian_operator, s, g);
         if (qs <= mu0 * slope) {
            search = false;
         }
         else {
            alpha /= 2.;
         }
      }
      if (alpha < std::min(1., break_points.min)) {
         alpha = break_points.min;
         project_step(x, alpha, d, s);
         hessian_operator(s, this->Hs);
      }
      project_step(x, alpha, d, s);
      x += s;
   }

   // compute Cauchy step s = P(x - α g) - x satisfying sufficient decrease.
   // returns true upon success, false upon failure
   bool TRONSolver::compute_cauchy_step(const MatrixOperator& hessian_operator, const Vector<double>& g, double& alpha,
         double radius) {
      // compute breakpoints along negative gradient direction
      s = -g;
      const BreakPoints break_points = compute_break_points(x, s);

      s.fill(0.);
      Hs.fill(0.);

      project_step(x, -alpha, g, s);

      // interpolate or extrapolate
      bool interpolate = true;
      if (norm_2(s) <= mu1 * radius) {
         auto [slope, qs] = compute_Hs_slope_qs(hessian_operator, s, g);
         interpolate = (qs >= mu0 * slope);
      }

      if (interpolate) {
         bool search = true;
         while (search) {
            alpha /= sigma;
            project_step(x, -alpha, g, s);
            if (norm_2(s) <= mu1 * radius) {
               auto [slope, qs] = compute_Hs_slope_qs(hessian_operator, s, g);
               search = (qs >= mu0 * slope);
            }
            // TODO: correctly assess why this fails
            if (alpha < std::sqrt(std::nextafter(0., 1.))) {
               return false;
            }
         }
      }
      else { // extrapolation
         double alpha_success = alpha;
         bool search = true;
         while (search && alpha <= break_points.max) {
            alpha *= sigma;
            project_step(x, -alpha, g, s);
            if (norm_2(s) <= mu1 * radius) {
               auto [slope, qs] = compute_Hs_slope_qs(hessian_operator, s, g);
               if (qs <= mu0 * slope) {
                  alpha_success = alpha;
               }
            }
            else {
               search = false;
            }
         }
         // recover the last successful step
         alpha = alpha_success;
         project_step(x, -alpha, g, s);
      }
      return true;
   }

   // find t ≥ 0 so that ‖d + t*p‖ = Delta (solve quadratic for t ≥ 0)
   double TRONSolver::compute_distance_to_trust_region(const Vector<double>& d, const Vector<double>& p, double radius) {
      const double dTd = dot(d, d);
      const double dTp = dot(d, p);
      const double pTp = dot(p, p);
      double discriminant = dTp * dTp - pTp * (dTd - radius * radius);
      if (discriminant < 0.0) {
         discriminant = 0.0;
      }
      return (-dTp + std::sqrt(discriminant)) / pTp;
   }

   // Conjugate Gradient (Steihaug-Toint style)
   CGStatus TRONSolver::CG(Vector<double>& d, const MatrixOperator& matrix_operator, const Vector<double>& rhs, double radius,
         double gfnorm_sqrt) {
      const size_t n = d.size();
      Vector<double> r(n), p(n), Hp(n);
      d.fill(0.);
      r = rhs;
      p = r;
      double norm_r = norm_2(r);

      for (size_t cg_it = 0; cg_it < max_cgiter; ++cg_it) {
         matrix_operator(p, Hp);
         double pHp = dot(p, Hp);
         if (pHp <= 0.0) {
            // Negative curvature: go to boundary
            const double t = compute_distance_to_trust_region(d, p, radius);
            d += t*p;
            return CGStatus::ON_TR_BOUNDARY;
         }
         double alpha_cg = norm_r / pHp;

         // trial step
         for (size_t i = 0; i < n; ++i) {
            w[i] = d[i] + alpha_cg * p[i];
         }
         if (norm_2(w) >= radius) {
            const double t = compute_distance_to_trust_region(d, p, radius);
            d += t*p;
            return CGStatus::ON_TR_BOUNDARY;
         }
         d = w;

         // Update residual r = r - alpha * Hp
         // TODO
         // r -= alpha_cg*Hp;
         double rr_new = dot(r, r);

         if (std::sqrt(rr_new) <= cgtol * gfnorm_sqrt) {
            return CGStatus::SUCCESS;
         }
         double beta = rr_new / norm_r;
         // p = r + beta*p
         for (size_t i = 0; i < n; ++i) {
            p[i] = r[i] + beta * p[i];
         }
         norm_r = rr_new;
      }
      return CGStatus::MAX_ITERATIONS;
   }

   std::string TRONSolver::projected_newton(const MatrixOperator& hessian_operator, const Vector<double>& g, double radius) {
      const size_t n = this->x.size();

      // Update Hs = H * s
      hessian_operator(s, Hs);

      // projected Newton step
      bool exit_optimal = false, exit_pcg = false, exit_itmax = false;
      std::string exit_status = "maximum number of iterations";
      size_t iters = 0;
      x += s;
      project_onto_bounds(x);
      std::vector<bool> active(n, false);
      while (!exit_optimal && !exit_pcg && !exit_itmax) {
         compute_active_set(active, x);
         // stop if all bounds are active
         if (std::all_of(active.begin(), active.end(), [](bool active) { return active; })) {
            exit_optimal = true;
            continue;
         }

         // build RHS = -(g + Hs) for free variables
         double gfnorm = 0.0;
         for (size_t i = 0; i < n; ++i) {
            this->quadratic_gradient[i] = active[i] ? 0. : -g[i];
            gfnorm += this->quadratic_gradient[i] * this->quadratic_gradient[i];
            this->quadratic_gradient[i] -= active[i] ? 0. : Hs[i];
         }
         double gfnorm_sqrt = std::sqrt(gfnorm);

         // define ZHZ, the Hessian-vector product with free-variable masking
         // (Hp)_i = ifix[i] ? 0 : (H * (ifix-masked d))_i
         const auto ZHZ = [&](const Vector<double>& d, Vector<double>& result) {
            // zero out fixed components of d, then apply H
            for (size_t i = 0; i < n; ++i) {
               this->d_masked[i] = active[i] ? 0. : d[i];
            }
            hessian_operator(this->d_masked, result);
            // zero out fixed components of the result
            for (size_t i = 0; i < n; ++i) {
               if (active[i]) {
                  result[i] = 0.;
               }
            }
         };

         const CGStatus cg_status = CG(this->d_cg, ZHZ, this->quadratic_gradient, radius, gfnorm_sqrt);
         iters++;

         // projected line search
         this->quadratic_gradient.scale(-1.);
         projected_line_search(ZHZ, x, d_cg, this->quadratic_gradient, w);
         s += w;
         hessian_operator(s, Hs);

         // Check optimality: ‖(g + Hs) restricted‖ ≤ cgtol * gfnorm_sqrt
         double new_norm = 0.0;
         for (size_t i = 0; i < n; ++i) {
            const double ri = active[i] ? 0. : Hs[i] + g[i];
            new_norm += ri * ri;
         }
         if (std::sqrt(new_norm) <= cgtol * gfnorm_sqrt) {
            exit_optimal = true;
         }
         else if (cg_status == CGStatus::ON_TR_BOUNDARY) {
            exit_pcg = true;
         }
         else if (iters >= max_cgiter) {
            exit_itmax = true;
         }
      }

      if (exit_optimal) return "stationary point found";
      if (exit_pcg) return "on trust-region boundary";
      if (exit_itmax) return "maximum number of iterations";
      return exit_status;
   }
} // namespace