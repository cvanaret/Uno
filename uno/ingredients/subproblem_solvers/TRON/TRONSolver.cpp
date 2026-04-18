// Copyright (c) 2026 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <cmath>
#include <optional>
#include <stdexcept>
#include "TRONSolver.hpp"
#include "ingredients/subproblem/Subproblem.hpp"
#include "linear_algebra/BLAS.hpp"
#include "linear_algebra/Vector.hpp"

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

   TRONStats TRONSolver::solve(const HessianOperator& hessian_operator, const std::vector<double>& x0) {
      /* TODO
      const double atol = (atol > 0) ? atol : std::sqrt(std::numeric_limits<double>::epsilon());
      const double rtol = (rtol > 0) ? rtol : std::sqrt(std::numeric_limits<double>::epsilon());
      */

      TRONStats stats;

      // ── Initialise x within bounds ───────────────────────────────────
      x_ = x0;
      project_bounds(x_);

      auto [fx, gx_new] = fg_(x_);
      gx_ = gx_new;
      stats.objective = fx;

      // Projected-gradient norm at x0
      project_step(gpx_, x_, gx_, -1.0); // gpx = P(x - g) - x
      double pi0 = norm_2(gpx_);
      double eps_tol = atol + rtol * pi0;
      double fmin = std::min(-1.0, fx) / std::numeric_limits<double>::epsilon();

      stats.dual_residual = pi0;

      // Trust-region radius
      double radius = std::min(std::max(1.0, pi0 / 10.0), max_radius_);

      // Cauchy step length (persistent across iterations)
      double alpha_c = 1.0;
      int num_success = 0;



      // ── Main loop ────────────────────────────────────────────────────
      while (stats.status == SolveStatus::Unknown) {
         bool optimal = (pi0 <= eps_tol);
         bool unbounded = (fx < fmin);
         if (optimal) {
            stats.status = SolveStatus::Optimal;
            break;
         }
         if (unbounded) {
            stats.status = SolveStatus::Unbounded;
            break;
         }
         if (stats.iter >= max_iter) {
            stats.status = SolveStatus::MaxIter;
            break;
         }

         // Save current point
         xc_ = x_;
         double fc = fx;

         // ── Cauchy step ──────────────────────────────────────────────
         auto cauchy_status = cauchy_step(hessian_operator, alpha_c, radius);
         if (cauchy_status != SolveStatus::Unknown) {
            stats.status = cauchy_status;
            break;
         }

         // ── Projected Newton refinement (CG) ────────────────────────
         std::string cg_info = projected_newton(hessian_operator, radius);

         // ── Ratio test ───────────────────────────────────────────────
         double slope = dot(gx_, s_);
         double qs = dot(s_, Hs_) / 2.0 + slope;

         auto [fx_new, dummy_g] = fg_(x_);
         (void) dummy_g;
         fx = fx_new;

         double ared = fc - fx;
         double pred = -qs;
         if (pred >= 0.0) {
            stats.status = SolveStatus::NegPred;
            break;
         }

         double ratio = ared / pred;

         if (ratio >= eta1_) {
            // accepted
            num_success++;
            auto [fx2, gx2] = fg_(x_);
            gx_ = gx2;
            project_step(gpx_, x_, gx_, -1.0);
            pi0 = norm_2(gpx_);
         }
         else {
            // rejected
            fx = fc;
            x_ = xc_;
         }

         // ── Update trust-region radius ───────────────────────────────
         double s_norm = norm_2(s_);
         if (num_success == 0)
            radius = std::min(radius, s_norm);

         if (ratio < eta1_)
            radius = std::max(min_radius_, radius / 4.0);
         else if (ratio >= eta2_)
            radius = std::min(max_radius_, radius * 4.0);

         // ── Update stats ─────────────────────────────────────────────
         stats.iter++;
         stats.objective = fx;
         stats.dual_residual = pi0;

         if (verbose > 0 && stats.iter % verbose == 0) {
            std::printf("iter %4d  f=%-14.6e  π=%-10.3e  Δ=%-10.3e  %s\n", stats.iter, fx, pi0, radius, cg_info.c_str());
         }
      }
      if (stats.status == SolveStatus::Unknown)
         stats.status = SolveStatus::MaxIter;

      // TODO
      // stats.solution = x_;
      return stats;
   }

   // protected member functions

   // Project v component-wise onto [ℓ, u]
   void TRONSolver::project_bounds(Vector<double> &v) const {
      for (size_t i = 0; i < v.size(); ++i)
         v[i] = std::max(this->lower_bounds[i], std::min(v[i], this->upper_bounds[i]));
   }

   /// Active-set indicator: ifix[i] = true if x[i] at a bound AND gradient pushes into it
   void TRONSolver::active_set(std::vector<bool>& ifix, const Vector<double>& x) const {
      for (size_t i = 0; i < x.size(); ++i) {
         ifix[i] = (x[i] <= this->lower_bounds[i] || x[i] >= this->upper_bounds[i]);
      }
   }

   // s = P(x + alpha*d) - x
   void TRONSolver::project_step(Vector<double>& s, const Vector<double>& x, const Vector<double>& d, double alpha) const {
      for (size_t i = 0; i < x.size(); ++i) {
         // TODO dinstiguish the cases
         s[i] = std::max(this->lower_bounds[i], std::min(x[i] + alpha * d[i], this->upper_bounds[i])) - x[i];
      }
   }

   /// Hs = H*s, slope = gᵀs, qs = ½sᵀHs + gᵀs
   std::pair<double, double> TRONSolver::compute_Hs_slope_qs(const HessianOperator& hessian_operator, const Vector<double> &s,
         const Vector<double> &g) {
      hessian_operator(s, Hs_); // at xc_
      double slope = dot(g, s);
      double qs = dot(s, Hs_) / 2.0 + slope;
      return {slope, qs};
   }

   // This subroutine computes the number of break-points, and the minimal and maximal break-points of the projection of
   // x + alpha*w on the n-dimensional interval [xl,xu].
   BreakPoints TRONSolver::compute_break_points(const Vector<double>& x, const Vector<double>& w) const {
      BreakPoints break_points{0, 0., 0.};

      for (size_t i = 0; i < x.size(); ++i) {
         std::optional<double> break_point = std::nullopt;
         if (x[i] < this->upper_bounds[i] && w[i] > 0.) {
            break_point = (this->upper_bounds[i] - x[i]) / w[i];
         }
         else if (x[i] > this->lower_bounds[i] && w[i] < 0.) {
            break_point = (this->lower_bounds[i] - x[i]) / w[i];
         }
         if (break_point.has_value()) {
            break_points.number++;
            if (break_points.number == 1) {
               break_points.min = *break_point;
               break_points.max = *break_point;
            }
            else {
               break_points.min = std::min(*break_point, break_points.min);
               break_points.max = std::max(*break_point, break_points.max);
            }
         }
      }

      // handle the exceptional case.
      if (break_points.number == 0) {
         break_points.min = 0.;
         break_points.max = 0.;
      }
      return break_points;
   }

   /**
    * Backtracking projected line search: find smallest t = 2^{-k} s.t.  q(s) ≤ μ₀ gᵀs
    * where s = P(x + t*d) - x.
    * x is updated in-place.
    */
   void TRONSolver::projected_line_search(const HessianOperator& hessian_operator, Vector<double> &x, const Vector<double> &d,
         const Vector<double> &g) const {
      const BreakPoints break_points = compute_break_points(x, d);
      double alpha = 1.0;
      Vector<double> projected_step(x.size()), Hs(x.size());

      bool search = true;
      while (search && alpha > break_points.min) {
         project_step(projected_step, x, d, alpha);
         hessian_operator(projected_step, Hs);
         double slope = dot(g, projected_step);
         double qs = dot(projected_step, Hs) / 2. + slope;
         if (qs <= mu0 * slope) {
            search = false;
         }
         else {
            alpha /= 2.0;
         }
      }
      if (alpha < std::min(1., break_points.min)) {
         alpha = break_points.min;
         project_step(projected_step, x, d, alpha);
      }
      project_step(projected_step, x, d, alpha);
      x += projected_step;
      project_bounds(x);

      // Update Hs_ for the full step s_
      // (caller recomputes via hv_ after returning)
   }

   /**
       * Computes s = P(x - α g) - x satisfying sufficient decrease.
       * Updates x_ = xc_ (caller is responsible for xc_ being current).
       * Modifies s_, Hs_, alpha_c in-place.
       */
   SolveStatus TRONSolver::cauchy_step(const HessianOperator& hessian_operator, double& alpha, double radius) {
      // Negative gradient direction for breakpoints
      for (size_t i = 0; i < xc_.size(); ++i) {
         temp_[i] = -gx_[i];
      }
      const BreakPoints break_points = compute_break_points(xc_, temp_);

      std::fill(s_.begin(), s_.end(), 0.0);
      std::fill(Hs_.begin(), Hs_.end(), 0.0);

      project_step(s_, xc_, gx_, -alpha);
      double s_norm = norm_2(s_);

      bool interp;
      if (s_norm > mu1 * radius) {
         interp = true;
      }
      else {
         auto [slope, qs] = compute_Hs_slope_qs(hessian_operator, s_, gx_);
         interp = (qs >= mu0 * slope);
      }

      if (interp) {
         bool search = true;
         while (search) {
            alpha /= sigma;
            project_step(s_, xc_, gx_, -alpha);
            s_norm = norm_2(s_);
            if (s_norm <= mu1 * radius) {
               auto [slope, qs] = compute_Hs_slope_qs(hessian_operator, s_, gx_);
               search = (qs >= mu0 * slope);
            }
            if (alpha < std::sqrt(std::numeric_limits<double>::min()))
               return SolveStatus::SmallStep;
         }
      }
      else {
         double alpha_s = alpha;
         bool search = true;
         while (search && alpha <= break_points.max) {
            alpha *= sigma;
            project_step(s_, xc_, gx_, -alpha);
            s_norm = norm_2(s_);
            if (s_norm <= mu1 * radius) {
               auto [slope, qs] = compute_Hs_slope_qs(hessian_operator, s_, gx_);
               if (qs <= mu0 * slope)
                  alpha_s = alpha;
            } else {
               search = false;
            }
         }
         alpha = alpha_s;
         project_step(s_, xc_, gx_, -alpha);
      }

      // Apply Cauchy step to x_
      for (size_t i = 0; i < x_.size(); ++i) {
         x_[i] = xc_[i] + s_[i];
      }
      project_bounds(x_);
      return SolveStatus::Unknown;
   }

   // Scale d in-place so that ‖d + t*p‖ = Delta (solve quadratic for t ≥ 0)
   double TRONSolver::compute_distance_to_trust_region(const Vector<double>& d, const Vector<double>& p, double radius) {
      const double dd = dot(d, d);
      const double dp = dot(d, p);
      const double pp = dot(p, p);
      double discriminant = dp * dp - pp * (dd - radius * radius);
      if (discriminant < 0.0) {
         discriminant = 0.0;
      }
      return (-dp + std::sqrt(discriminant)) / pp;
      // axpy(t, p, d);
   }

   std::string TRONSolver::projected_newton(const HessianOperator& hessian_operator, double Delta) {
         const size_t n = this->x_.size();

         std::vector<bool> ifix(n, false);
         Vector<double> rhs(n), r(n), p(n), Hp(n), d(n);

         // Hessian-vector product with free-variable masking
         // ZHZ: (Hp)_i = ifix[i] ? 0 : (H * (ifix-masked d))_i
         auto masked_hv = [&](const Vector<double>& v, Vector<double>& Hv) {
            // zero out fixed components, then apply H
            for (size_t i = 0; i < n; ++i) {
               d[i] = ifix[i] ? 0.0 : v[i];
            }
            hessian_operator(d, Hv);
            for (size_t i = 0; i < n; ++i) {
               if (ifix[i]) {
                  Hv[i] = 0.0;
               }
            }
         };

         // Update Hs_ = H * s_
         hessian_operator(s_, Hs_);

         // x_ = xc_ + s_ projected
         for (size_t i = 0; i < n; ++i) {
            x_[i] = xc_[i] + s_[i];
         }
         project_bounds(x_);

         std::string exit_status = "maximum number of iterations";
         int iters = 0;

         bool exit_optimal = false, exit_pcg = false, exit_itmax = false;

         while (!(exit_optimal || exit_pcg || exit_itmax)) {
            active_set(ifix, x_);
            int n_free = 0;
            for (size_t i = 0; i < n; ++i) {
               if (!ifix[i]) {
                  n_free++;
               }
            }
            if (n_free == 0) {
               exit_optimal = true;
               continue;
            }

            // Build RHS = -(g + Hs) for free variables
            double gfnorm = 0.0;
            for (size_t i = 0; i < n; ++i) {
               rhs[i] = ifix[i] ? 0.0 : -(gx_[i] + Hs_[i]);
               gfnorm += rhs[i] * rhs[i];
            }
            double gfnorm_sqrt = std::sqrt(gfnorm);

            // ── Conjugate Gradient (Steihaug-Toint style) ────────────────
            d.fill(0.);
            r = rhs;
            p = r;
            double norm_r = norm_2(r);
            bool cg_done = false;
            std::string cg_flag = "maximum number of iterations";

            for (int cg_it = 0; cg_it < max_cgiter && !cg_done; ++cg_it) {
               masked_hv(p, Hp);
               double pHp = dot(p, Hp);
               if (pHp <= 0.0) {
                  // Negative curvature: go to boundary
                  const double t = compute_distance_to_trust_region(d, p, Delta);
                  d += t*p;
                  cg_flag = "on trust-region boundary";
                  cg_done = true;
                  break;
               }
               double alpha_cg = norm_r / pHp;

               // Trial step
               for (size_t i = 0; i < n; ++i) {
                  w_[i] = d[i] + alpha_cg * p[i];
               }
               if (norm_2(w_) >= Delta) {
                  const double t = compute_distance_to_trust_region(d, p, Delta);
                  d += t*p;
                  cg_flag = "on trust-region boundary";
                  cg_done = true;
                  break;
               }
               d = w_;

               // Update residual r = r - alpha * Hp
               // TODO
               // r -= alpha_cg*Hp;
               double rr_new = dot(r, r);

               if (std::sqrt(rr_new) <= cgtol * gfnorm_sqrt) {
                  cg_flag = "converged";
                  cg_done = true;
                  break;
               }
               double beta = rr_new / norm_r;
               // p = r + beta*p
               for (size_t i = 0; i < n; ++i) {
                  p[i] = r[i] + beta * p[i];
               }
               norm_r = rr_new;
            }

            // Projected line search along d from current x_
            // First negate rhs (= g + Hs) for the line-search gradient arg
            for (size_t i = 0; i < n; ++i) {
               rhs[i] = -rhs[i]; // now = (g+Hs) restricted
            }
            projected_line_search(hessian_operator, x_, d, rhs);

            // w_ = x_ - (xc_ + s_)  →  delta_s
            for (size_t i = 0; i < n; ++i) {
               w_[i] = x_[i] - xc_[i] - s_[i];
            }
            s_ += w_;

            hessian_operator(s_, Hs_);

            // Check optimality: ‖(g + Hs) restricted‖ ≤ cgtol * gfnorm_sqrt
            double newnorm = 0.0;
            for (size_t i = 0; i < n; ++i) {
               double ri = ifix[i] ? 0.0 : (gx_[i] + Hs_[i]);
               newnorm += ri * ri;
            }
            if (std::sqrt(newnorm) <= cgtol * gfnorm_sqrt)
               exit_optimal = true;
            else if (cg_flag == "on trust-region boundary")
               exit_pcg = true;

            iters++;
            if (iters >= max_cgiter) {
               exit_itmax = true;
            }
         }

         if (exit_optimal) return "stationary point found";
         if (exit_pcg) return "on trust-region boundary";
         if (exit_itmax) return "maximum number of iterations";
         return exit_status;
      }
} // namespace