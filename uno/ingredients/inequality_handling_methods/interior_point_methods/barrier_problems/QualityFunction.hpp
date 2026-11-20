// Copyright (c) 2026 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_QUALITYFUNCTION_H
#define UNO_QUALITYFUNCTION_H

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <utility>
#include <vector>

// Linear quality function of Nocedal, Waechter & Waltz ("Adaptive barrier update strategies for nonlinear interior
// methods", SIOPT 2009), in the 2-norm-squared variant of IPOPT's QualityFunctionMuOracle (no centrality term):
//   q(sigma) = Dn (1 - alpha_d)^2 + Pn (1 - alpha_p)^2 + Cn sum_k [(c_k + alpha_p dx_k) (e_k + alpha_d dz_k)]^2
// along the step d(sigma) = d_aff + sigma d_cen, where alpha_p and alpha_d are the fraction-to-boundary steplengths.
// The reciprocal steplengths 1/alpha_p and 1/alpha_d are upper envelopes of affine functions of sigma, hence piecewise
// affine. On each interval where both are affine, q = P(sigma) / (L_p L_d)^2 with deg P <= 4: its stationary points are the
// real roots of a polynomial of degree <= 5. minimize() computes the global minimizer of q on [sigma_lo, sigma_up].
namespace uno::quality_function {
   // primal distance to the bound c (signed: negative for an upper bound) and bound dual e (same sign as c), with their
   // affine (dx_aff, dz_aff) and centering (dx_cen, dz_cen) directions
   struct ComplementarityPair {
      double distance, affine_primal, centering_primal;
      double dual, affine_dual, centering_dual;
   };

   struct Minimizer {
      double sigma;
      double quality;
      size_t number_intervals; // intervals on which both reciprocal steplengths are affine
   };

   // polynomials with coefficients in increasing degree
   using Polynomial = std::vector<double>;

   inline Polynomial add(const Polynomial& p, const Polynomial& q) {
      Polynomial result(std::max(p.size(), q.size()), 0.);
      for (size_t index = 0; index < p.size(); ++index) result[index] += p[index];
      for (size_t index = 0; index < q.size(); ++index) result[index] += q[index];
      return result;
   }

   inline Polynomial multiply(const Polynomial& p, const Polynomial& q) {
      if (p.empty() || q.empty()) return {};
      Polynomial result(p.size() + q.size() - 1, 0.);
      for (size_t i = 0; i < p.size(); ++i) {
         for (size_t j = 0; j < q.size(); ++j) {
            result[i + j] += p[i] * q[j];
         }
      }
      return result;
   }

   inline Polynomial scale(Polynomial p, double factor) {
      for (double& coefficient: p) coefficient *= factor;
      return p;
   }

   inline Polynomial derivative(const Polynomial& p) {
      if (p.size() <= 1) return {0.};
      Polynomial result(p.size() - 1);
      for (size_t index = 1; index < p.size(); ++index) result[index - 1] = static_cast<double>(index) * p[index];
      return result;
   }

   inline double evaluate(const Polynomial& p, double x) {
      double value = 0.;
      for (size_t index = p.size(); index-- > 0;) value = value * x + p[index];
      return value;
   }

   // real roots of p in the open interval (lo, hi): isolation by the critical points (roots of p'), then bisection
   inline void roots_in_interval(Polynomial p, double lo, double hi, std::vector<double>& roots) {
      while (p.size() > 1 && p.back() == 0.) p.pop_back();
      if (p.size() <= 1 || !(lo < hi)) return;
      if (p.size() == 2) {
         const double root = -p[0] / p[1];
         if (lo < root && root < hi) roots.push_back(root);
         return;
      }
      std::vector<double> nodes{lo};
      roots_in_interval(derivative(p), lo, hi, nodes);
      std::sort(nodes.begin() + 1, nodes.end());
      nodes.push_back(hi);
      for (size_t k = 0; k + 1 < nodes.size(); ++k) {
         double a = nodes[k], b = nodes[k + 1];
         double pa = evaluate(p, a);
         const double pb = evaluate(p, b);
         if (pa == 0.) {
            if (0 < k) roots.push_back(a); // interior node
            continue;
         }
         if ((pa < 0.) == (pb < 0.)) continue; // p is monotone on [a, b]: no sign change, no root
         while (4. * std::numeric_limits<double>::epsilon() * std::max(std::abs(a), std::abs(b)) < b - a) {
            const double midpoint = 0.5 * (a + b);
            if (midpoint <= a || b <= midpoint) break;
            const double pm = evaluate(p, midpoint);
            if (pm == 0.) {
               a = b = midpoint;
               break;
            }
            if ((pm < 0.) == (pa < 0.)) {
               a = midpoint;
               pa = pm;
            }
            else {
               b = midpoint;
            }
         }
         roots.push_back(0.5 * (a + b));
      }
   }

   // y = intercept + slope * sigma
   struct Line {
      double intercept;
      double slope;
   };

   // affine piece [start, end] of an upper envelope
   struct Piece {
      double start;
      double end;
      Line line;
   };

   // upper envelope of the lines and of the constant line y = 1 on [lo, hi] (convex hull trick)
   inline std::vector<Piece> upper_envelope(std::vector<Line> lines, double lo, double hi) {
      lines.push_back({1., 0.});
      std::sort(lines.begin(), lines.end(), [](const Line& line1, const Line& line2) {
         return (line1.slope < line2.slope) || (line1.slope == line2.slope && line1.intercept < line2.intercept);
      });
      const auto intersection = [](const Line& line1, const Line& line2) { // precondition: line1.slope < line2.slope
         return (line1.intercept - line2.intercept) / (line2.slope - line1.slope);
      };
      std::vector<Line> hull;
      for (size_t index = 0; index < lines.size(); ++index) {
         // keep the largest intercept per slope
         if (index + 1 < lines.size() && lines[index + 1].slope == lines[index].slope) continue;
         const Line& line = lines[index];
         while (2 <= hull.size() && intersection(hull[hull.size() - 2], line) <= intersection(hull[hull.size() - 2], hull.back())) {
            hull.pop_back();
         }
         hull.push_back(line);
      }
      std::vector<Piece> pieces;
      double start = lo;
      for (size_t index = 0; index < hull.size() && start < hi; ++index) {
         const double end = (index + 1 < hull.size()) ? intersection(hull[index], hull[index + 1]) : hi;
         if (start < end) {
            pieces.push_back({start, std::min(end, hi), hull[index]});
            start = std::min(end, hi);
         }
      }
      return pieces;
   }

   class QualityFunction {
   public:
      // dual_weight Dn, primal_weight Pn, complementarity_weight Cn
      QualityFunction(std::vector<ComplementarityPair> pairs, double tau, double dual_weight, double primal_weight,
            double complementarity_weight):
            pairs(std::move(pairs)), tau(tau), dual_weight(dual_weight), primal_weight(primal_weight),
            complementarity_weight(complementarity_weight) {
      }

      // fraction-to-boundary steplengths (primal, dual) of the step d_aff + sigma d_cen
      [[nodiscard]] std::pair<double, double> step_lengths(double sigma) const {
         double primal_step_length = 1., dual_step_length = 1.;
         for (const ComplementarityPair& pair: this->pairs) {
            const double primal_step = pair.affine_primal + sigma * pair.centering_primal;
            const double dual_step = pair.affine_dual + sigma * pair.centering_dual;
            if (pair.distance * primal_step < 0.) {
               primal_step_length = std::min(primal_step_length, -this->tau * pair.distance / primal_step);
            }
            if (pair.dual * dual_step < 0.) {
               dual_step_length = std::min(dual_step_length, -this->tau * pair.dual / dual_step);
            }
         }
         return {primal_step_length, dual_step_length};
      }

      [[nodiscard]] double operator()(double sigma) const {
         const auto [primal_step_length, dual_step_length] = this->step_lengths(sigma);
         double complementarity = 0.;
         for (const ComplementarityPair& pair: this->pairs) {
            const double product = (pair.distance + primal_step_length * (pair.affine_primal + sigma * pair.centering_primal)) *
               (pair.dual + dual_step_length * (pair.affine_dual + sigma * pair.centering_dual));
            complementarity += product * product;
         }
         return this->dual_weight * (1. - dual_step_length) * (1. - dual_step_length) +
            this->primal_weight * (1. - primal_step_length) * (1. - primal_step_length) +
            this->complementarity_weight * complementarity;
      }

      // global minimizer on [sigma_lo, sigma_up]. Cost: O(n log n) for the envelopes, O(n) for the complementarity sums,
      // O(1) per interval, then a few O(n) evaluations of q
      [[nodiscard]] Minimizer minimize(double sigma_lo, double sigma_up) const;

   private:
      const std::vector<ComplementarityPair> pairs;
      const double tau;
      const double dual_weight, primal_weight, complementarity_weight;

      // 1/alpha = max(1, max_k -(a_k + sigma b_k) / (tau c_k)), restricted to the pairs that may block
      [[nodiscard]] std::vector<Line> reciprocal_primal_step_lines() const;
      [[nodiscard]] std::vector<Line> reciprocal_dual_step_lines() const;
      // N_k M_k = A_k u w + B_k u + C_k w + D_k with u = L_p(sigma), w = L_d(sigma): the sums over k of the products of
      // A, B, C, D are polynomials in sigma independent of the interval
      [[nodiscard]] std::array<Polynomial, 10> complementarity_sums() const;
   };

   inline std::vector<Line> QualityFunction::reciprocal_primal_step_lines() const {
      std::vector<Line> lines;
      lines.reserve(this->pairs.size());
      for (const ComplementarityPair& pair: this->pairs) {
         if (pair.distance != 0.) {
            lines.push_back({-pair.affine_primal / (this->tau * pair.distance), -pair.centering_primal / (this->tau * pair.distance)});
         }
      }
      return lines;
   }

   inline std::vector<Line> QualityFunction::reciprocal_dual_step_lines() const {
      std::vector<Line> lines;
      lines.reserve(this->pairs.size());
      for (const ComplementarityPair& pair: this->pairs) {
         if (pair.dual != 0.) {
            lines.push_back({-pair.affine_dual / (this->tau * pair.dual), -pair.centering_dual / (this->tau * pair.dual)});
         }
      }
      return lines;
   }

   // sums AA BB CC DD AB AC AD BC BD CD
   inline std::array<Polynomial, 10> QualityFunction::complementarity_sums() const {
      std::array<Polynomial, 10> sums;
      for (Polynomial& sum: sums) sum.assign(5, 0.);
      for (const ComplementarityPair& pair: this->pairs) {
         const double c = pair.distance, a = pair.affine_primal, b = pair.centering_primal;
         const double e = pair.dual, g = pair.affine_dual, h = pair.centering_dual;
         const double A = c * e;
         const double B0 = c * g, B1 = c * h; // B = c (g + sigma h)
         const double C0 = e * a, C1 = e * b; // C = e (a + sigma b)
         const double D0 = a * g, D1 = a * h + b * g, D2 = b * h; // D = (a + sigma b) (g + sigma h)
         sums[0][0] += A * A;
         sums[1][0] += B0 * B0; sums[1][1] += 2. * B0 * B1; sums[1][2] += B1 * B1;
         sums[2][0] += C0 * C0; sums[2][1] += 2. * C0 * C1; sums[2][2] += C1 * C1;
         sums[3][0] += D0 * D0; sums[3][1] += 2. * D0 * D1; sums[3][2] += D1 * D1 + 2. * D0 * D2; sums[3][3] += 2. * D1 * D2;
         sums[3][4] += D2 * D2;
         sums[4][0] += A * B0; sums[4][1] += A * B1;
         sums[5][0] += A * C0; sums[5][1] += A * C1;
         sums[6][0] += A * D0; sums[6][1] += A * D1; sums[6][2] += A * D2;
         sums[7][0] += B0 * C0; sums[7][1] += B0 * C1 + B1 * C0; sums[7][2] += B1 * C1;
         sums[8][0] += B0 * D0; sums[8][1] += B0 * D1 + B1 * D0; sums[8][2] += B0 * D2 + B1 * D1; sums[8][3] += B1 * D2;
         sums[9][0] += C0 * D0; sums[9][1] += C0 * D1 + C1 * D0; sums[9][2] += C0 * D2 + C1 * D1; sums[9][3] += C1 * D2;
      }
      return sums;
   }

   inline Minimizer QualityFunction::minimize(double sigma_lo, double sigma_up) const {
      Minimizer best{sigma_lo, (*this)(sigma_lo), 0};
      if (!(sigma_lo < sigma_up)) {
         return best;
      }
      const auto consider = [&](double sigma) {
         const double value = (*this)(sigma);
         if (value < best.quality) {
            best.sigma = sigma;
            best.quality = value;
         }
      };
      consider(sigma_up);

      // the candidates are ranked with the rational form of their interval (O(1)). The best ones are then evaluated with q,
      // which absorbs the rounding of the breakpoints
      constexpr size_t number_verified_candidates = 8;
      std::vector<std::pair<double, double>> candidates; // max-heap of (rational form, sigma)
      const auto add_candidate = [&](double value, double sigma) {
         if (!std::isfinite(value)) return;
         if (candidates.size() < number_verified_candidates) {
            candidates.emplace_back(value, sigma);
            std::push_heap(candidates.begin(), candidates.end());
         }
         else if (value < candidates.front().first) {
            std::pop_heap(candidates.begin(), candidates.end());
            candidates.back() = {value, sigma};
            std::push_heap(candidates.begin(), candidates.end());
         }
      };

      const std::vector<Piece> primal_pieces = upper_envelope(this->reciprocal_primal_step_lines(), sigma_lo, sigma_up);
      const std::vector<Piece> dual_pieces = upper_envelope(this->reciprocal_dual_step_lines(), sigma_lo, sigma_up);
      const std::array<Polynomial, 10> sums = this->complementarity_sums();
      std::vector<double> roots;
      size_t primal_index = 0, dual_index = 0;
      double start = sigma_lo;
      while (primal_index < primal_pieces.size() && dual_index < dual_pieces.size()) {
         const double end = std::min(primal_pieces[primal_index].end, dual_pieces[dual_index].end);
         if (start < end) {
            ++best.number_intervals;
            const Line& Lp = primal_pieces[primal_index].line;
            const Line& Ld = dual_pieces[dual_index].line;
            const Polynomial u{Lp.intercept, Lp.slope}, w{Ld.intercept, Ld.slope};
            const Polynomial u_minus_1{Lp.intercept - 1., Lp.slope}, w_minus_1{Ld.intercept - 1., Ld.slope};
            const Polynomial uu = multiply(u, u), ww = multiply(w, w), uw = multiply(u, w);
            // complementarity term: sum_k (A_k u w + B_k u + C_k w + D_k)^2 (degree <= 4)
            Polynomial complementarity = add(add(multiply(sums[0], multiply(uu, ww)), multiply(sums[1], uu)), add(multiply(sums[2], ww), sums[3]));
            complementarity = add(complementarity, scale(add(add(multiply(sums[4], multiply(uu, w)), multiply(sums[5], multiply(u, ww))),
               add(multiply(add(sums[6], sums[7]), uw), add(multiply(sums[8], u), multiply(sums[9], w)))), 2.));
            complementarity.resize(5); // the higher coefficients are structural zeros
            // q = numerator / (L_p L_d)^2 with numerator = Dn (L_d - 1)^2 L_p^2 + Pn (L_p - 1)^2 L_d^2 + Cn complementarity
            const Polynomial numerator = add(add(scale(multiply(multiply(w_minus_1, w_minus_1), uu), this->dual_weight),
               scale(multiply(multiply(u_minus_1, u_minus_1), ww), this->primal_weight)), scale(complementarity, this->complementarity_weight));
            const auto rational_form = [&](double sigma) {
               const double denominator = evaluate(uw, sigma);
               return evaluate(numerator, sigma) / (denominator * denominator);
            };
            // stationary points: numerator' L_p L_d - 2 numerator (L_p L_d)' = 0
            const Polynomial stationarity = add(multiply(derivative(numerator), uw), scale(multiply(numerator, derivative(uw)), -2.));
            roots.clear();
            roots_in_interval(stationarity, start, end, roots);
            for (const double root: roots) {
               add_candidate(rational_form(root), root);
            }
            // the minimum is often at a kink (a steplength reaches 1), where q may be very steep on one side (e.g. a large
            // Dn (1 - alpha_d)^2): the computed breakpoint carries a rounding error, so also consider points just inside
            const double offset = 1e-10 * (end - start);
            if (start != sigma_lo) {
               add_candidate(rational_form(start), start);
            }
            add_candidate(rational_form(start + offset), start + offset);
            add_candidate(rational_form(end - offset), end - offset);
         }
         start = std::max(start, end);
         if (primal_pieces[primal_index].end <= end) ++primal_index;
         if (dual_pieces[dual_index].end <= end) ++dual_index;
      }
      for (const auto& candidate: candidates) {
         consider(candidate.second);
      }
      return best;
   }
} // namespace

#endif // UNO_QUALITYFUNCTION_H
