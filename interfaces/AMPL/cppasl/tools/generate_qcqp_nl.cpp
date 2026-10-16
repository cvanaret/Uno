// Generator of instances of Mittelmann's qcqp.mod (random QCQP with a known KKT point) as .nl files.
// The model (AMPL):
//    minimize 0.5 * sum_{i,j} Q[i,j] x_i x_j + g^T x
//    subject to  A_e x (+ 0.5 x^T P_l x)  = b_e   (ml linear, mq quadratic equalities)
//                A_i x (+ 0.5 x^T P_l x) <= b_i   (pl linear, pq quadratic inequalities)
// Q = LQ LQ^T (sd = 1, convex) or the symmetrization of LQ (sd = 0, nonconvex); LQ is lower triangular with diagonal
// ~ U(0,1) and off-diagonal ~ U(-10,10) with probability sq. P_l = LP LP^T with density sp. Row i of A keeps its
// largest N(0,1) entry and the others with probability sq. x* ~ N(0,1), equality multipliers ~ N(0,1), inequality
// multipliers ~ U(0,10) with probability plf (linear) / pqf (quadratic), and g, b chosen so that the KKT conditions
// hold at x*. (AMPL's random generator is not reproduced: instances are statistically, not bitwise, equivalent.)
//
// Encoding follows AMPL: nonlinear part  o2 n0.5 o54 <count> {o2 o2 n<coef> v<i> v<j>}, linear parts in J/G.
//
// usage: generate_qcqp_nl n ml mq pl pq sd sq sp plf pqf seed text|binary output.nl
#include <cstdio>
#include <cstdlib>
#include <random>
#include <string>
#include <vector>
#include "nl_writer.hpp"

namespace {
   using Dense = std::vector<double>; // row-major n x n

   struct Generator {
      int n;
      std::mt19937_64 random;
      std::uniform_real_distribution<double> uniform01{0., 1.};
      std::normal_distribution<double> normal01{0., 1.};

      /// lower-triangular factor: sparse columns
      std::vector<std::vector<std::pair<int, double>>> lower_factor(double density) {
         std::vector<std::vector<std::pair<int, double>>> columns(static_cast<std::size_t>(this->n));
         for (int i = 0; i < this->n; ++i) {
            for (int j = 0; j < i; ++j) {
               if (this->uniform01(this->random) < density) {
                  columns[static_cast<std::size_t>(j)].emplace_back(i, -10. + 20. * this->uniform01(this->random));
               }
            }
            columns[static_cast<std::size_t>(i)].emplace_back(i, this->uniform01(this->random));
         }
         return columns;
      }

      /// L L^T (dense) or the symmetrization of L
      Dense quadratic_matrix(double density, bool convex) {
         const auto columns = this->lower_factor(density);
         Dense matrix(static_cast<std::size_t>(this->n) * static_cast<std::size_t>(this->n), 0.);
         for (const auto& column: columns) {
            if (convex) {
               for (const auto& [i, li]: column) {
                  for (const auto& [j, lj]: column) matrix[static_cast<std::size_t>(i) * static_cast<std::size_t>(this->n) + static_cast<std::size_t>(j)] += li * lj;
               }
            }
         }
         if (!convex) {
            for (std::size_t j = 0; j < columns.size(); ++j) {
               for (const auto& [i, lij]: columns[j]) {
                  matrix[static_cast<std::size_t>(i) * static_cast<std::size_t>(this->n) + j] = lij;
                  matrix[j * static_cast<std::size_t>(this->n) + static_cast<std::size_t>(i)] = lij;
               }
            }
         }
         return matrix;
      }

      std::vector<std::pair<int, double>> random_row(double density) {
         std::vector<double> dense(static_cast<std::size_t>(this->n));
         int largest = 0;
         for (int j = 0; j < this->n; ++j) {
            dense[static_cast<std::size_t>(j)] = this->normal01(this->random);
            if (std::fabs(dense[static_cast<std::size_t>(j)]) > std::fabs(dense[static_cast<std::size_t>(largest)])) largest = j;
         }
         std::vector<std::pair<int, double>> row;
         for (int j = 0; j < this->n; ++j) {
            if (j == largest || this->uniform01(this->random) < density) row.emplace_back(j, dense[static_cast<std::size_t>(j)]);
         }
         return row;
      }
   };

   /// o2 n0.5 o54 k {o2 o2 n<M_ij> v_i v_j}
   void write_quadratic_expression(nlwriter::NlStreamWriter& writer, const Dense& matrix, int n) {
      std::size_t count = 0;
      for (double value: matrix) count += (value != 0.);
      if (count == 0) { writer.number(0.); return; }
      writer.operation(2);
      writer.number(0.5);
      if (count >= 3) { writer.operation(54); writer.count(static_cast<int>(count)); }
      else if (count == 2) writer.operation(0);
      for (int i = 0; i < n; ++i) {
         for (int j = 0; j < n; ++j) {
            const double value = matrix[static_cast<std::size_t>(i) * static_cast<std::size_t>(n) + static_cast<std::size_t>(j)];
            if (value == 0.) continue;
            writer.operation(2);
            writer.operation(2);
            writer.number(value);
            writer.variable(i);
            writer.variable(j);
         }
      }
   }
} // namespace

int main(int argc, char** argv) {
   if (argc != 14) {
      std::fprintf(stderr, "usage: %s n ml mq pl pq sd sq sp plf pqf seed text|binary output.nl\n", argv[0]);
      return 1;
   }
   const int n = std::atoi(argv[1]), ml = std::atoi(argv[2]), mq = std::atoi(argv[3]), pl = std::atoi(argv[4]), pq = std::atoi(argv[5]);
   const bool convex = std::atoi(argv[6]) != 0;
   const double sq = std::atof(argv[7]), sp = std::atof(argv[8]), plf = std::atof(argv[9]), pqf = std::atof(argv[10]);
   const auto format = (std::string(argv[12]) == "binary") ? nlwriter::Format::Binary : nlwriter::Format::Text;
   Generator generator{n, std::mt19937_64(std::strtoull(argv[11], nullptr, 10))};
   const auto un = static_cast<std::size_t>(n);

   std::vector<double> xstar(un);
   for (double& value: xstar) value = generator.normal01(generator.random);
   const Dense Q = generator.quadratic_matrix(sq, convex);

   // constraints: equalities (linear, then quadratic), inequalities (linear, then quadratic)
   const int m = ml + mq + pl + pq;
   std::vector<std::vector<std::pair<int, double>>> rows(static_cast<std::size_t>(m));
   std::vector<Dense> P(static_cast<std::size_t>(m));
   std::vector<double> lower(static_cast<std::size_t>(m)), upper(static_cast<std::size_t>(m));
   // gradient of the Lagrangian at x*: g = -(Q x* + sum_i multiplier_i grad c_i(x*))
   std::vector<double> g(un, 0.);
   for (std::size_t i = 0; i < un; ++i) {
      double sum = 0.;
      for (std::size_t j = 0; j < un; ++j) sum += Q[i * un + j] * xstar[j];
      g[i] = -sum;
   }
   for (int i = 0; i < m; ++i) {
      const bool is_equality = i < ml + mq;
      const bool is_quadratic = (i >= ml && i < ml + mq) || i >= ml + mq + pl;
      auto& row = rows[static_cast<std::size_t>(i)];
      row = generator.random_row(sq);
      double value = 0.;
      std::vector<double> gradient(un, 0.);
      for (const auto& [j, a]: row) { value += a * xstar[static_cast<std::size_t>(j)]; gradient[static_cast<std::size_t>(j)] += a; }
      if (is_quadratic) {
         P[static_cast<std::size_t>(i)] = generator.quadratic_matrix(sp, true);
         const Dense& Pi = P[static_cast<std::size_t>(i)];
         for (std::size_t r = 0; r < un; ++r) {
            double product = 0.;
            for (std::size_t c = 0; c < un; ++c) product += Pi[r * un + c] * xstar[c];
            gradient[r] += product;
            value += 0.5 * xstar[r] * product;
         }
      }
      double multiplier = 0.;
      if (is_equality) {
         multiplier = generator.normal01(generator.random);
         lower[static_cast<std::size_t>(i)] = upper[static_cast<std::size_t>(i)] = value;
      }
      else {
         const double active_probability = is_quadratic ? pqf : plf;
         const bool active = generator.uniform01(generator.random) < active_probability;
         multiplier = active ? 10. * generator.uniform01(generator.random) : 0.;
         lower[static_cast<std::size_t>(i)] = -INFINITY;
         upper[static_cast<std::size_t>(i)] = active ? value : value + generator.uniform01(generator.random);
      }
      if (multiplier != 0.) for (std::size_t j = 0; j < un; ++j) g[j] -= multiplier * gradient[j];
      if (is_quadratic) { // the linear part of a quadratic constraint lists all its variables
         std::vector<double> dense(un, 0.);
         std::vector<char> present(un, 0);
         for (const auto& [j, a]: row) { dense[static_cast<std::size_t>(j)] += a; present[static_cast<std::size_t>(j)] = 1; }
         const Dense& Pi = P[static_cast<std::size_t>(i)];
         for (std::size_t r = 0; r < un; ++r) for (std::size_t c = 0; c < un; ++c) if (Pi[r * un + c] != 0.) present[r] = present[c] = 1;
         row.clear();
         for (std::size_t j = 0; j < un; ++j) if (present[j]) row.emplace_back(static_cast<int>(j), dense[j]);
      }
   }

   // write the .nl file; AMPL (and ASL) require the nonlinear constraints to come first
   std::vector<int> order;
   for (int i = 0; i < m; ++i) if (!P[static_cast<std::size_t>(i)].empty()) order.push_back(i);
   for (int i = 0; i < m; ++i) if (P[static_cast<std::size_t>(i)].empty()) order.push_back(i);
   std::size_t jacobian_nonzeros = 0;
   std::vector<long long> column_counts(un, 0);
   int nonlinear_constraints = 0;
   for (int i = 0; i < m; ++i) {
      jacobian_nonzeros += rows[static_cast<std::size_t>(i)].size();
      for (const auto& entry: rows[static_cast<std::size_t>(i)]) ++column_counts[static_cast<std::size_t>(entry.first)];
      if (!P[static_cast<std::size_t>(i)].empty()) ++nonlinear_constraints;
   }
   int equations = ml + mq;
   nlwriter::NlStreamWriter writer(format);
   auto line = [&](const std::string& text) { writer.header_line(" " + text); };
   writer.header_line(std::string(format == nlwriter::Format::Text ? "g" : "b") + "3 1 1 0\t# problem qcqp");
   line(std::to_string(n) + " " + std::to_string(m) + " 1 0 " + std::to_string(equations) + " 0");
   line(std::to_string(nonlinear_constraints) + " 1");
   line("0 0");
   line(std::to_string(nonlinear_constraints > 0 ? n : 0) + " " + std::to_string(n) + " " + std::to_string(nonlinear_constraints > 0 ? n : 0));
   line("0 0 0 1");
   line("0 0 0 0 0");
   line(std::to_string(jacobian_nonzeros) + " " + std::to_string(n));
   line("0 0");
   line("0 0 0 0 0");
   for (int k = 0; k < m; ++k) {
      const int i = order[static_cast<std::size_t>(k)];
      writer.segment('C', {k});
      if (P[static_cast<std::size_t>(i)].empty()) writer.number(0.);
      else write_quadratic_expression(writer, P[static_cast<std::size_t>(i)], n);
   }
   writer.segment('O', {0, 0});
   write_quadratic_expression(writer, Q, n);
   writer.letter('r');
   writer.end_of_segment_header();
   for (int i: order) writer.bound(lower[static_cast<std::size_t>(i)], upper[static_cast<std::size_t>(i)]);
   writer.letter('b');
   writer.end_of_segment_header();
   for (int j = 0; j < n; ++j) writer.bound(-INFINITY, INFINITY);
   writer.segment('k', {n - 1});
   long long cumulative = 0;
   for (int j = 0; j + 1 < n; ++j) writer.count(static_cast<int>(cumulative += column_counts[static_cast<std::size_t>(j)]));
   for (int k = 0; k < m; ++k) {
      const int i = order[static_cast<std::size_t>(k)];
      writer.segment('J', {k, static_cast<long long>(rows[static_cast<std::size_t>(i)].size())});
      for (const auto& [j, a]: rows[static_cast<std::size_t>(i)]) writer.index_value(j, a);
   }
   writer.segment('G', {0, n});
   for (int j = 0; j < n; ++j) writer.index_value(j, g[static_cast<std::size_t>(j)]);
   nlwriter::write_file(argv[13], writer.buffer());
   std::printf("wrote %s: n = %d, m = %d, %zu bytes\n", argv[13], n, m, writer.buffer().size());
   return 0;
}
