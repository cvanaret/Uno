// Generator of large nonlinear (non-quadratic) test instances, to compare cppasl and ASL2 on the AD path.
// usage: generate_nlp_nl family n text|binary output.nl
// families:
//   chain     partially separable, small elements: min sum exp(x_i - x_{i+1}) + x_i^4 + sum log(1 + x_i^2)
//             s.t. sin(x_i) x_{i+1}^2 + cos(x_{i+2}) x_i = 0.5, i < n-2   (m = n-2)
//   defined   common expressions shared by functions: v_i = x_i x_{i+1} + sin(x_i) + x_{i+2};
//             min sum (v_i^2 + exp(-v_i)), s.t. v_i * v_{i+1} + atan(v_i) <= 1   (m = n/2)
//   blocks    medium elements: min sum_k log(sum_{j in block k} exp(x_j)) + sum x_j^4 (blocks of 20),
//             s.t. sqrt(1 + sum_{j in block k} (x_j - 0.1)^2) * cos(x_first) <= 2 for each block
//   shared    defined variables shared by many functions: v_k = log(sum_{j in block k} exp(x_j)) (blocks of 20),
//             min sum_k v_k^2, s.t. v_k * x_r + sin(v_k) <= 1 for 20 constraints per block (r in the next block)
//   dense     one element over all variables: min sqrt(1 + sum_i sin(x_i)^2) + sum_i x_i^4, s.t. tanh(x_i x_{i+1}) >= -1
#include <cstdio>
#include <cstdlib>
#include <string>
#include "nl_writer.hpp"

using namespace nlwriter;

int main(int argc, char** argv) {
   if (argc != 5) {
      std::fprintf(stderr, "usage: %s chain|defined|blocks|dense n text|binary output.nl\n", argv[0]);
      return 1;
   }
   const std::string family = argv[1];
   const int n = std::atoi(argv[2]);
   const Format format = (std::string(argv[3]) == "binary") ? Format::Binary : Format::Text;
   NlProblem p;
   p.number_variables = n;
   std::vector<Expression> terms;
   auto add_constraint = [&](const Expression& body, std::vector<std::pair<int, double>> linear, double lower, double upper) {
      p.constraints.push_back({body, std::move(linear)});
      p.constraint_lower.push_back(lower);
      p.constraint_upper.push_back(upper);
   };
   if (family == "chain") {
      for (int i = 0; i + 1 < n; ++i) terms.push_back(op(44, {var(i) - var(i + 1)}));
      for (int i = 0; i < n; ++i) terms.push_back(pow(var(i), num(4.)));
      for (int i = 0; i < n; ++i) terms.push_back(op(43, {num(1.) + pow(var(i), num(2.))}));
      for (int i = 0; i + 2 < n; ++i) {
         add_constraint(op(41, {var(i)}) * pow(var(i + 1), num(2.)) + op(46, {var(i + 2)}) * var(i), {}, 0.5, 0.5);
      }
   }
   else if (family == "defined") {
      const int k = n - 2;
      for (int i = 0; i < k; ++i) p.defined_variables.push_back({{{i + 2, 1.}}, var(i) * var(i + 1) + op(41, {var(i)})});
      for (int i = 0; i < k; ++i) terms.push_back(pow(var(n + i), num(2.)) + op(44, {neg(var(n + i))}));
      for (int i = 0; i + 1 < k; i += 2) add_constraint(var(n + i) * var(n + i + 1) + op(49, {var(n + i)}), {}, -INFINITY, 1.);
   }
   else if (family == "blocks") {
      const int block = 20;
      for (int start = 0; start + block <= n; start += block) {
         std::vector<Expression> exponentials, squares{num(1.)};
         for (int j = start; j < start + block; ++j) {
            exponentials.push_back(op(44, {var(j)}));
            squares.push_back(pow(var(j) - num(0.1), num(2.)));
         }
         terms.push_back(op(43, {op_list(54, exponentials)}));
         add_constraint(op(39, {op_list(54, squares)}) * op(46, {var(start)}), {}, -INFINITY, 2.);
      }
      for (int i = 0; i < n; ++i) terms.push_back(pow(var(i), num(4.)));
   }
   else if (family == "shared") {
      const int block = 20, number_blocks = n / block;
      for (int k = 0; k < number_blocks; ++k) {
         std::vector<Expression> exponentials;
         for (int j = k * block; j < (k + 1) * block; ++j) exponentials.push_back(op(44, {var(j)}));
         p.defined_variables.push_back({{}, op(43, {op_list(54, exponentials)})});
      }
      for (int k = 0; k < number_blocks; ++k) terms.push_back(pow(var(n + k), num(2.)));
      for (int k = 0; k < number_blocks; ++k) {
         for (int t = 0; t < block; ++t) {
            const int r = ((k + 1) % number_blocks) * block + t;
            add_constraint(var(n + k) * var(r) + op(41, {var(n + k)}), {}, -INFINITY, 1.);
         }
      }
   }
   else if (family == "dense") {
      std::vector<Expression> sines{num(1.)};
      for (int i = 0; i < n; ++i) sines.push_back(pow(op(41, {var(i)}), num(2.)));
      terms.push_back(op(39, {op_list(54, sines)}));
      for (int i = 0; i < n; ++i) terms.push_back(pow(var(i), num(4.)));
      for (int i = 0; i + 1 < n; ++i) add_constraint(op(37, {var(i) * var(i + 1)}), {}, -1., INFINITY);
   }
   else {
      std::fprintf(stderr, "unknown family %s\n", family.c_str());
      return 1;
   }
   p.objectives.push_back({op_list(54, terms), {}});
   p.objective_senses = {0};
   for (int j = 0; j < n; ++j) p.initial_point.emplace_back(j, 0.1 + 0.3 * ((j * 7919) % 11) / 11.);
   const std::string contents = p.write(format);
   write_file(argv[4], contents);
   std::printf("wrote %s: n = %d, m = %zu, %zu bytes\n", argv[4], n, p.constraints.size(), contents.size());
   return 0;
}
