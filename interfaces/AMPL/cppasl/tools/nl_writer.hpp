#pragma once
// Minimal writer of AMPL .nl files (text 'g' and binary 'b'), used to generate test and benchmark instances
// without an AMPL license. NlStreamWriter writes tokens; NlProblem is a convenient in-memory description.
#include <algorithm>
#include <charconv>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <fstream>
#include <initializer_list>
#include <set>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace nlwriter {

   enum class Format { Text, Binary };

   /// Streaming token writer: each call writes one .nl "item" (a line in text mode).
   class NlStreamWriter {
   public:
      explicit NlStreamWriter(Format format): format(format) {}

      std::string& buffer() { return this->output; }
      [[nodiscard]] Format encoding() const { return this->format; }

      /// raw header line (always text)
      void header_line(const std::string& line) { this->output += line; this->output += '\n'; }

      // expression tokens
      void number(double value) { this->letter('n'); this->put_double(value); this->end(); }
      void variable(int index) { this->letter('v'); this->put_int(index); this->end(); }
      void operation(int opcode) { this->letter('o'); this->put_int(opcode); this->end(); }
      void count(int value) { this->put_int(value); this->end(); }

      // segment headers
      void segment(char letter, std::initializer_list<long long> integers) {
         this->letter(letter);
         bool first = true;
         for (long long value: integers) {
            if (this->format == Format::Text && !first) this->output += ' ';
            this->put_int(value);
            first = false;
         }
         this->end();
      }
      /// "index value" line (J, G, x, d, V linear terms)
      void index_value(int index, double value) {
         this->put_int(index);
         if (this->format == Format::Text) this->output += ' ';
         this->put_double(value);
         this->end();
      }
      /// one bound line: type 0..4 (r/b segments)
      void bound(double lower, double upper) {
         // text: "<type> <values>", as AMPL writes it
         const bool has_lower = std::isfinite(lower), has_upper = std::isfinite(upper);
         if (has_lower && has_upper && lower == upper) { this->bound_type('4'); this->put_double(lower); }
         else if (has_lower && has_upper) {
            this->bound_type('0'); this->put_double(lower);
            if (this->format == Format::Text) this->output += ' ';
            this->put_double(upper);
         }
         else if (has_upper) { this->bound_type('1'); this->put_double(upper); }
         else if (has_lower) { this->bound_type('2'); this->put_double(lower); }
         else this->letter('3');
         this->end();
      }
      void end_of_segment_header() { this->end(); } // e.g. after 'r' / 'b'
      void letter(char c) { this->output += c; }
      void suffix_header(int kind, int count, const std::string& name) {
         this->letter('S');
         this->put_int(kind);
         if (this->format == Format::Text) this->output += ' ';
         this->put_int(count);
         if (this->format == Format::Text) { this->output += ' '; this->output += name; }
         else { this->put_raw(static_cast<std::int32_t>(name.size())); this->output += name; }
         this->end();
      }
      void index_integer(int index, int value) {
         this->put_int(index);
         if (this->format == Format::Text) this->output += ' ';
         this->put_int(value);
         this->end();
      }

   private:
      Format format;
      std::string output;

      void bound_type(char type) {
         this->output += type;
         if (this->format == Format::Text) this->output += ' ';
      }
      void end() { if (this->format == Format::Text) this->output += '\n'; }
      template <typename T>
      void put_raw(T value) {
         char bytes[sizeof(T)];
         std::memcpy(bytes, &value, sizeof(T));
         this->output.append(bytes, sizeof(T));
      }
      void put_int(long long value) {
         if (this->format == Format::Binary) { this->put_raw(static_cast<std::int32_t>(value)); return; }
         char text[24];
         const auto result = std::to_chars(text, text + sizeof(text), value);
         this->output.append(text, result.ptr);
      }
      void put_double(double value) {
         if (this->format == Format::Binary) { this->put_raw(value); return; }
         char text[32];
         const auto result = std::to_chars(text, text + sizeof(text), value); // shortest round-trip representation
         this->output.append(text, result.ptr);
      }
   };

   // ---------------------------------------------------------------------------------------------------------------
   // expression trees in prefix form
   // ---------------------------------------------------------------------------------------------------------------

   struct Token {
      enum Kind : std::uint8_t { Number, Variable, Operation, Count } kind;
      double value;
      int integer;
   };

   /// An expression in .nl prefix notation
   struct Expression {
      std::vector<Token> tokens;
      [[nodiscard]] bool empty() const { return this->tokens.empty(); }
   };

   inline Expression num(double value) { return {{{Token::Number, value, 0}}}; }
   inline Expression var(int index) { return {{{Token::Variable, 0., index}}}; }
   inline Expression op(int opcode, std::initializer_list<Expression> operands) {
      Expression result{{{Token::Operation, 0., opcode}}};
      for (const Expression& operand: operands) result.tokens.insert(result.tokens.end(), operand.tokens.begin(), operand.tokens.end());
      return result;
   }
   /// n-ary operation with an explicit count (sumlist 54, min 11, max 12, ...)
   inline Expression op_list(int opcode, const std::vector<Expression>& operands) {
      Expression result{{{Token::Operation, 0., opcode}, {Token::Count, 0., static_cast<int>(operands.size())}}};
      for (const Expression& operand: operands) result.tokens.insert(result.tokens.end(), operand.tokens.begin(), operand.tokens.end());
      return result;
   }
   /// piecewise-linear term: slopes s_0..s_{k-1}, breakpoints b_0..b_{k-2}
   inline Expression piecewise_linear(const std::vector<double>& slopes, const std::vector<double>& breakpoints, const Expression& argument) {
      Expression result{{{Token::Operation, 0., 64}, {Token::Count, 0., static_cast<int>(slopes.size())}}};
      for (std::size_t j = 0; j < slopes.size(); ++j) {
         result.tokens.push_back({Token::Number, slopes[j], 0});
         if (j < breakpoints.size()) result.tokens.push_back({Token::Number, breakpoints[j], 0});
      }
      result.tokens.insert(result.tokens.end(), argument.tokens.begin(), argument.tokens.end());
      return result;
   }
   inline Expression operator+(const Expression& a, const Expression& b) { return op(0, {a, b}); }
   inline Expression operator-(const Expression& a, const Expression& b) { return op(1, {a, b}); }
   inline Expression operator*(const Expression& a, const Expression& b) { return op(2, {a, b}); }
   inline Expression operator/(const Expression& a, const Expression& b) { return op(3, {a, b}); }
   inline Expression pow(const Expression& a, const Expression& b) { return op(5, {a, b}); }
   inline Expression neg(const Expression& a) { return op(16, {a}); }

   inline void write_expression(NlStreamWriter& writer, const Expression& expression) {
      if (expression.empty()) { writer.number(0.); return; }
      for (const Token& token: expression.tokens) {
         switch (token.kind) {
            case Token::Number: writer.number(token.value); break;
            case Token::Variable: writer.variable(token.integer); break;
            case Token::Operation: writer.operation(token.integer); break;
            case Token::Count: writer.count(token.integer); break;
         }
      }
   }

   // ---------------------------------------------------------------------------------------------------------------
   // whole problems
   // ---------------------------------------------------------------------------------------------------------------

   struct DefinedVariable {
      std::vector<std::pair<int, double>> linear_terms;
      Expression expression;
   };

   struct Function {
      Expression expression;                           ///< nonlinear part (empty: none)
      std::vector<std::pair<int, double>> linear_terms; ///< linear part (J or G entries)
   };

   /// An optimization problem written in AMPL's variable order convention: nonlinear variables first (all declared
   /// nonlinear in both constraints and objectives), then the binary and integer variables (assumed linear).
   struct NlProblem {
      int number_variables{0};
      int number_binary_variables{0};
      int number_integer_variables{0};
      std::vector<DefinedVariable> defined_variables; ///< indices number_variables + k
      std::vector<Function> objectives;
      std::vector<int> objective_senses; ///< 0 minimize, 1 maximize
      std::vector<Function> constraints;
      std::vector<double> constraint_lower, constraint_upper;
      std::vector<double> variable_lower, variable_upper;
      std::vector<std::pair<int, double>> initial_point, initial_dual;
      /// (kind, name, [(index, value)])
      std::vector<std::tuple<int, std::string, std::vector<std::pair<int, double>>>> suffixes;

      /// variables appearing in an expression (through the defined variables too)
      void collect_variables(const Expression& expression, std::set<int>& variables) const {
         for (const Token& token: expression.tokens) {
            if (token.kind != Token::Variable) continue;
            if (token.integer < this->number_variables) variables.insert(token.integer);
            else {
               const DefinedVariable& defined = this->defined_variables[static_cast<std::size_t>(token.integer - this->number_variables)];
               for (const auto& term: defined.linear_terms) variables.insert(term.first);
               this->collect_variables(defined.expression, variables);
            }
         }
      }

      [[nodiscard]] std::string write(Format format) const {
         const int n = this->number_variables;
         const auto m = static_cast<int>(this->constraints.size());
         const auto number_objectives = static_cast<int>(this->objectives.size());
         // complete the linear parts with the variables of the nonlinear parts (coefficient 0), as AMPL does
         auto gradient_entries = [&](const Function& function) {
            std::set<int> variables;
            this->collect_variables(function.expression, variables);
            std::vector<std::pair<int, double>> entries;
            for (const auto& term: function.linear_terms) {
               variables.erase(term.first);
            }
            std::vector<std::pair<int, double>> merged(function.linear_terms);
            for (int v: variables) merged.emplace_back(v, 0.);
            std::sort(merged.begin(), merged.end());
            for (const auto& term: merged) {
               if (!entries.empty() && entries.back().first == term.first) entries.back().second += term.second;
               else entries.push_back(term);
            }
            return entries;
         };
         std::vector<std::vector<std::pair<int, double>>> jacobian_rows, gradient_rows;
         std::size_t jacobian_nonzeros = 0, gradient_nonzeros = 0;
         std::vector<long long> column_counts(static_cast<std::size_t>(n), 0);
         int nonlinear_constraints = 0, nonlinear_objectives = 0;
         for (const Function& constraint: this->constraints) {
            jacobian_rows.push_back(gradient_entries(constraint));
            jacobian_nonzeros += jacobian_rows.back().size();
            for (const auto& entry: jacobian_rows.back()) ++column_counts[static_cast<std::size_t>(entry.first)];
            if (!constraint.expression.empty()) ++nonlinear_constraints;
         }
         for (const Function& objective: this->objectives) {
            gradient_rows.push_back(gradient_entries(objective));
            gradient_nonzeros += gradient_rows.back().size();
            if (!objective.expression.empty()) ++nonlinear_objectives;
         }
         int ranges = 0, equations = 0;
         for (int i = 0; i < m; ++i) {
            const double l = this->constraint_lower[static_cast<std::size_t>(i)], u = this->constraint_upper[static_cast<std::size_t>(i)];
            if (std::isfinite(l) && std::isfinite(u)) (l == u) ? ++equations : ++ranges;
         }
         const bool nonlinear = nonlinear_constraints + nonlinear_objectives > 0;
         const int nonlinear_variables = nonlinear ? n - this->number_binary_variables - this->number_integer_variables : 0;
         const auto defined = static_cast<int>(this->defined_variables.size());

         NlStreamWriter writer(format);
         auto line = [&](std::initializer_list<long long> values, const char* comment) {
            std::string text = " ";
            bool first = true;
            for (long long value: values) {
               if (!first) text += ' ';
               text += std::to_string(value);
               first = false;
            }
            writer.header_line(text + "\t# " + comment);
         };
         writer.header_line(std::string(format == Format::Text ? "g" : "b") + "3 1 1 0\t# problem cppasl_generated");
         line({n, m, number_objectives, ranges, equations, 0}, "vars, constraints, objectives, ranges, eqns, lcons");
         line({nonlinear_constraints, nonlinear_objectives}, "nonlinear constraints, objectives");
         line({0, 0}, "network constraints: nonlinear, linear");
         line({nonlinear_constraints > 0 ? nonlinear_variables : 0, nonlinear_objectives > 0 ? nonlinear_variables : 0,
            (nonlinear_constraints > 0 && nonlinear_objectives > 0) ? nonlinear_variables : 0}, "nonlinear vars in constraints, objectives, both");
         line({0, 0, 0, 1}, "linear network variables; functions; arith, flags");
         line({this->number_binary_variables, this->number_integer_variables, 0, 0, 0}, "discrete variables: binary, integer, nonlinear (b,c,o)");
         line({static_cast<long long>(jacobian_nonzeros), static_cast<long long>(gradient_nonzeros)}, "nonzeros in Jacobian, gradients");
         line({0, 0}, "max name lengths: constraints, variables");
         line({defined, 0, 0, 0, 0}, "common exprs: b,c,o,c1,o1");

         for (int k = 0; k < defined; ++k) {
            const DefinedVariable& variable = this->defined_variables[static_cast<std::size_t>(k)];
            writer.segment('V', {n + k, static_cast<long long>(variable.linear_terms.size()), 0});
            for (const auto& term: variable.linear_terms) writer.index_value(term.first, term.second);
            write_expression(writer, variable.expression);
         }
         for (int i = 0; i < m; ++i) {
            writer.segment('C', {i});
            write_expression(writer, this->constraints[static_cast<std::size_t>(i)].expression);
         }
         for (int k = 0; k < number_objectives; ++k) {
            writer.segment('O', {k, this->objective_senses.empty() ? 0 : this->objective_senses[static_cast<std::size_t>(k)]});
            write_expression(writer, this->objectives[static_cast<std::size_t>(k)].expression);
         }
         if (!this->initial_dual.empty()) {
            writer.segment('d', {static_cast<long long>(this->initial_dual.size())});
            for (const auto& entry: this->initial_dual) writer.index_value(entry.first, entry.second);
         }
         if (!this->initial_point.empty()) {
            writer.segment('x', {static_cast<long long>(this->initial_point.size())});
            for (const auto& entry: this->initial_point) writer.index_value(entry.first, entry.second);
         }
         if (m > 0) {
            writer.letter('r');
            writer.end_of_segment_header();
            for (int i = 0; i < m; ++i) writer.bound(this->constraint_lower[static_cast<std::size_t>(i)], this->constraint_upper[static_cast<std::size_t>(i)]);
         }
         writer.letter('b');
         writer.end_of_segment_header();
         for (int j = 0; j < n; ++j) {
            const double l = this->variable_lower.empty() ? -INFINITY : this->variable_lower[static_cast<std::size_t>(j)];
            const double u = this->variable_upper.empty() ? INFINITY : this->variable_upper[static_cast<std::size_t>(j)];
            writer.bound(l, u);
         }
         if (m > 0) {
            writer.segment('k', {n - 1});
            long long cumulative = 0;
            for (int j = 0; j + 1 < n; ++j) {
               cumulative += column_counts[static_cast<std::size_t>(j)];
               writer.count(static_cast<int>(cumulative));
            }
            for (int i = 0; i < m; ++i) {
               const auto& row = jacobian_rows[static_cast<std::size_t>(i)];
               if (row.empty()) continue;
               writer.segment('J', {i, static_cast<long long>(row.size())});
               for (const auto& entry: row) writer.index_value(entry.first, entry.second);
            }
         }
         for (int k = 0; k < number_objectives; ++k) {
            const auto& row = gradient_rows[static_cast<std::size_t>(k)];
            if (row.empty()) continue;
            writer.segment('G', {k, static_cast<long long>(row.size())});
            for (const auto& entry: row) writer.index_value(entry.first, entry.second);
         }
         for (const auto& [kind, name, values]: this->suffixes) {
            writer.suffix_header(kind, static_cast<int>(values.size()), name);
            for (const auto& entry: values) {
               if (kind & 4) writer.index_value(entry.first, entry.second);
               else writer.index_integer(entry.first, static_cast<int>(entry.second));
            }
         }
         return std::move(writer.buffer());
      }
   };

   inline void write_file(const std::string& file_name, const std::string& contents) {
      std::ofstream file(file_name, std::ios::binary);
      if (!file) throw std::runtime_error("cannot write " + file_name);
      file.write(contents.data(), static_cast<std::streamsize>(contents.size()));
   }

} // namespace nlwriter
