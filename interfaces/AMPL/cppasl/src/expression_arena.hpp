#pragma once

#include <cstdint>
#include <vector>

namespace cppasl::detail {

   /// Operation codes of the .nl format (ASL's opcode.hd); only codes that can appear in an .nl file.
   namespace nl_opcode {
      constexpr int plus = 0, minus = 1, multiply = 2, divide = 3, remainder = 4, power = 5, less = 6;
      constexpr int minimum_list = 11, maximum_list = 12, floor = 13, ceil = 14, abs = 15, negate = 16;
      constexpr int logical_or = 20, logical_and = 21, less_than = 22, less_equal = 23, equal = 24;
      constexpr int greater_equal = 28, greater_than = 29, not_equal = 30, logical_not = 34, if_then_else = 35;
      constexpr int tanh = 37, tan = 38, sqrt = 39, sinh = 40, sin = 41, log10 = 42, log = 43, exp = 44;
      constexpr int cosh = 45, cos = 46, atanh = 47, atan2 = 48, atan = 49, asinh = 50, asin = 51, acosh = 52;
      constexpr int acos = 53, sum_list = 54, integer_divide = 55, precision = 56, round = 57, trunc = 58;
      constexpr int count = 59, number_of = 60, number_of_strings = 61, at_least = 62, at_most = 63;
      constexpr int piecewise_linear = 64, symbolic_if = 65, exactly = 66, not_at_least = 67, not_at_most = 68;
      constexpr int not_exactly = 69, and_list = 70, or_list = 71, implies = 72, iff = 73, all_different = 74;
      constexpr int some_same = 75, logistic = 79, sign_power = 80;
      constexpr int number_opcodes = 81;
   } // namespace nl_opcode

   /// How the operands of an opcode are laid out in the file
   enum class OperandLayout : std::uint8_t {
      Invalid,
      Unary,
      Binary,
      Ternary,          ///< if-then-else, implies
      CountedList,      ///< count on its own line, then the operands (min, max, count, alldiff, ...)
      SumList,          ///< same as CountedList (count >= 3 for sums)
      PiecewiseLinear   ///< count k, then 2k-1 numbers (slopes and breakpoints), then the argument
   };

   OperandLayout operand_layout(int opcode);

   using NodeIndex = std::uint32_t;

   /// Special opcodes for leaves
   constexpr int constant_opcode = -1;
   constexpr int variable_opcode = -2; ///< index >= number of variables denotes a defined variable

   struct ExpressionNode {
      int opcode;
      std::uint32_t argument_count;
      std::uint32_t first_argument; ///< index into ExpressionArena::arguments; variable index for variables
      double value;                 ///< constants only
   };

   struct LinearTerm {
      int variable;
      double coefficient;
   };

   /// A common expression "V" segment: sum of linear terms + nonlinear expression
   struct DefinedVariable {
      std::uint32_t linear_begin{0};
      std::uint32_t linear_end{0};
      NodeIndex expression{0};
      bool is_defined{false};
   };

   /// Flat storage for expression trees (DAG through defined variables). Function expressions are appended,
   /// analyzed, then truncated away: only the defined variables stay resident.
   struct ExpressionArena {
      std::vector<ExpressionNode> nodes;
      std::vector<NodeIndex> arguments;
      std::vector<DefinedVariable> defined_variables;
      std::vector<LinearTerm> defined_variable_linear_terms;

      struct Mark {
         std::size_t node_count;
         std::size_t argument_count;
      };

      [[nodiscard]] Mark mark() const { return {this->nodes.size(), this->arguments.size()}; }
      void truncate(const Mark& mark) {
         this->nodes.resize(mark.node_count);
         this->arguments.resize(mark.argument_count);
      }
      NodeIndex add_constant(double value) {
         this->nodes.push_back({constant_opcode, 0, 0, value});
         return static_cast<NodeIndex>(this->nodes.size() - 1);
      }
      NodeIndex add_variable(int index) {
         this->nodes.push_back({variable_opcode, 0, static_cast<std::uint32_t>(index), 0.});
         return static_cast<NodeIndex>(this->nodes.size() - 1);
      }
      NodeIndex add_operation(int opcode, const NodeIndex* operands, std::uint32_t count) {
         const auto first = static_cast<std::uint32_t>(this->arguments.size());
         this->arguments.insert(this->arguments.end(), operands, operands + count);
         this->nodes.push_back({opcode, count, first, 0.});
         return static_cast<NodeIndex>(this->nodes.size() - 1);
      }
      [[nodiscard]] NodeIndex operand(const ExpressionNode& node, std::uint32_t k) const {
         return this->arguments[node.first_argument + k];
      }
   };

} // namespace cppasl::detail
