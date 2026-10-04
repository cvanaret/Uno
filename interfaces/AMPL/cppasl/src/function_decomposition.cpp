#include "function_decomposition.hpp"
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>

namespace cppasl::detail {

   OperandLayout operand_layout(int opcode) {
      namespace op = nl_opcode;
      switch (opcode) {
         case op::floor: case op::ceil: case op::abs: case op::negate: case op::logical_not: case op::tanh: case op::tan:
         case op::sqrt: case op::sinh: case op::sin: case op::log10: case op::log: case op::exp: case op::cosh:
         case op::cos: case op::atanh: case op::atan: case op::asinh: case op::asin: case op::acosh: case op::acos:
         case op::logistic:
            return OperandLayout::Unary;
         case op::plus: case op::minus: case op::multiply: case op::divide: case op::remainder: case op::power:
         case op::less: case op::logical_or: case op::logical_and: case op::less_than: case op::less_equal:
         case op::equal: case op::greater_equal: case op::greater_than: case op::not_equal: case op::atan2:
         case op::integer_divide: case op::precision: case op::round: case op::trunc: case op::at_least:
         case op::at_most: case op::exactly: case op::not_at_least: case op::not_at_most: case op::not_exactly:
         case op::iff: case op::sign_power:
            return OperandLayout::Binary;
         case op::if_then_else: case op::symbolic_if: case op::implies:
            return OperandLayout::Ternary;
         case op::minimum_list: case op::maximum_list: case op::count: case op::number_of: case op::number_of_strings:
         case op::all_different: case op::some_same:
            return OperandLayout::CountedList;
         case op::sum_list: case op::and_list: case op::or_list:
            return OperandLayout::SumList;
         case op::piecewise_linear:
            return OperandLayout::PiecewiseLinear;
         default:
            return OperandLayout::Invalid;
      }
   }

   namespace {
      [[noreturn]] void unsupported_operation(int opcode) {
         throw std::runtime_error("cppasl: .nl operation o" + std::to_string(opcode) +
            " (logical/counting/symbolic) is not supported in algebraic functions");
      }

      /// Scalar semantics of an .nl operation; only used to fold constant subexpressions.
      double evaluate_operation(int opcode, const double* a, std::uint32_t count) {
         namespace op = nl_opcode;
         switch (opcode) {
            case op::plus: return a[0] + a[1];
            case op::minus: return a[0] - a[1];
            case op::multiply: return a[0] * a[1];
            case op::divide: return a[0] / a[1];
            case op::remainder: return std::fmod(a[0], a[1]);
            case op::power: return std::pow(a[0], a[1]);
            case op::less: return (a[0] > a[1]) ? a[0] - a[1] : 0.;
            case op::minimum_list: return *std::min_element(a, a + count);
            case op::maximum_list: return *std::max_element(a, a + count);
            case op::floor: return std::floor(a[0]);
            case op::ceil: return std::ceil(a[0]);
            case op::abs: return std::fabs(a[0]);
            case op::negate: return -a[0];
            case op::logical_or: return (a[0] != 0. || a[1] != 0.) ? 1. : 0.;
            case op::logical_and: return (a[0] != 0. && a[1] != 0.) ? 1. : 0.;
            case op::less_than: return (a[0] < a[1]) ? 1. : 0.;
            case op::less_equal: return (a[0] <= a[1]) ? 1. : 0.;
            case op::equal: return (a[0] == a[1]) ? 1. : 0.;
            case op::greater_equal: return (a[0] >= a[1]) ? 1. : 0.;
            case op::greater_than: return (a[0] > a[1]) ? 1. : 0.;
            case op::not_equal: return (a[0] != a[1]) ? 1. : 0.;
            case op::logical_not: return (a[0] == 0.) ? 1. : 0.;
            case op::if_then_else: case op::implies: return (a[0] != 0.) ? a[1] : a[2];
            case op::tanh: return std::tanh(a[0]);
            case op::tan: return std::tan(a[0]);
            case op::sqrt: return std::sqrt(a[0]);
            case op::sinh: return std::sinh(a[0]);
            case op::sin: return std::sin(a[0]);
            case op::log10: return std::log10(a[0]);
            case op::log: return std::log(a[0]);
            case op::exp: return std::exp(a[0]);
            case op::cosh: return std::cosh(a[0]);
            case op::cos: return std::cos(a[0]);
            case op::atanh: return std::atanh(a[0]);
            case op::atan2: return std::atan2(a[0], a[1]);
            case op::atan: return std::atan(a[0]);
            case op::asinh: return std::asinh(a[0]);
            case op::asin: return std::asin(a[0]);
            case op::acosh: return std::acosh(a[0]);
            case op::acos: return std::acos(a[0]);
            case op::sum_list: { double sum = 0.; for (std::uint32_t k = 0; k < count; ++k) sum += a[k]; return sum; }
            case op::integer_divide: return std::trunc(a[0] / a[1]);
            case op::round: { const double s = std::pow(10., a[1]); return std::round(a[0] * s) / s; }
            case op::trunc: { const double s = std::pow(10., a[1]); return std::trunc(a[0] * s) / s; }
            case op::precision: {
               if (a[0] == 0.) return 0.;
               const double s = std::pow(10., a[1] - 1. - std::floor(std::log10(std::fabs(a[0]))));
               return std::round(a[0] * s) / s;
            }
            case op::and_list: { for (std::uint32_t k = 0; k < count; ++k) if (a[k] == 0.) return 0.; return 1.; }
            case op::or_list: { for (std::uint32_t k = 0; k < count; ++k) if (a[k] != 0.) return 1.; return 0.; }
            case op::logistic: return 1. / (1. + std::exp(-a[0]));
            case op::sign_power: return std::copysign(std::pow(std::fabs(a[0]), a[1]), a[0]);
            default: unsupported_operation(opcode);
         }
      }
   } // namespace

   // --------------------------------------------------------------------------------------------------------------
   // ExpressionAnalysis
   // --------------------------------------------------------------------------------------------------------------

   ExpressionAnalysis::ExpressionAnalysis(const ExpressionArena& arena, int number_variables):
         arena(arena), number_variables(number_variables) {
   }

   void ExpressionAnalysis::forget_nodes_from(std::size_t node_count) {
      if (this->degree_cache.size() > node_count) this->degree_cache.resize(node_count);
   }

   int ExpressionAnalysis::degree(NodeIndex node) {
      if (node >= this->degree_cache.size()) this->degree_cache.resize(this->arena.nodes.size(), -1);
      std::int8_t& cached = this->degree_cache[node];
      if (cached < 0) {
         const int computed = this->compute_degree(this->arena.nodes[node]);
         this->degree_cache[node] = static_cast<std::int8_t>(computed); // the recursion may have resized the cache
      }
      return this->degree_cache[node];
   }

   int ExpressionAnalysis::compute_degree(const ExpressionNode& node) {
      namespace op = nl_opcode;
      constexpr int np = non_polynomial;
      switch (node.opcode) {
         case constant_opcode: return 0;
         case variable_opcode: {
            if (static_cast<int>(node.first_argument) < this->number_variables) return 1;
            const DefinedVariable& defined = this->arena.defined_variables[node.first_argument - this->number_variables];
            const int linear_degree = (defined.linear_end > defined.linear_begin) ? 1 : 0;
            return std::max(linear_degree, this->degree(defined.expression));
         }
         case op::plus: case op::minus:
            return std::max(this->degree(this->arena.operand(node, 0)), this->degree(this->arena.operand(node, 1)));
         case op::sum_list: {
            int result = 0;
            for (std::uint32_t k = 0; k < node.argument_count && result < np; ++k) {
               result = std::max(result, this->degree(this->arena.operand(node, k)));
            }
            return result;
         }
         case op::negate: return this->degree(this->arena.operand(node, 0));
         case op::multiply:
            return std::min(np, this->degree(this->arena.operand(node, 0)) + this->degree(this->arena.operand(node, 1)));
         case op::divide: {
            const int numerator = this->degree(this->arena.operand(node, 0));
            return (this->degree(this->arena.operand(node, 1)) == 0) ? numerator : np;
         }
         case op::power: {
            const NodeIndex base = this->arena.operand(node, 0), exponent = this->arena.operand(node, 1);
            const int base_degree = this->degree(base), exponent_degree = this->degree(exponent);
            if (exponent_degree > 0) return np;
            if (base_degree == 0) return 0;
            const double p = this->constant_value(exponent);
            if (p == 0.) return 0;
            if (p == 1.) return base_degree;
            if (p == 2.) return std::min(np, 2 * base_degree);
            return np;
         }
         default: {
            // any other operation is constant iff all its operands are, and non-polynomial otherwise
            for (std::uint32_t k = 0; k < node.argument_count; ++k) {
               if (this->degree(this->arena.operand(node, k)) != 0) return np;
            }
            return 0;
         }
      }
   }

   double ExpressionAnalysis::constant_value(NodeIndex node_index) {
      const ExpressionNode& node = this->arena.nodes[node_index];
      if (node.opcode == constant_opcode) return node.value;
      if (node.opcode == variable_opcode) {
         if (static_cast<int>(node.first_argument) < this->number_variables) {
            throw std::logic_error("cppasl: constant_value called on a variable");
         }
         // a constant defined variable has no linear part
         return this->constant_value(this->arena.defined_variables[node.first_argument - this->number_variables].expression);
      }
      std::vector<double> operand_values(node.argument_count);
      if (node.opcode == nl_opcode::piecewise_linear) {
         // operands: 2k-1 constants (slopes and breakpoints), then the argument
         operand_values[0] = static_cast<double>((node.argument_count) / 2);
         for (std::uint32_t k = 0; k + 1 < node.argument_count; ++k) {
            operand_values[k + 1] = this->constant_value(this->arena.operand(node, k));
         }
         return evaluate_piecewise_linear(operand_values.data(),
            this->constant_value(this->arena.operand(node, node.argument_count - 1)), nullptr);
      }
      for (std::uint32_t k = 0; k < node.argument_count; ++k) {
         operand_values[k] = this->constant_value(this->arena.operand(node, k));
      }
      return evaluate_operation(node.opcode, operand_values.data(), node.argument_count);
   }

   // --------------------------------------------------------------------------------------------------------------
   // FunctionDecomposer
   // --------------------------------------------------------------------------------------------------------------

   FunctionDecomposer::FunctionDecomposer(const ExpressionArena& arena, ExpressionAnalysis& analysis,
         int number_variables, bool detect_quadratic_structure):
         arena(arena), analysis(analysis), number_variables(number_variables),
         maximum_polynomial_degree(detect_quadratic_structure ? 2 : 1) {
   }

   void FunctionDecomposer::decompose(NodeIndex root, DecomposedFunction& result) {
      result.clear();
      this->split(root, 1., result);
   }

   void FunctionDecomposer::split(NodeIndex node_index, double scale, DecomposedFunction& result) {
      namespace op = nl_opcode;
      if (this->analysis.degree(node_index) <= this->maximum_polynomial_degree) {
         this->expand_polynomial(node_index, scale, result);
         return;
      }
      const ExpressionNode& node = this->arena.nodes[node_index];
      auto operand = [&](std::uint32_t k) { return this->arena.operand(node, k); };
      switch (node.opcode) {
         case op::plus: this->split(operand(0), scale, result); this->split(operand(1), scale, result); return;
         case op::minus: this->split(operand(0), scale, result); this->split(operand(1), -scale, result); return;
         case op::sum_list:
            for (std::uint32_t k = 0; k < node.argument_count; ++k) this->split(operand(k), scale, result);
            return;
         case op::negate: this->split(operand(0), -scale, result); return;
         case op::multiply:
            if (this->analysis.degree(operand(0)) == 0) {
               this->split(operand(1), scale * this->analysis.constant_value(operand(0)), result);
               return;
            }
            if (this->analysis.degree(operand(1)) == 0) {
               this->split(operand(0), scale * this->analysis.constant_value(operand(1)), result);
               return;
            }
            break;
         case op::divide:
            if (this->analysis.degree(operand(1)) == 0) {
               this->split(operand(0), scale / this->analysis.constant_value(operand(1)), result);
               return;
            }
            break;
         case variable_opcode: { // non-polynomial defined variable: split its linear part and expression
            const DefinedVariable& defined = this->arena.defined_variables[node.first_argument - this->number_variables];
            for (std::uint32_t k = defined.linear_begin; k < defined.linear_end; ++k) {
               const LinearTerm& term = this->arena.defined_variable_linear_terms[k];
               result.linear_terms.push_back({term.variable, scale * term.coefficient});
            }
            this->split(defined.expression, scale, result);
            return;
         }
         default: break;
      }
      result.elements.push_back({node_index, scale});
   }

   void FunctionDecomposer::expand_polynomial(NodeIndex node_index, double scale, DecomposedFunction& result) {
      namespace op = nl_opcode;
      const int degree = this->analysis.degree(node_index);
      if (degree <= 1) {
         this->collect_affine(node_index, scale, result.constant, result.linear_terms);
         return;
      }
      const ExpressionNode& node = this->arena.nodes[node_index];
      auto operand = [&](std::uint32_t k) { return this->arena.operand(node, k); };
      switch (node.opcode) {
         case op::plus: this->expand_polynomial(operand(0), scale, result); this->expand_polynomial(operand(1), scale, result); return;
         case op::minus: this->expand_polynomial(operand(0), scale, result); this->expand_polynomial(operand(1), -scale, result); return;
         case op::sum_list:
            for (std::uint32_t k = 0; k < node.argument_count; ++k) this->expand_polynomial(operand(k), scale, result);
            return;
         case op::negate: this->expand_polynomial(operand(0), -scale, result); return;
         case op::multiply: {
            if (this->analysis.degree(operand(0)) == 0) {
               this->expand_polynomial(operand(1), scale * this->analysis.constant_value(operand(0)), result);
            }
            else if (this->analysis.degree(operand(1)) == 0) {
               this->expand_polynomial(operand(0), scale * this->analysis.constant_value(operand(1)), result);
            }
            else {
               this->expand_product(operand(0), operand(1), scale, result); // two affine factors
            }
            return;
         }
         case op::divide:
            this->expand_polynomial(operand(0), scale / this->analysis.constant_value(operand(1)), result);
            return;
         case op::power: // exponent 1 or 2 (see degree)
            if (this->analysis.constant_value(operand(1)) == 1.) this->expand_polynomial(operand(0), scale, result);
            else this->expand_product(operand(0), operand(0), scale, result);
            return;
         case variable_opcode: {
            const DefinedVariable& defined = this->arena.defined_variables[node.first_argument - this->number_variables];
            for (std::uint32_t k = defined.linear_begin; k < defined.linear_end; ++k) {
               const LinearTerm& term = this->arena.defined_variable_linear_terms[k];
               result.linear_terms.push_back({term.variable, scale * term.coefficient});
            }
            this->expand_polynomial(defined.expression, scale, result);
            return;
         }
         default:
            throw std::logic_error("cppasl: unexpected operation in a quadratic expression");
      }
   }

   bool FunctionDecomposer::as_scaled_variable(NodeIndex node_index, double& coefficient, int& variable) const {
      const ExpressionNode* node = &this->arena.nodes[node_index];
      coefficient = 1.;
      if (node->opcode == nl_opcode::negate) {
         coefficient = -1.;
         node = &this->arena.nodes[this->arena.operand(*node, 0)];
      }
      else if (node->opcode == nl_opcode::multiply) {
         const ExpressionNode& left = this->arena.nodes[this->arena.operand(*node, 0)];
         const ExpressionNode& right = this->arena.nodes[this->arena.operand(*node, 1)];
         if (left.opcode == constant_opcode) { coefficient = left.value; node = &right; }
         else if (right.opcode == constant_opcode) { coefficient = right.value; node = &left; }
         else return false;
      }
      if (node->opcode != variable_opcode || static_cast<int>(node->first_argument) >= this->number_variables) return false;
      variable = static_cast<int>(node->first_argument);
      return true;
   }

   void FunctionDecomposer::expand_product(NodeIndex left, NodeIndex right, double scale, DecomposedFunction& result) {
      // fast path (the vast majority of AMPL quadratic terms): (a x_i) * (b x_j)
      double left_coefficient = 0., right_coefficient = 0.;
      int i = 0, j = 0;
      if (this->as_scaled_variable(left, left_coefficient, i) && this->as_scaled_variable(right, right_coefficient, j)) {
         result.quadratic_terms.push_back({std::max(i, j), std::min(i, j), scale * left_coefficient * right_coefficient});
         return;
      }
      // general case: (a0 + sum a_p x_p) * (b0 + sum b_q x_q), both affine forms stored on the scratch stack
      std::vector<LinearTerm>& scratch = this->affine_scratch;
      const std::size_t left_begin = scratch.size();
      double left_constant = 0., right_constant = 0.;
      this->collect_affine(left, 1., left_constant, scratch);
      const std::size_t right_begin = scratch.size();
      this->collect_affine(right, 1., right_constant, scratch);
      const std::size_t right_end = scratch.size();

      result.constant += scale * left_constant * right_constant;
      if (right_constant != 0.) {
         for (std::size_t p = left_begin; p < right_begin; ++p) {
            result.linear_terms.push_back({scratch[p].variable, scale * right_constant * scratch[p].coefficient});
         }
      }
      if (left_constant != 0.) {
         for (std::size_t q = right_begin; q < right_end; ++q) {
            result.linear_terms.push_back({scratch[q].variable, scale * left_constant * scratch[q].coefficient});
         }
      }
      for (std::size_t p = left_begin; p < right_begin; ++p) {
         for (std::size_t q = right_begin; q < right_end; ++q) {
            const int i = scratch[p].variable, j = scratch[q].variable;
            result.quadratic_terms.push_back({std::max(i, j), std::min(i, j), scale * scratch[p].coefficient * scratch[q].coefficient});
         }
      }
      scratch.resize(left_begin);
   }

   void FunctionDecomposer::collect_affine(NodeIndex node_index, double scale, double& constant, std::vector<LinearTerm>& terms) {
      namespace op = nl_opcode;
      if (this->analysis.degree(node_index) == 0) {
         constant += scale * this->analysis.constant_value(node_index);
         return;
      }
      const ExpressionNode& node = this->arena.nodes[node_index];
      auto operand = [&](std::uint32_t k) { return this->arena.operand(node, k); };
      switch (node.opcode) {
         case variable_opcode: {
            if (static_cast<int>(node.first_argument) < this->number_variables) {
               terms.push_back({static_cast<int>(node.first_argument), scale});
               return;
            }
            const DefinedVariable& defined = this->arena.defined_variables[node.first_argument - this->number_variables];
            for (std::uint32_t k = defined.linear_begin; k < defined.linear_end; ++k) {
               const LinearTerm& term = this->arena.defined_variable_linear_terms[k];
               terms.push_back({term.variable, scale * term.coefficient});
            }
            this->collect_affine(defined.expression, scale, constant, terms);
            return;
         }
         case op::plus: this->collect_affine(operand(0), scale, constant, terms); this->collect_affine(operand(1), scale, constant, terms); return;
         case op::minus: this->collect_affine(operand(0), scale, constant, terms); this->collect_affine(operand(1), -scale, constant, terms); return;
         case op::sum_list:
            for (std::uint32_t k = 0; k < node.argument_count; ++k) this->collect_affine(operand(k), scale, constant, terms);
            return;
         case op::negate: this->collect_affine(operand(0), -scale, constant, terms); return;
         case op::multiply:
            if (this->analysis.degree(operand(0)) == 0) {
               this->collect_affine(operand(1), scale * this->analysis.constant_value(operand(0)), constant, terms);
            }
            else {
               this->collect_affine(operand(0), scale * this->analysis.constant_value(operand(1)), constant, terms);
            }
            return;
         case op::divide:
            this->collect_affine(operand(0), scale / this->analysis.constant_value(operand(1)), constant, terms);
            return;
         case op::power: // exponent 1
            this->collect_affine(operand(0), scale, constant, terms);
            return;
         default:
            throw std::logic_error("cppasl: unexpected operation in an affine expression");
      }
   }

   // --------------------------------------------------------------------------------------------------------------
   // TapeCompiler
   // --------------------------------------------------------------------------------------------------------------

   TapeCompiler::TapeCompiler(const ExpressionArena& arena, ExpressionAnalysis& analysis, int number_variables,
         TapeStorage& storage):
         arena(arena), analysis(analysis), number_variables(number_variables), storage(storage), target(&storage),
         variable_instruction(static_cast<std::size_t>(number_variables), -1) {
   }

   CompiledTape TapeCompiler::compile(NodeIndex root, TapeStorage& destination, bool use_shared_defined_variables,
         bool is_last_pass) {
      this->target = &destination;
      this->use_shared = use_shared_defined_variables;
      if (!this->element_open) { // the passes of one element share its local variables
         this->variable_begin = this->storage.element_variables.size();
         this->element_open = true;
      }
      this->instruction_begin = this->target->instructions.size();
      if (this->defined_variable_instruction.size() < this->arena.defined_variables.size()) {
         this->defined_variable_instruction.resize(this->arena.defined_variables.size(), -1);
      }
      if (this->local_index_of_variable.size() < static_cast<std::size_t>(this->number_variables)) {
         this->local_index_of_variable.resize(static_cast<std::size_t>(this->number_variables), -1);
      }
      this->has_group = this->find_group_sum(root, this->group_sum_node);
      this->group_sum_instruction = -1;
      const std::size_t range_begin = this->target->group_operand_ranges.size();
      const std::int32_t result = this->compile_node(root);
      // the output of the tape must be its last instruction
      if (static_cast<std::size_t>(result) + 1 != this->target->instructions.size() - this->function_begin) {
         this->append(TapeOperation::MultiplyConstant, result, -1, 1.);
      }
      this->forget_instruction_caches();
      CompiledTape tape{this->instruction_begin, this->target->instructions.size() - this->instruction_begin,
         this->variable_begin, this->storage.element_variables.size() - this->variable_begin};
      if (this->group_sum_instruction >= 0) {
         tape.group_sum = this->group_sum_instruction;
         tape.group_range_begin = range_begin;
         tape.group_range_end = this->target->group_operand_ranges.size();
      }
      else this->target->group_operand_ranges.resize(range_begin);
      if (is_last_pass) {
         for (int variable: this->element_variables_seen) this->local_index_of_variable[static_cast<std::size_t>(variable)] = -1;
         this->element_variables_seen.clear();
         this->element_open = false;
      }
      return tape;
   }

   void TapeCompiler::forget_instruction_caches() {
      for (int variable: this->touched_variables) this->variable_instruction[static_cast<std::size_t>(variable)] = -1;
      for (int defined: this->touched_defined_variables) this->defined_variable_instruction[static_cast<std::size_t>(defined)] = -1;
      this->touched_variables.clear();
      this->touched_defined_variables.clear();
   }

   bool TapeCompiler::find_group_sum(NodeIndex root, NodeIndex& sum_node) {
      namespace op = nl_opcode;
      constexpr std::uint32_t minimum_group_size = 3;
      NodeIndex node_index = root;
      for (;;) {
         const ExpressionNode& node = this->arena.nodes[node_index];
         if (node.opcode == op::sum_list) {
            if (node.argument_count < minimum_group_size) return false;
            sum_node = node_index;
            return true;
         }
         if (node.opcode == variable_opcode) return false; // defined variables are not followed
         const OperandLayout layout = operand_layout(node.opcode);
         if (layout == OperandLayout::Unary && node.opcode != op::logical_not && node.opcode != op::floor &&
               node.opcode != op::ceil) {
            node_index = this->arena.operand(node, 0);
            continue;
         }
         // binary operations with one constant operand compile to unary instructions
         if (node.opcode == op::plus || node.opcode == op::minus || node.opcode == op::multiply || node.opcode == op::divide ||
               node.opcode == op::power) {
            const NodeIndex left = this->arena.operand(node, 0), right = this->arena.operand(node, 1);
            const bool left_constant = this->analysis.degree(left) == 0, right_constant = this->analysis.degree(right) == 0;
            if (right_constant && !left_constant) { node_index = left; continue; }
            if (left_constant && !right_constant && node.opcode != op::divide) { node_index = right; continue; }
         }
         return false;
      }
   }

   /// the operands of a group sum are compiled as self-contained, contiguous instruction ranges (their variables and
   /// defined variables are not shared), so that the Hessian of each can be computed on its range only
   std::int32_t TapeCompiler::compile_group_sum(const ExpressionNode& node) {
      const std::size_t scratch_begin = this->operand_scratch.size();
      bool self_contained = true;
      for (std::uint32_t k = 0; k < node.argument_count; ++k) {
         this->forget_instruction_caches();
         const std::size_t begin = this->target->instructions.size();
         const std::int32_t operand = this->compile_node(this->arena.operand(node, k));
         this->operand_scratch.push_back(operand);
         // an operand may be an instruction compiled earlier (e.g. a repeated constant): only proper ranges count
         if (static_cast<std::size_t>(operand) + 1 == this->target->instructions.size() - this->function_begin &&
               this->target->instructions.size() > begin) {
            this->target->group_operand_ranges.emplace_back(begin, this->target->instructions.size());
         }
         else self_contained = false;
      }
      this->forget_instruction_caches();
      const auto first = static_cast<std::int32_t>(this->target->arguments.size());
      this->target->arguments.insert(this->target->arguments.end(), this->operand_scratch.begin() + static_cast<std::ptrdiff_t>(scratch_begin),
         this->operand_scratch.end());
      this->operand_scratch.resize(scratch_begin);
      const std::int32_t sum = this->append(TapeOperation::Sum, first, static_cast<std::int32_t>(node.argument_count), 0.);
      if (self_contained) this->group_sum_instruction = sum;
      return sum;
   }

   std::int32_t TapeCompiler::append(TapeOperation operation, std::int32_t first, std::int32_t second, double constant) {
      if (this->target->instructions.size() >= static_cast<std::size_t>(std::numeric_limits<std::int32_t>::max())) {
         throw std::runtime_error("cppasl: more than 2^31-1 tape instructions are not supported");
      }
      this->target->instructions.push_back({operation, first, second, constant});
      return static_cast<std::int32_t>(this->target->instructions.size() - 1 - this->function_begin);
   }

   /// the local index of a variable in the current element (registered on first use)
   std::int32_t TapeCompiler::register_local(int variable) {
      std::int32_t& local_index = this->local_index_of_variable[static_cast<std::size_t>(variable)];
      if (local_index < 0) {
         local_index = static_cast<std::int32_t>(this->storage.element_variables.size() - this->variable_begin);
         this->element_variables_seen.push_back(variable);
         this->storage.element_variables.push_back(variable);
      }
      return local_index;
   }

   /// a new leaf x + offset (not shared through the cache: other uses of x read the plain leaf)
   std::int32_t TapeCompiler::shifted_variable(int variable, double offset) {
      std::int32_t& local_index = this->local_index_of_variable[static_cast<std::size_t>(variable)];
      const bool is_new = local_index < 0;
      if (is_new) {
         local_index = static_cast<std::int32_t>(this->storage.element_variables.size() - this->variable_begin);
         this->element_variables_seen.push_back(variable);
      }
      const std::int32_t instruction = this->append(TapeOperation::Variable, variable, local_index, offset);
      if (is_new) this->storage.element_variables.push_back(variable);
      return instruction;
   }

   std::int32_t TapeCompiler::variable_reference(int variable) {
      std::int32_t& instruction = this->variable_instruction[static_cast<std::size_t>(variable)];
      if (instruction < 0) {
         // a variable has one local index per element, but may have several instructions (group operands)
         std::int32_t& local_index = this->local_index_of_variable[static_cast<std::size_t>(variable)];
         const bool is_new = local_index < 0;
         if (is_new) {
            local_index = static_cast<std::int32_t>(this->storage.element_variables.size() - this->variable_begin);
            this->element_variables_seen.push_back(variable);
         }
         instruction = this->append(TapeOperation::Variable, variable, local_index, 0.);
         if (is_new) this->storage.element_variables.push_back(variable);
         this->touched_variables.push_back(variable);
      }
      return instruction;
   }

   std::int32_t TapeCompiler::append_nary(TapeOperation operation, const ExpressionNode& node) {
      const std::size_t scratch_begin = this->operand_scratch.size();
      for (std::uint32_t k = 0; k < node.argument_count; ++k) {
         const std::int32_t operand = this->compile_node(this->arena.operand(node, k));
         this->operand_scratch.push_back(operand);
      }
      const auto first = static_cast<std::int32_t>(this->target->arguments.size());
      this->target->arguments.insert(this->target->arguments.end(), this->operand_scratch.begin() + static_cast<std::ptrdiff_t>(scratch_begin),
         this->operand_scratch.end());
      this->operand_scratch.resize(scratch_begin);
      return this->append(operation, first, static_cast<std::int32_t>(node.argument_count), 0.);
   }

   std::int32_t TapeCompiler::compile_node(NodeIndex node_index) {
      namespace op = nl_opcode;
      using O = TapeOperation;
      if (this->analysis.degree(node_index) == 0) {
         return this->append(O::Constant, -1, -1, this->analysis.constant_value(node_index));
      }
      const ExpressionNode& node = this->arena.nodes[node_index];
      auto operand = [&](std::uint32_t k) { return this->compile_node(this->arena.operand(node, k)); };
      auto constant_operand = [&](std::uint32_t k, double& value) {
         const NodeIndex operand_node = this->arena.operand(node, k);
         if (this->analysis.degree(operand_node) != 0) return false;
         value = this->analysis.constant_value(operand_node);
         return true;
      };
      // x + c with x a variable is a shifted variable leaf (one instruction), otherwise an AddConstant instruction
      auto shifted = [&](std::uint32_t k, double c) {
         const ExpressionNode& operand_node = this->arena.nodes[this->arena.operand(node, k)];
         if (operand_node.opcode == variable_opcode && static_cast<int>(operand_node.first_argument) < this->number_variables) {
            return this->shifted_variable(static_cast<int>(operand_node.first_argument), c);
         }
         return this->append(O::AddConstant, operand(k), -1, c);
      };
      auto unary = [&](O operation) { return this->append(operation, operand(0), -1, 0.); };
      auto binary = [&](O operation) {
         const std::int32_t first = operand(0);
         const std::int32_t second = operand(1);
         return this->append(operation, first, second, 0.);
      };

      switch (node.opcode) {
         case variable_opcode: {
            if (static_cast<int>(node.first_argument) < this->number_variables) {
               return this->variable_reference(static_cast<int>(node.first_argument));
            }
            const std::size_t defined_index = node.first_argument - static_cast<std::size_t>(this->number_variables);
            if (this->defined_variable_instruction[defined_index] >= 0) return this->defined_variable_instruction[defined_index];
            if (this->use_shared && this->shared_defined_outputs != nullptr && (*this->shared_defined_outputs)[defined_index] >= 0) {
               // evaluated once per point by its own tape: a leaf reading its value
               const std::int32_t leaf = this->append(O::DefinedValue,
                  static_cast<std::int32_t>((*this->shared_defined_outputs)[defined_index]), static_cast<std::int32_t>(defined_index), 0.);
               // the variables of the defined variable are variables of the element (Jacobian and Hessian patterns)
               const auto [variables_begin, variables_end] = (*this->shared_defined_variable_ranges)[defined_index];
               for (std::size_t v = variables_begin; v < variables_end; ++v) this->register_local(this->storage.element_variables[v]);
               this->defined_variable_instruction[defined_index] = leaf;
               this->touched_defined_variables.push_back(static_cast<int>(defined_index));
               return leaf;
            }
            const DefinedVariable& defined = this->arena.defined_variables[defined_index];
            std::int32_t result = this->compile_node(defined.expression);
            for (std::uint32_t k = defined.linear_begin; k < defined.linear_end; ++k) {
               const LinearTerm& term = this->arena.defined_variable_linear_terms[k];
               const std::int32_t scaled = this->append(O::MultiplyConstant, this->variable_reference(term.variable), -1, term.coefficient);
               result = this->append(O::Add, result, scaled, 0.);
            }
            this->defined_variable_instruction[defined_index] = result;
            this->touched_defined_variables.push_back(static_cast<int>(defined_index));
            return result;
         }
         // binary operations with a constant operand become unary instructions (the constant is not emitted)
         case op::plus: {
            double c = 0.;
            if (constant_operand(0, c)) return shifted(1, c);
            if (constant_operand(1, c)) return shifted(0, c);
            const std::int32_t a = operand(0), b = operand(1);
            return this->append(O::Add, a, b, 0.);
         }
         case op::minus: {
            double c = 0.;
            if (constant_operand(1, c)) return shifted(0, -c);
            if (constant_operand(0, c)) return this->append(O::AddConstant, this->append(O::Negate, operand(1), -1, 0.), -1, c);
            const std::int32_t a = operand(0), b = operand(1);
            return this->append(O::Subtract, a, b, 0.);
         }
         case op::multiply: {
            double c = 0.;
            if (constant_operand(0, c)) return this->append(O::MultiplyConstant, operand(1), -1, c);
            if (constant_operand(1, c)) return this->append(O::MultiplyConstant, operand(0), -1, c);
            const std::int32_t a = operand(0), b = operand(1);
            return this->append(O::Multiply, a, b, 0.);
         }
         case op::divide: {
            double c = 0.;
            if (constant_operand(1, c)) return this->append(O::MultiplyConstant, operand(0), -1, 1. / c);
            const std::int32_t a = operand(0), b = operand(1);
            return this->append(O::Divide, a, b, 0.);
         }
         case op::power: {
            double c = 0.;
            if (constant_operand(1, c)) {
               const std::int32_t a = operand(0);
               if (c == 1.) return a;
               if (c == 2.) return this->append(O::Square, a, -1, 0.);
               if (c == std::trunc(c) && std::fabs(c) <= 16.) return this->append(O::PowerInteger, a, -1, c);
               if (c == 0.5) return this->append(O::Sqrt, a, -1, 0.);
               return this->append(O::PowerConstantExponent, a, -1, c);
            }
            if (constant_operand(0, c)) return this->append(O::PowerConstantBase, operand(1), -1, c);
            const std::int32_t a = operand(0), b = operand(1);
            return this->append(O::Power, a, b, 0.);
         }
         case op::sign_power: {
            double c = 0.;
            if (!constant_operand(1, c)) throw std::runtime_error("cppasl: signpow with a variable exponent is not supported");
            return this->append(O::SignPowerConstant, operand(0), -1, c);
         }
         case op::negate: return unary(O::Negate);
         case op::abs: return unary(O::Abs);
         case op::sqrt: return unary(O::Sqrt);
         case op::exp: return unary(O::Exp);
         case op::log: return unary(O::Log);
         case op::log10: return unary(O::Log10);
         case op::sin: return unary(O::Sin);
         case op::cos: return unary(O::Cos);
         case op::tan: return unary(O::Tan);
         case op::sinh: return unary(O::Sinh);
         case op::cosh: return unary(O::Cosh);
         case op::tanh: return unary(O::Tanh);
         case op::asin: return unary(O::Asin);
         case op::acos: return unary(O::Acos);
         case op::atan: return unary(O::Atan);
         case op::asinh: return unary(O::Asinh);
         case op::acosh: return unary(O::Acosh);
         case op::atanh: return unary(O::Atanh);
         case op::logistic: return unary(O::Logistic);
         case op::floor: return unary(O::Floor);
         case op::ceil: return unary(O::Ceil);
         case op::logical_not: return unary(O::Not);
         case op::atan2: return binary(O::Atan2);
         case op::less: return binary(O::Less);
         case op::remainder: return binary(O::Remainder);
         case op::round: return binary(O::Round);
         case op::trunc: return binary(O::Trunc);
         case op::precision: return binary(O::Precision);
         case op::integer_divide: return binary(O::IntegerDivide);
         case op::less_than: return binary(O::LessThan);
         case op::less_equal: return binary(O::LessEqual);
         case op::equal: return binary(O::Equal);
         case op::greater_equal: return binary(O::GreaterEqual);
         case op::greater_than: return binary(O::GreaterThan);
         case op::not_equal: return binary(O::NotEqual);
         case op::logical_and: return binary(O::And);
         case op::logical_or: return binary(O::Or);
         case op::sum_list:
            if (this->has_group && node_index == this->group_sum_node) return this->compile_group_sum(node);
            return this->append_nary(O::Sum, node);
         case op::minimum_list: return this->append_nary(O::Minimum, node);
         case op::maximum_list: return this->append_nary(O::Maximum, node);
         case op::and_list: return this->append_nary(O::AndList, node);
         case op::or_list: return this->append_nary(O::OrList, node);
         case op::if_then_else: case op::implies: return this->append_nary(O::IfThenElse, node);
         case op::piecewise_linear: {
            // operands: 2k-1 constant slopes/breakpoints, then the argument
            const std::int32_t argument = operand(node.argument_count - 1);
            const auto data_offset = static_cast<std::int32_t>(this->target->constants.size());
            this->target->constants.push_back(static_cast<double>(node.argument_count / 2));
            for (std::uint32_t k = 0; k + 1 < node.argument_count; ++k) {
               this->target->constants.push_back(this->analysis.constant_value(this->arena.operand(node, k)));
            }
            return this->append(O::PiecewiseLinear, argument, data_offset, 0.);
         }
         default:
            unsupported_operation(node.opcode);
      }
   }

} // namespace cppasl::detail
