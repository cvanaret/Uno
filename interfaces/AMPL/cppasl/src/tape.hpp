#pragma once

#include <cstddef>
#include <array>
#include <cstdint>
#include <utility>

namespace cppasl::detail {

   /// Operations of the AD tapes. Tapes are flat, topologically sorted arrays of instructions; the result of
   /// instruction i is value[i]. Operand fields are local instruction indices, except where noted.
   enum class TapeOperation : std::uint8_t {
      // leaves
      Constant,          ///< value = constant
      Variable,          ///< value = x[first] + constant; second = local variable index
      DefinedValue,      ///< value = values[first], the output of a shared defined variable's tape; second = its index
      // unary: operand = first; parameter = constant
      Negate, AddConstant, MultiplyConstant, Square,
      PowerInteger,      ///< u^constant, constant integer in [-16, 16] (repeated multiplications instead of pow)
      PowerConstantExponent, PowerConstantBase,
      Abs, Sqrt, Exp, Log, Log10, Sin, Cos, Tan, Sinh, Cosh, Tanh, Asin, Acos, Atan, Asinh, Acosh, Atanh,
      Logistic,
      SignPowerConstant, ///< sign(u) |u|^constant
      PiecewiseLinear,   ///< second = offset of {number of slopes, s0, b0, s1, b1, ..., s_last} in tape constants
      // binary: operands = first, second
      Add, Subtract, Multiply, Divide, Power, Atan2, Less, Remainder,
      // n-ary: operands = tape_arguments[first, first + second)
      Sum, Minimum, Maximum,
      IfThenElse,        ///< n-ary with 3 operands: condition, then, else
      // piecewise constant (zero derivatives): unary, binary or n-ary as indicated
      Floor, Ceil,                                                   // unary
      Round, Trunc, Precision, IntegerDivide,                        // binary
      LessThan, LessEqual, Equal, GreaterEqual, GreaterThan, NotEqual, And, Or, // binary
      Not,                                                           // unary
      AndList, OrList,                                               // n-ary
      NumberOperations
   };

   /// Shape of an operation for the derivative sweeps
   enum class DerivativeCategory : std::uint8_t { Leaf, Unary, Binary, Sum, Select, Flat };

   namespace tables {
      constexpr DerivativeCategory category_of(TapeOperation operation) {
         using O = TapeOperation;
         switch (operation) {
            case O::Constant: case O::Variable: case O::DefinedValue: return DerivativeCategory::Leaf;
            case O::Add: case O::Subtract: case O::Multiply: case O::Divide: case O::Power: case O::Atan2: case O::Less:
            case O::Remainder: return DerivativeCategory::Binary;
            case O::Sum: return DerivativeCategory::Sum;
            case O::Minimum: case O::Maximum: case O::IfThenElse: return DerivativeCategory::Select;
            case O::Floor: case O::Ceil: case O::Round: case O::Trunc: case O::Precision: case O::IntegerDivide:
            case O::LessThan: case O::LessEqual: case O::Equal: case O::GreaterEqual: case O::GreaterThan:
            case O::NotEqual: case O::And: case O::Or: case O::Not: case O::AndList: case O::OrList:
            case O::NumberOperations: return DerivativeCategory::Flat;
            default: return DerivativeCategory::Unary;
         }
      }
      constexpr int arity_of(TapeOperation operation) {
         using O = TapeOperation;
         switch (category_of(operation)) {
            case DerivativeCategory::Leaf: return 0;
            case DerivativeCategory::Unary: return 1;
            case DerivativeCategory::Binary: return 2;
            case DerivativeCategory::Sum: case DerivativeCategory::Select: return -1;
            default:
               switch (operation) {
                  case O::Floor: case O::Ceil: case O::Not: return 1;
                  case O::AndList: case O::OrList: return -1;
                  default: return 2;
               }
         }
      }
      template <typename T, typename F, std::size_t... I>
      constexpr auto make_table(F function, std::index_sequence<I...>) {
         return std::array<T, sizeof...(I)>{function(static_cast<TapeOperation>(I))...};
      }
      constexpr std::size_t number_operations = static_cast<std::size_t>(TapeOperation::NumberOperations) + 1;
      inline constexpr auto categories = make_table<DerivativeCategory>(category_of, std::make_index_sequence<number_operations>{});
      inline constexpr auto arities = make_table<int>(arity_of, std::make_index_sequence<number_operations>{});
   } // namespace tables

   /// shape of an operation (table lookup)
   inline DerivativeCategory derivative_category(TapeOperation operation) {
      return tables::categories[static_cast<std::size_t>(operation)];
   }
   /// number of instruction operands: 0 (leaf), 1, 2, or -1 (n-ary, operands in the tape arguments)
   inline int operand_arity(TapeOperation operation) { return tables::arities[static_cast<std::size_t>(operation)]; }

   /// position of (a, b), a >= b, in the packed lower triangle of a k x k matrix stored column by column
   inline std::size_t packed_lower_index(std::size_t a, std::size_t b, std::size_t k) {
      return b * (2 * k - b + 1) / 2 + (a - b);
   }

   /// u^n for a small integer n by binary exponentiation
   inline double integer_power(double u, int n) {
      const bool inverse = n < 0;
      unsigned exponent = static_cast<unsigned>(inverse ? -n : n);
      double result = 1., base = u;
      while (exponent != 0) {
         if (exponent & 1u) result *= base;
         base *= base;
         exponent >>= 1u;
      }
      return inverse ? 1. / result : result;
   }

   struct TapeInstruction {
      TapeOperation operation;
      std::int32_t first;
      std::int32_t second;
      double constant;
   };

   /// First and second partial derivatives of an instruction w.r.t. its (at most two) operands
   struct InstructionPartials {
      double first;         ///< d/du
      double second;        ///< d/dw
      double first_first;   ///< d2/du2
      double first_second;  ///< d2/dudw
      double second_second; ///< d2/dw2
   };

   /// A range of instructions [begin, end) of the model's tape storage. Operand fields and the n-ary arguments are
   /// global instruction indices, so that consecutive elements (of one function, or of many constraints) are swept
   /// in one loop for values and gradients, or one element at a time for Hessians. The result of an element is its
   /// last instruction.
   struct TapeView {
      const TapeInstruction* instructions; ///< all the instructions of the model
      std::size_t begin;
      std::size_t end;
      const std::int32_t* arguments;  ///< operands of n-ary instructions
      const double* constants;        ///< piecewise-linear data
   };

   // Forward sweeps: values (and first partials, 2 per instruction) indexed like the instructions.
   void evaluate_tape_values(const TapeView& tape, const double* x, double* values);
   /// also writes the partials of the selection instructions (min, max, if) w.r.t. each operand, indexed like the
   /// tape arguments, in operand_partials
   void evaluate_tape_values_and_partials(const TapeView& tape, const double* x, double* values, double* partials,
      double* operand_partials, const double* defined_values = nullptr);

   /// One step of the reverse (adjoint) program: adjoints[target] += adjoints[source] * partials[partial].
   /// The partial index addresses [instruction partials (2 per instruction) | operand partials (1 per tape
   /// argument) | the constant 1]. A function's program lists its edges by decreasing source: the reverse sweep is a
   /// branch-free loop of multiply-adds (like ASL's derp lists), with no dispatch on the operations.
   /// The first edge into a target assigns instead of accumulating (top bit of `partial`), so the adjoints need no
   /// zeroing; edges from instructions that cannot receive an adjoint are not generated.
   struct ReverseEdge {
      std::uint32_t source;
      std::uint32_t target;
      std::uint32_t partial; ///< index into the partials | assign_flag
      static constexpr std::uint32_t assign_flag = 0x80000000u;
   };
   /// a variable leaf and the gradient position its adjoint is added to
   struct VariableOccurrence {
      std::uint32_t instruction;
      std::uint32_t position;
   };
   inline void propagate_adjoints(const ReverseEdge* edge, const ReverseEdge* end, double* adjoints, const double* partials) {
      for (; edge != end; ++edge) {
         const double contribution = adjoints[edge->source] * partials[edge->partial & ~ReverseEdge::assign_flag];
         double& target = adjoints[edge->target];
         target = ((edge->partial & ReverseEdge::assign_flag) ? 0. : target) + contribution; // select, not a branch
      }
   }
   // Second-order sweeps of one element: scratch arrays are indexed from tape.begin.
   /// first and second partial derivatives
   void evaluate_tape_partials(const TapeView& tape, const double* values, InstructionPartials* partials, bool second_order);
   /// adjoints w.r.t. the last instruction
   void evaluate_tape_adjoints(const TapeView& tape, const double* values, const InstructionPartials* partials, double* adjoints);
   /// number of directions processed simultaneously by the vector forward-over-reverse Hessian sweeps
   constexpr std::size_t maximum_hessian_directions = 8;
   /// A DefinedValue leaf of an element (its index is stored in the leaf's constant field): the pairs (element
   /// local variable, index of the defined variable's partial derivative in the gradients) of the defined variable.
   struct DefinedLeafLink {
      std::size_t pair_begin;
      std::uint32_t pair_count;
   };
   /// Chain rule through the shared defined variables in the element Hessians: H = J^T H_phi J + sum_v a_v H_v, where
   /// J holds the gradients of the defined variables and the curvature terms a_v H_v are accumulated in weights
   /// (indexed by the defined variable's output instruction) and added once for all functions.
   struct DefinedHessianLinks {
      const DefinedLeafLink* leaves;
      const std::uint32_t* pair_locals;
      const std::uint32_t* pair_gradients;
      const double* gradients;
      double* weights;
   };
   /// packed_lower += factor * Hessian w.r.t. the local variables (packed lower triangle column by column, see
   /// packed_lower_index); tangents and second_order_adjoints need maximum_hessian_directions entries per instruction
   void evaluate_tape_hessian(const TapeView& tape, const double* values, InstructionPartials* partials,
      double* adjoints, double* tangents, double* second_order_adjoints, std::size_t number_variables, double factor,
      double* packed_lower, const DefinedHessianLinks* links = nullptr);
   /// Same for a group element phi(sum_k t_k) whose operands t_k occupy self-contained instruction ranges:
   /// H = phi'' g g^T + phi' sum_k H_k (g = gradient of the sum), where each H_k costs sweeps over its own range only
   /// (ASL's group partial separability). direction_of_local must be -1 on entry (restored on exit). Returns false
   /// if the tape does not have the group shape (then nothing is accumulated).
   bool evaluate_group_hessian(const TapeView& tape, std::size_t group_sum,
      const std::pair<std::size_t, std::size_t>* operand_ranges, std::size_t number_ranges, const double* values,
      InstructionPartials* partials, double* adjoints, double* tangents, double* second_order_adjoints,
      std::int32_t* direction_of_local, std::int32_t* local_of_direction, double* sum_gradient,
      std::size_t number_variables, double factor, double* packed_lower);
   /// Links of an evaluation tape to the shared defined variables (indexed by their output instruction): the
   /// tangents of the defined variables along the direction, and the first- and second-order adjoint seeds the
   /// element contributes to them.
   struct DefinedVariableLinks {
      const double* tangents;
      double* adjoint_seeds;
      double* second_order_seeds;
   };
   /// result[variable] += factor * (H d)_variable, with d indexed by local variable. With links, the DefinedValue
   /// leaves read their tangents and accumulate the seeds of the defined variables.
   void evaluate_tape_hessian_vector_product(const TapeView& tape, const double* values, InstructionPartials* partials,
      double* adjoints, double* tangents, double* second_order_adjoints, const double* local_direction, double factor,
      double* result, const DefinedVariableLinks* links = nullptr);
   /// result += second-order adjoints at the variable leaves, for seeds given on entry in adjoints and
   /// second_order_adjoints (indexed from tape.begin), along the direction d (indexed by variable)
   void evaluate_seeded_hessian_vector_product(const TapeView& tape, const double* values, InstructionPartials* partials,
      double* adjoints, double* tangents, double* second_order_adjoints, const double* direction, double* result);

   /// value and slope of an AMPL piecewise-linear term (zero at 0)
   double evaluate_piecewise_linear(const double* data, double x, double* slope);

} // namespace cppasl::detail
