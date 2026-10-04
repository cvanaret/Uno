#pragma once

#include <cstdint>
#include <utility>
#include <vector>
#include "expression_arena.hpp"
#include "tape.hpp"

namespace cppasl::detail {

   /// coefficient * x[row] * x[column]
   struct QuadraticTerm {
      int row;
      int column;
      double coefficient;
   };

   /// f = constant + sum linear + sum quadratic + sum_e scale_e * element_e
   struct DecomposedFunction {
      struct Element {
         NodeIndex root;
         double scale;
      };
      double constant{0.};
      std::vector<LinearTerm> linear_terms;       // may contain duplicates
      std::vector<QuadraticTerm> quadratic_terms; // may contain duplicates, any orientation
      std::vector<Element> elements;

      void clear() {
         this->constant = 0.;
         this->linear_terms.clear();
         this->quadratic_terms.clear();
         this->elements.clear();
      }
   };

   /// Polynomial degree and constant folding of expression nodes (memoized per node).
   class ExpressionAnalysis {
   public:
      static constexpr int non_polynomial = 3; ///< degree > 2 or not a polynomial

      ExpressionAnalysis(const ExpressionArena& arena, int number_variables);
      /// degree in {0, 1, 2} or non_polynomial
      int degree(NodeIndex node);
      /// value of a node of degree 0
      double constant_value(NodeIndex node);
      /// to be called when the arena is truncated
      void forget_nodes_from(std::size_t node_count);

   private:
      const ExpressionArena& arena;
      const int number_variables;
      std::vector<std::int8_t> degree_cache; // -1: unknown
      int compute_degree(const ExpressionNode& node);
   };

   /// Splits a function into its linear, quadratic and non-quadratic (element) parts, following the top-level
   /// sums, negations and scalings. Degree <= 2 subexpressions are expanded symbolically.
   class FunctionDecomposer {
   public:
      FunctionDecomposer(const ExpressionArena& arena, ExpressionAnalysis& analysis, int number_variables,
         bool detect_quadratic_structure);
      void decompose(NodeIndex root, DecomposedFunction& result);

   private:
      const ExpressionArena& arena;
      ExpressionAnalysis& analysis;
      const int number_variables;
      const int maximum_polynomial_degree;
      std::vector<LinearTerm> affine_scratch; // stack of affine forms for products

      void split(NodeIndex node, double scale, DecomposedFunction& result);
      void expand_polynomial(NodeIndex node, double scale, DecomposedFunction& result);
      void expand_product(NodeIndex left, NodeIndex right, double scale, DecomposedFunction& result);
      void collect_affine(NodeIndex node, double scale, double& constant, std::vector<LinearTerm>& terms);
      /// recognizes c * x_i, x_i * c, x_i and -x_i (the operands of AMPL's products), without the general machinery
      bool as_scaled_variable(NodeIndex node, double& coefficient, int& variable) const;
   };

   /// Destination of compiled tapes
   struct TapeStorage {
      std::vector<TapeInstruction> instructions;
      std::vector<std::int32_t> arguments;
      std::vector<double> constants;
      std::vector<int> element_variables;                     ///< global index of each local variable
      std::vector<std::pair<std::size_t, std::size_t>> group_operand_ranges;
   };

   struct CompiledTape {
      std::size_t instruction_begin;
      std::size_t instruction_count;
      std::size_t variable_begin;
      std::size_t variable_count;
      /// group elements phi(sum_k t_k): instruction of the sum (-1 if the element is not a group) and the
      /// [begin, end) instruction ranges of its self-contained operands t_k, in TapeStorage::group_operand_ranges
      std::int64_t group_sum{-1};
      std::size_t group_range_begin{0};
      std::size_t group_range_end{0};
   };

   /// Compiles expression trees into flat tapes (post-order, constants folded, defined variables inlined once
   /// per tape, simple peephole specializations such as u*c, u+c, u^2).
   class TapeCompiler {
   public:
      TapeCompiler(const ExpressionArena& arena, ExpressionAnalysis& analysis, int number_variables, TapeStorage& storage);
      /// Compiles an element at the end of `destination` (operands are global instruction indices of it). An element
      /// may be compiled in several passes (e.g. a Hessian tape and an evaluation tape) that share its local
      /// variables; the local variables of the element are recorded in the main storage. With
      /// use_shared_defined_variables, the defined variables listed in shared_defined_outputs are DefinedValue leaves.
      CompiledTape compile(NodeIndex root, TapeStorage& destination, bool use_shared_defined_variables, bool is_last_pass);
      /// output instruction (main storage) of each defined variable evaluated by its own tape, -1 otherwise
      const std::vector<std::int64_t>* shared_defined_outputs{nullptr};
      /// [begin, end) of the variables of each shared defined variable's tape in the main storage's element_variables
      const std::vector<std::pair<std::size_t, std::size_t>>* shared_defined_variable_ranges{nullptr};

   private:
      const ExpressionArena& arena;
      ExpressionAnalysis& analysis;
      const int number_variables;
      TapeStorage& storage;  ///< main storage: also receives the local variables of the elements
      TapeStorage* target;   ///< where the instructions of the current pass go
      bool use_shared{false};
      bool element_open{false};
      const std::size_t function_begin{0}; ///< indices are global
      std::size_t instruction_begin{0};
      std::size_t variable_begin{0};
      std::vector<std::int32_t> variable_instruction;         // per variable, -1 if not in the current tape
      std::vector<std::int32_t> defined_variable_instruction; // per defined variable
      std::vector<int> touched_variables;
      std::vector<std::int32_t> local_index_of_variable;       // per variable, -1 if not in the current element
      std::vector<int> element_variables_seen;
      NodeIndex group_sum_node{0};
      bool has_group{false};
      std::int64_t group_sum_instruction{-1};
      std::vector<int> touched_defined_variables;
      std::vector<std::int32_t> operand_scratch;

      std::int32_t compile_node(NodeIndex node);
      std::int32_t shifted_variable(int variable, double offset);
      std::int32_t register_local(int variable);
      /// the sum_list node S if root = unary chain(S) with S large enough to benefit from the group structure
      bool find_group_sum(NodeIndex root, NodeIndex& sum_node);
      std::int32_t compile_group_sum(const ExpressionNode& node);
      void forget_instruction_caches();
      std::int32_t append(TapeOperation operation, std::int32_t first, std::int32_t second, double constant);
      std::int32_t append_nary(TapeOperation operation, const ExpressionNode& node);
      std::int32_t variable_reference(int variable);
   };

} // namespace cppasl::detail
