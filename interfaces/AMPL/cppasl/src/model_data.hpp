#pragma once

#include <cstddef>
#include <cstdint>
#include <limits>
#include <mutex>
#include <vector>
#include "cppasl/nl_header.hpp"
#include "cppasl/nl_model.hpp"
#include "function_decomposition.hpp"
#include "tape.hpp"

namespace cppasl::detail {

   constexpr std::uint32_t no_hessian_slot = std::numeric_limits<std::uint32_t>::max();

   /// Location of the parts of one objective or constraint inside the flat model arrays.
   /// Linear part: for an objective, [linear_begin, linear_end) in objective_linear_*; for a constraint,
   /// the row of the CSR Jacobian (the linear coefficients are aligned with the Jacobian pattern).
   /// Quadratic part: [quadratic_begin, quadratic_off_diagonal_begin) diagonal, then strictly lower entries.
   struct FunctionBlock {
      double constant{0.};
      std::size_t linear_begin{0}, linear_end{0};
      std::size_t quadratic_begin{0}, quadratic_off_diagonal_begin{0}, quadratic_end{0};
      std::size_t element_begin{0}, element_end{0};
      std::size_t tape_begin{0}, tape_end{0}; ///< instructions of all the elements (one contiguous block)
      std::size_t edge_begin{0}, edge_end{0};             ///< reverse program, in ModelData::reverse_edges
      std::size_t defined_occurrence_begin{0}, defined_occurrence_end{0}; ///< in ModelData::defined_occurrences
      std::size_t occurrence_begin{0}, occurrence_end{0}; ///< variable leaves, in ModelData::variable_occurrences

      [[nodiscard]] bool has_quadratic_part() const { return this->quadratic_end > this->quadratic_begin; }
      [[nodiscard]] bool has_elements() const { return this->element_end > this->element_begin; }
   };

   /// A non-quadratic term scale * tape(x) of the top-level sum of a function
   /// A non-quadratic term scale * tape(x) of the top-level sum of a function. 32-bit fields (the tapes are limited
   /// to 2^31 instructions, the Hessian to 2^32 slots): the elements are streamed through by several loops.
   struct NonlinearElement {
      std::uint32_t function_index;     ///< into ModelData::functions
      std::uint32_t instruction_begin;  ///< into the tape arrays and the workspace tape values
      std::uint32_t instruction_count;
      std::uint32_t variable_begin;     ///< into tapes.element_variables
      std::uint32_t variable_count;
      std::uint32_t hessian_slot_begin; ///< into element_hessian_slots (variable_count * (variable_count + 1) / 2 slots)
      double scale;
      /// group elements phi(sum_k t_k): the sum instruction (none: not a group) and its operand ranges in the
      /// group_operand_ranges of the storage of the Hessian tape
      std::uint32_t group_sum{no_group};
      std::uint32_t group_range_begin{0}, group_range_end{0};
      static constexpr std::uint32_t no_group = static_cast<std::uint32_t>(-1);
      /// the element reads shared defined variables (DefinedValue leaves): its derivatives use the chain rule through
      /// their tapes
      bool has_defined_leaves{false};
   };

   /// a nonlinear defined variable used by several functions, evaluated once per point by its own tape (in
   /// ModelData::defined_storage, where all of them are adjacent: one sweep evaluates them all), with its gradient
   struct DefinedTape {
      std::size_t instruction_begin, instruction_count;
      std::size_t variable_begin, variable_count; ///< in tapes.element_variables
      std::size_t gradient_begin;                 ///< in the workspace's defined_gradients
      std::size_t hessian_slot_begin{0};          ///< in defined_hessian_slots (packed lower triangle)
      bool contributes_to_hessian{false};         ///< read by a function of the Lagrangian
   };

   /// a DefinedValue leaf of a function: adjoint * gradient(defined variable) is added at position_begin...
   struct DefinedOccurrence {
      std::uint32_t instruction;
      std::uint32_t variable_count;  ///< of the defined variable
      std::size_t position_begin;    ///< in defined_occurrence_positions, one per variable of the defined variable
      std::size_t gradient_begin;    ///< in the defined gradients
   };

   /// consecutive nonlinear constraints whose element tapes are adjacent in the storage
   struct ConstraintRun {
      std::size_t tape_begin, tape_end;
      std::size_t index_begin, index_end; ///< into nonlinear_constraint_indices
   };

   struct ModelData {
      NlHeader header;
      int number_variables{0};
      int number_constraints{0};
      int number_objectives{0};

      std::vector<double> variable_lower_bounds, variable_upper_bounds;
      std::vector<double> constraint_lower_bounds, constraint_upper_bounds;
      std::vector<double> initial_primal_point, initial_dual_point;
      std::vector<int> complementary_variables;
      std::vector<char> variable_is_integer;
      std::vector<char> objective_is_maximization;
      std::vector<Suffix> suffixes;

      /// objectives first, then constraints
      std::vector<FunctionBlock> functions;
      std::vector<int> quadratic_constraint_indices; ///< constraints with a quadratic part
      std::vector<int> nonlinear_constraint_indices; ///< constraints with elements
      std::vector<ConstraintRun> constraint_runs;
      std::vector<int> objective_linear_variables;
      std::vector<double> objective_linear_coefficients;

      // Jacobian (CSR)
      std::vector<std::size_t> jacobian_row_starts;
      /// COO row indices of the Jacobian, COO column indices of the Hessian: convenience copies of the compressed
      /// structures, built on first request (thread-safe)
      mutable std::vector<int> jacobian_row_indices;
      mutable std::once_flag jacobian_row_indices_built;
      std::vector<int> jacobian_column_indices;
      std::vector<double> jacobian_linear_coefficients;

      // quadratic parts: Hessians of 1/2 x^T H x, lower COO, all functions concatenated
      std::vector<int> quadratic_rows, quadratic_columns;
      std::vector<double> quadratic_values;
      /// position of the row/column variable in the function's gradient: the variable itself for objectives,
      /// the Jacobian nonzero index for constraints
      std::vector<std::uint32_t> quadratic_row_positions, quadratic_column_positions;
      std::vector<std::uint32_t> quadratic_hessian_slots;

      // non-quadratic elements
      std::vector<NonlinearElement> elements;
      TapeStorage tapes;
      std::vector<std::uint32_t> element_variable_positions; ///< same convention as quadratic_row_positions
      std::vector<std::uint32_t> element_hessian_slots;
      // compact copies of the fields read by the evaluation loops (the structs above are large: streaming through
      // them would waste memory bandwidth)
      std::vector<double> constraint_constants;
      std::vector<std::size_t> constraint_element_begin, constraint_element_end;
      std::vector<std::uint32_t> element_outputs; ///< last instruction of each element
      std::vector<double> element_scales;
      std::vector<ReverseEdge> reverse_edges;
      // shared defined variables: tapes, reverse program (all seeded at their outputs) and variable leaves
      // (position = slot in the defined gradients)
      TapeStorage defined_storage;
      std::vector<DefinedTape> defined_tapes;
      std::vector<ReverseEdge> defined_reverse_edges;
      std::vector<VariableOccurrence> defined_variable_occurrences;
      std::size_t number_defined_partials{0};
      std::vector<std::int64_t> shared_defined_outputs; ///< per defined variable: output instruction of its tape, or -1
      std::vector<DefinedOccurrence> defined_occurrences;
      std::vector<std::uint32_t> defined_occurrence_positions;
      /// per defined occurrence (index stored in the DefinedValue leaf's constant): its pairs (element local
      /// variable, index of the partial derivative in the defined gradients), aligned with defined_occurrence_positions
      std::vector<DefinedLeafLink> defined_leaf_links;
      std::vector<std::uint32_t> defined_pair_locals, defined_pair_gradients;
      std::vector<std::uint32_t> defined_hessian_slots;
      std::size_t number_defined_gradients{0};
      std::vector<VariableOccurrence> variable_occurrences;
      /// size of the partials array: 2 per instruction, 1 per tape argument, and the constant 1 at the end
      std::size_t number_partials{0};
      std::size_t total_tape_length{0};
      std::size_t maximum_element_tape_length{0};
      std::size_t maximum_function_tape_length{0};
      std::size_t maximum_element_variables{0};

      // Lagrangian Hessian (lower triangle, CSC)
      int hessian_objective_index{0};
      std::vector<std::size_t> hessian_column_starts;
      std::vector<int> hessian_row_indices;
      mutable std::vector<int> hessian_column_indices;
      mutable std::once_flag hessian_column_indices_built;

      [[nodiscard]] const FunctionBlock& objective_block(int objective_index) const {
         return this->functions[static_cast<std::size_t>(objective_index)];
      }
      [[nodiscard]] const FunctionBlock& constraint_block(int constraint_index) const {
         return this->functions[static_cast<std::size_t>(this->number_objectives + constraint_index)];
      }
      /// the tape of an element (may contain DefinedValue leaves)
      [[nodiscard]] TapeView tape_of(const NonlinearElement& element) const {
         return {this->tapes.instructions.data(), element.instruction_begin, element.instruction_begin + element.instruction_count,
            this->tapes.arguments.data(), this->tapes.constants.data()};
      }
      /// all the shared defined variables (adjacent in their storage)
      [[nodiscard]] TapeView defined_tapes_view() const {
         return {this->defined_storage.instructions.data(), 0, this->defined_storage.instructions.size(),
            this->defined_storage.arguments.data(), this->defined_storage.constants.data()};
      }
      [[nodiscard]] TapeView tape_of(const DefinedTape& defined) const {
         return {this->defined_storage.instructions.data(), defined.instruction_begin, defined.instruction_begin + defined.instruction_count,
            this->defined_storage.arguments.data(), this->defined_storage.constants.data()};
      }
      [[nodiscard]] TapeView tape_of(const ConstraintRun& run) const {
         return {this->tapes.instructions.data(), run.tape_begin, run.tape_end, this->tapes.arguments.data(),
            this->tapes.constants.data()};
      }
      /// all the elements of a function
      [[nodiscard]] TapeView tape_of(const FunctionBlock& block) const {
         return {this->tapes.instructions.data(), block.tape_begin, block.tape_end, this->tapes.arguments.data(),
            this->tapes.constants.data()};
      }
   };

} // namespace cppasl::detail
