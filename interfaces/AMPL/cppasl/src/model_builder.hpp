#pragma once

#include <memory>
#include <unordered_map>
#include <vector>
#include "cppasl/nl_model.hpp"
#include "expression_arena.hpp"
#include "function_decomposition.hpp"
#include "model_data.hpp"

namespace cppasl::detail {

   /// Receives the segments of an .nl file (in any order allowed by the format) and assembles the ModelData.
   /// Function expressions are processed as soon as they are read (decomposed, quadratic parts merged, elements
   /// compiled into tapes) and then discarded, which keeps the memory footprint proportional to the output.
   class ModelBuilder {
   public:
      ModelBuilder(const NlHeader& header, const ReaderOptions& options);

      ExpressionArena& expression_arena() { return this->arena; }
      ModelData& model() { return *this->data; }

      void define_variable(int defined_index, const LinearTerm* linear_terms, std::size_t number_linear_terms,
         NodeIndex expression);
      /// the expression (and everything parsed after `mark`) is discarded after processing
      void set_objective(int objective_index, bool is_maximization, NodeIndex root, const ExpressionArena::Mark& mark);
      void set_constraint_body(int constraint_index, NodeIndex root, const ExpressionArena::Mark& mark);
      void add_jacobian_entry(int constraint_index, int variable, double coefficient);
      void add_objective_gradient_entry(int objective_index, int variable, double coefficient);

      std::shared_ptr<ModelData> finalize();

   private:
      struct SparseEntry {
         int owner; ///< constraint, objective or function index
         int variable;
         double coefficient;
      };
      struct HessianRecord {
         int row;
         int column;
         std::uint32_t source;
      };

      std::shared_ptr<ModelData> data;
      const ReaderOptions options;
      ExpressionArena arena;
      ExpressionAnalysis analysis;
      FunctionDecomposer decomposer;
      TapeCompiler compiler;
      DecomposedFunction decomposed;
      std::vector<QuadraticTerm> quadratic_sort_buffer;
      std::vector<std::size_t> bucket_counts;
      std::vector<SparseEntry> jacobian_entries;
      std::vector<SparseEntry> objective_gradient_entries;
      std::vector<SparseEntry> expression_linear_entries; ///< linear terms found inside expressions (owner = function)

      void process_function(std::size_t function_index, NodeIndex root, const ExpressionArena::Mark& mark);
      /// merges the elements of the current function that use the same nonlinear defined variables, so that each
      /// defined variable is evaluated once per function
      void merge_elements_sharing_defined_variables();
      std::vector<std::uint32_t> defined_variable_stamp;
      std::uint32_t stamp{0};
      std::vector<NodeIndex> dfs_stack;
      void append_quadratic_part(FunctionBlock& block);
      void sort_local_variables(TapeStorage& storage, std::size_t instruction_begin, std::size_t instruction_count,
         std::size_t variable_begin, std::size_t variable_count);
      void prepare_shared_defined_variables();
      std::vector<char> element_uses_shared_defined;
      std::vector<char> defined_sharing_decided;
      std::vector<std::pair<std::size_t, std::size_t>> shared_defined_variable_ranges;
      std::unordered_map<std::size_t, std::size_t> defined_tape_index;
      void assemble_constraint_rows();
      void assemble_objectives();
      void assemble_hessian_structure();
      void compute_variable_integrality();
      void build_reverse_programs();
   };

} // namespace cppasl::detail
