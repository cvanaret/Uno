#include "model_builder.hpp"
#include <algorithm>
#include <limits>
#include <memory>
#include <numeric>
#include <memory>
#include <stdexcept>
#include <string>
#include <unordered_map>

namespace cppasl::detail {

   namespace {
      constexpr double infinity = std::numeric_limits<double>::infinity();

      std::uint32_t checked_position(std::size_t position) {
         if (position >= std::numeric_limits<std::uint32_t>::max()) {
            throw std::runtime_error("cppasl: more than 2^32-1 Jacobian nonzeros are not supported");
         }
         return static_cast<std::uint32_t>(position);
      }

      /// Stable counting sort of `entries` by `key(entry)` in [0, number_keys); `counts` and `output` are scratch.
      template <typename Entry, typename Key>
      void counting_sort(std::vector<Entry>& entries, std::vector<Entry>& output, std::vector<std::size_t>& counts,
            std::size_t number_keys, Key key) {
         counts.assign(number_keys + 1, 0);
         for (const Entry& entry: entries) ++counts[static_cast<std::size_t>(key(entry)) + 1];
         for (std::size_t k = 0; k < number_keys; ++k) counts[k + 1] += counts[k];
         output.resize(entries.size());
         for (const Entry& entry: entries) output[counts[static_cast<std::size_t>(key(entry))]++] = entry;
         entries.swap(output);
      }
   } // namespace

   ModelBuilder::ModelBuilder(const NlHeader& header, const ReaderOptions& options):
         data(std::make_shared<ModelData>()), options(options),
         analysis(this->arena, header.number_variables),
         decomposer(this->arena, this->analysis, header.number_variables, options.detect_quadratic_structure),
         compiler(this->arena, this->analysis, header.number_variables, this->data->tapes) {
      ModelData& model = *this->data;
      model.header = header;
      model.number_variables = header.number_variables;
      model.number_constraints = header.number_constraints;
      model.number_objectives = header.number_objectives;
      const auto n = static_cast<std::size_t>(header.number_variables);
      const auto m = static_cast<std::size_t>(header.number_constraints);
      model.variable_lower_bounds.assign(n, -infinity);
      model.variable_upper_bounds.assign(n, infinity);
      model.constraint_lower_bounds.assign(m, -infinity);
      model.constraint_upper_bounds.assign(m, infinity);
      model.initial_primal_point.assign(n, 0.);
      model.initial_dual_point.assign(m, 0.);
      model.complementary_variables.assign(m, -1);
      model.objective_is_maximization.assign(static_cast<std::size_t>(header.number_objectives), 0);
      model.functions.resize(static_cast<std::size_t>(header.number_objectives) + m);
      model.elements.reserve(static_cast<std::size_t>(header.number_nonlinear_constraints + header.number_nonlinear_objectives));
      this->arena.defined_variables.resize(static_cast<std::size_t>(header.number_defined_variables()));
      // the default (0) silently becomes -1 for feasibility problems without objective
      const int hessian_objective_index = (header.number_objectives == 0 && options.hessian_objective_index == 0) ? -1 :
         options.hessian_objective_index;
      if (hessian_objective_index >= header.number_objectives) {
         throw std::invalid_argument("cppasl: hessian_objective_index is out of range");
      }
      model.hessian_objective_index = hessian_objective_index;
   }

   void ModelBuilder::define_variable(int defined_index, const LinearTerm* linear_terms, std::size_t number_linear_terms,
         NodeIndex expression) {
      DefinedVariable& defined = this->arena.defined_variables.at(static_cast<std::size_t>(defined_index));
      defined.linear_begin = static_cast<std::uint32_t>(this->arena.defined_variable_linear_terms.size());
      this->arena.defined_variable_linear_terms.insert(this->arena.defined_variable_linear_terms.end(), linear_terms,
         linear_terms + number_linear_terms);
      defined.linear_end = static_cast<std::uint32_t>(this->arena.defined_variable_linear_terms.size());
      defined.expression = expression;
      defined.is_defined = true;
   }

   void ModelBuilder::set_objective(int objective_index, bool is_maximization, NodeIndex root,
         const ExpressionArena::Mark& mark) {
      this->data->objective_is_maximization[static_cast<std::size_t>(objective_index)] = is_maximization ? 1 : 0;
      this->process_function(static_cast<std::size_t>(objective_index), root, mark);
   }

   void ModelBuilder::set_constraint_body(int constraint_index, NodeIndex root, const ExpressionArena::Mark& mark) {
      this->process_function(static_cast<std::size_t>(this->data->number_objectives + constraint_index), root, mark);
   }

   void ModelBuilder::add_jacobian_entry(int constraint_index, int variable, double coefficient) {
      this->jacobian_entries.push_back({constraint_index, variable, coefficient});
   }

   void ModelBuilder::add_objective_gradient_entry(int objective_index, int variable, double coefficient) {
      this->objective_gradient_entries.push_back({objective_index, variable, coefficient});
   }

   // ------------------------------------------------------------------------------------------------------------------
   // per-function processing (streaming)
   // ------------------------------------------------------------------------------------------------------------------

   void ModelBuilder::process_function(std::size_t function_index, NodeIndex root, const ExpressionArena::Mark& mark) {
      ModelData& model = *this->data;
      FunctionBlock& block = model.functions[function_index];
      this->decomposer.decompose(root, this->decomposed);

      block.constant = this->decomposed.constant;
      for (const LinearTerm& term: this->decomposed.linear_terms) {
         this->expression_linear_entries.push_back({static_cast<int>(function_index), term.variable, term.coefficient});
      }
      this->append_quadratic_part(block);
      if (this->decomposed.elements.size() > 1 && !this->arena.defined_variables.empty()) {
         this->merge_elements_sharing_defined_variables();
      }

      // shared defined variables used by the elements: their tapes go before the function's block
      this->prepare_shared_defined_variables();

      block.element_begin = model.elements.size();
      block.tape_begin = model.tapes.instructions.size();
      for (std::size_t e = 0; e < this->decomposed.elements.size(); ++e) {
         const DecomposedFunction::Element& element = this->decomposed.elements[e];
         const bool uses_shared = this->element_uses_shared_defined[e] != 0;
         const CompiledTape tape = this->compiler.compile(element.root, model.tapes, uses_shared, true);
         auto u32 = [](std::size_t value) { return static_cast<std::uint32_t>(value); };
         NonlinearElement compiled{u32(function_index), u32(tape.instruction_begin), u32(tape.instruction_count),
            u32(tape.variable_begin), u32(tape.variable_count), 0, element.scale};
         compiled.has_defined_leaves = uses_shared;
         this->sort_local_variables(model.tapes, compiled.instruction_begin, compiled.instruction_count, compiled.variable_begin,
            compiled.variable_count);
         // the group structure pays off when the element has more variables than one vector Hessian sweep handles
         // (not with DefinedValue leaves: their contributions go through the chain rule)
         if (!uses_shared && tape.group_sum >= 0 && tape.variable_count > maximum_hessian_directions) {
            compiled.group_sum = u32(static_cast<std::size_t>(tape.group_sum));
            compiled.group_range_begin = u32(tape.group_range_begin);
            compiled.group_range_end = u32(tape.group_range_end);
         }
         model.elements.push_back(compiled);
      }
      block.element_end = model.elements.size();
      block.tape_end = model.tapes.instructions.size();

      // the expression is no longer needed: only the defined variables (parsed before `mark`) stay resident
      this->arena.truncate(mark);
      this->analysis.forget_nodes_from(mark.node_count);
   }

   void ModelBuilder::merge_elements_sharing_defined_variables() {
      constexpr std::size_t maximum_shared_defined_variables = 4; // bounds the size of merged elements
      const int n = this->data->number_variables;
      std::vector<DecomposedFunction::Element>& elements = this->decomposed.elements;
      if (this->defined_variable_stamp.size() < this->arena.defined_variables.size()) {
         this->defined_variable_stamp.resize(this->arena.defined_variables.size(), 0);
      }
      // defined variables (transitively) used by each element
      std::vector<std::vector<int>> used(elements.size());
      for (std::size_t e = 0; e < elements.size(); ++e) {
         ++this->stamp;
         this->dfs_stack.assign(1, elements[e].root);
         while (!this->dfs_stack.empty()) {
            const ExpressionNode& node = this->arena.nodes[this->dfs_stack.back()];
            this->dfs_stack.pop_back();
            if (node.opcode == variable_opcode) {
               if (static_cast<int>(node.first_argument) < n) continue;
               const std::size_t defined = node.first_argument - static_cast<std::size_t>(n);
               if (this->defined_variable_stamp[defined] == this->stamp) continue;
               this->defined_variable_stamp[defined] = this->stamp;
               used[e].push_back(static_cast<int>(defined));
               this->dfs_stack.push_back(this->arena.defined_variables[defined].expression);
               continue;
            }
            for (std::uint32_t k = 0; k < node.argument_count; ++k) this->dfs_stack.push_back(this->arena.operand(node, k));
         }
         std::sort(used[e].begin(), used[e].end());
      }
      // greedy grouping: an element joins the groups it shares defined variables with, if the union stays small
      std::vector<int> group_of_element(elements.size(), -1);
      std::vector<std::vector<int>> group_defined(elements.size());
      std::unordered_map<int, int> group_of_defined;
      int number_groups = 0;
      std::vector<int> candidates, merged;
      for (std::size_t e = 0; e < elements.size(); ++e) {
         candidates.clear();
         for (int defined: used[e]) {
            const auto found = group_of_defined.find(defined);
            if (found != group_of_defined.end()) candidates.push_back(found->second);
         }
         std::sort(candidates.begin(), candidates.end());
         candidates.erase(std::unique(candidates.begin(), candidates.end()), candidates.end());
         merged = used[e];
         for (int g: candidates) merged.insert(merged.end(), group_defined[static_cast<std::size_t>(g)].begin(), group_defined[static_cast<std::size_t>(g)].end());
         std::sort(merged.begin(), merged.end());
         merged.erase(std::unique(merged.begin(), merged.end()), merged.end());
         int group;
         if (!candidates.empty() && merged.size() <= maximum_shared_defined_variables) {
            group = candidates.front();
            for (std::size_t k = 1; k < candidates.size(); ++k) { // absorb the other groups
               for (int& g: group_of_element) if (g == candidates[k]) g = group;
               group_defined[static_cast<std::size_t>(candidates[k])].clear();
            }
            group_defined[static_cast<std::size_t>(group)] = merged;
         }
         else {
            group = number_groups++;
            group_defined[static_cast<std::size_t>(group)] = used[e];
         }
         group_of_element[e] = group;
         for (int defined: group_defined[static_cast<std::size_t>(group)]) group_of_defined[defined] = group;
      }
      // one element per group: sum of the scaled roots
      std::vector<std::vector<std::size_t>> members(static_cast<std::size_t>(number_groups));
      for (std::size_t e = 0; e < elements.size(); ++e) members[static_cast<std::size_t>(group_of_element[e])].push_back(e);
      std::vector<DecomposedFunction::Element> result;
      std::vector<NodeIndex> operands;
      for (const std::vector<std::size_t>& group: members) {
         if (group.empty()) continue;
         if (group.size() == 1) { result.push_back(elements[group.front()]); continue; }
         operands.clear();
         for (std::size_t e: group) {
            NodeIndex root = elements[e].root;
            if (elements[e].scale != 1.) {
               const NodeIndex factors[2] = {this->arena.add_constant(elements[e].scale), root};
               root = this->arena.add_operation(nl_opcode::multiply, factors, 2);
            }
            operands.push_back(root);
         }
         result.push_back({this->arena.add_operation(nl_opcode::sum_list, operands.data(), static_cast<std::uint32_t>(operands.size())), 1.});
      }
      elements.swap(result);
   }

   /// Renumbers the local variables of a compiled tape in increasing global order (the Hessian structure relies on
   /// it: the rows of local column b are the local variables b, b + 1, ...)
   void ModelBuilder::sort_local_variables(TapeStorage& storage, std::size_t instruction_begin, std::size_t instruction_count,
         std::size_t variable_begin, std::size_t variable_count) {
      int* variables = this->data->tapes.element_variables.data() + variable_begin;
      if (std::is_sorted(variables, variables + variable_count)) return;
      std::vector<std::size_t> order(variable_count);
      std::iota(order.begin(), order.end(), 0);
      std::sort(order.begin(), order.end(), [&](std::size_t a, std::size_t b) { return variables[a] < variables[b]; });
      std::vector<std::int32_t> new_local(variable_count);
      std::vector<int> sorted_variables(variable_count);
      for (std::size_t k = 0; k < order.size(); ++k) {
         new_local[order[k]] = static_cast<std::int32_t>(k);
         sorted_variables[k] = variables[order[k]];
      }
      std::copy(sorted_variables.begin(), sorted_variables.end(), variables);
      for (std::size_t i = instruction_begin; i < instruction_begin + instruction_count; ++i) {
         TapeInstruction& instruction = storage.instructions[i];
         if (instruction.operation == TapeOperation::Variable) instruction.second = new_local[static_cast<std::size_t>(instruction.second)];
      }
   }

   /// Decides which nonlinear defined variables used by the current elements are evaluated once per point by their
   /// own tape (shared by all functions), compiles the new ones, and flags the elements that read them.
   void ModelBuilder::prepare_shared_defined_variables() {
      ModelData& model = *this->data;
      const int n = model.number_variables;
      const std::size_t number_defined = this->arena.defined_variables.size();
      this->element_uses_shared_defined.assign(this->decomposed.elements.size(), 0);
      if (number_defined == 0) return;
      if (model.shared_defined_outputs.size() < number_defined) {
         model.shared_defined_outputs.resize(number_defined, -1);
         this->defined_sharing_decided.resize(number_defined, 0);
         this->compiler.shared_defined_outputs = &model.shared_defined_outputs;
         this->shared_defined_variable_ranges.resize(number_defined, {0, 0});
         this->compiler.shared_defined_variable_ranges = &this->shared_defined_variable_ranges;
      }
      // a defined variable is shared if AMPL reports it used by several functions (the c1/o1 ones, used by a single
      // function, come last) and its expression is not trivial (otherwise inlining is cheaper)
      constexpr std::size_t minimum_shared_size = 4;
      const NlHeader& h = model.header;
      const auto number_multiple_use = static_cast<std::size_t>(h.number_common_expressions_in_both +
         h.number_common_expressions_in_constraints + h.number_common_expressions_in_objectives);
      auto expression_size = [&](std::size_t defined) {
         const DefinedVariable& variable = this->arena.defined_variables[defined];
         std::size_t size = variable.linear_end - variable.linear_begin;
         this->dfs_stack.assign(1, variable.expression);
         while (!this->dfs_stack.empty() && size < minimum_shared_size) {
            const ExpressionNode& node = this->arena.nodes[this->dfs_stack.back()];
            this->dfs_stack.pop_back();
            ++size;
            for (std::uint32_t k = 0; k < node.argument_count; ++k) this->dfs_stack.push_back(this->arena.operand(node, k));
         }
         return size;
      };
      std::vector<NodeIndex> stack;
      for (std::size_t e = 0; e < this->decomposed.elements.size(); ++e) {
         stack.assign(1, this->decomposed.elements[e].root);
         while (!stack.empty()) {
            const ExpressionNode& node = this->arena.nodes[stack.back()];
            stack.pop_back();
            if (node.opcode == variable_opcode) {
               if (static_cast<int>(node.first_argument) < n) continue;
               const std::size_t defined = node.first_argument - static_cast<std::size_t>(n);
               if (!this->defined_sharing_decided[defined]) {
                  this->defined_sharing_decided[defined] = 1;
                  if (defined < number_multiple_use && expression_size(defined) >= minimum_shared_size) {
                     const NodeIndex reference = this->arena.add_variable(n + static_cast<int>(defined));
                     const CompiledTape tape = this->compiler.compile(reference, model.defined_storage, false, true);
                     this->sort_local_variables(model.defined_storage, tape.instruction_begin, tape.instruction_count,
                        tape.variable_begin, tape.variable_count);
                     model.defined_tapes.push_back({tape.instruction_begin, tape.instruction_count, tape.variable_begin,
                        tape.variable_count, model.number_defined_gradients});
                     model.number_defined_gradients += tape.variable_count;
                     model.shared_defined_outputs[defined] = static_cast<std::int64_t>(tape.instruction_begin + tape.instruction_count - 1);
                     this->shared_defined_variable_ranges[defined] = {tape.variable_begin, tape.variable_begin + tape.variable_count};
                     this->defined_tape_index[defined] = model.defined_tapes.size() - 1;
                  }
               }
               if (model.shared_defined_outputs[defined] >= 0) {
                  this->element_uses_shared_defined[e] = 1;
                  continue; // a leaf of the evaluation tape
               }
               stack.push_back(this->arena.defined_variables[defined].expression);
               continue;
            }
            for (std::uint32_t k = 0; k < node.argument_count; ++k) stack.push_back(this->arena.operand(node, k));
         }
      }
   }

   /// Merges the duplicate quadratic terms of the current function and appends the Hessian of its quadratic part:
   /// diagonal entries first, then the strictly lower entries sorted by (column, row). Exact zeros are dropped.
   void ModelBuilder::append_quadratic_part(FunctionBlock& block) {
      ModelData& model = *this->data;
      std::vector<QuadraticTerm>& terms = this->decomposed.quadratic_terms;
      block.quadratic_begin = block.quadratic_off_diagonal_begin = block.quadratic_end = model.quadratic_values.size();
      if (terms.empty()) return;

      // sort by (column, row): LSD radix (two stable counting sorts) for large functions, comparison sort otherwise
      const auto n = static_cast<std::size_t>(model.number_variables);
      if (terms.size() > 64 && terms.size() >= n / 4) {
         counting_sort(terms, this->quadratic_sort_buffer, this->bucket_counts, n, [](const QuadraticTerm& t) { return t.row; });
         counting_sort(terms, this->quadratic_sort_buffer, this->bucket_counts, n, [](const QuadraticTerm& t) { return t.column; });
      }
      else {
         std::sort(terms.begin(), terms.end(), [](const QuadraticTerm& a, const QuadraticTerm& b) {
            return (a.column != b.column) ? a.column < b.column : a.row < b.row;
         });
      }
      // merge duplicates in place
      std::size_t merged_count = 0;
      for (std::size_t k = 0; k < terms.size();) {
         QuadraticTerm merged = terms[k++];
         while (k < terms.size() && terms[k].row == merged.row && terms[k].column == merged.column) {
            merged.coefficient += terms[k++].coefficient;
         }
         if (merged.coefficient != 0.) terms[merged_count++] = merged;
      }
      terms.resize(merged_count);

      // H_ii = 2 c for c x_i^2, H_ij = c for c x_i x_j (i != j), so that the part equals 1/2 x^T H x
      for (const QuadraticTerm& term: terms) {
         if (term.row == term.column) {
            model.quadratic_rows.push_back(term.row);
            model.quadratic_columns.push_back(term.column);
            model.quadratic_values.push_back(2. * term.coefficient);
         }
      }
      block.quadratic_off_diagonal_begin = model.quadratic_values.size();
      for (const QuadraticTerm& term: terms) {
         if (term.row != term.column) {
            model.quadratic_rows.push_back(term.row);
            model.quadratic_columns.push_back(term.column);
            model.quadratic_values.push_back(term.coefficient);
         }
      }
      block.quadratic_end = model.quadratic_values.size();
   }

   // ------------------------------------------------------------------------------------------------------------------
   // final assembly
   // ------------------------------------------------------------------------------------------------------------------

   std::shared_ptr<ModelData> ModelBuilder::finalize() {
      ModelData& model = *this->data;
      model.quadratic_row_positions.resize(model.quadratic_values.size());
      model.quadratic_column_positions.resize(model.quadratic_values.size());
      model.element_variable_positions.resize(model.tapes.element_variables.size());

      // group the linear entries by owner (stable: keeps the file order within an owner)
      std::vector<SparseEntry> buffer;
      counting_sort(this->jacobian_entries, buffer, this->bucket_counts, static_cast<std::size_t>(model.number_constraints),
         [](const SparseEntry& e) { return e.owner; });
      counting_sort(this->objective_gradient_entries, buffer, this->bucket_counts,
         static_cast<std::size_t>(model.number_objectives), [](const SparseEntry& e) { return e.owner; });
      counting_sort(this->expression_linear_entries, buffer, this->bucket_counts, model.functions.size(),
         [](const SparseEntry& e) { return e.owner; });

      this->assemble_constraint_rows();
      this->assemble_objectives();
      this->assemble_hessian_structure();
      this->compute_variable_integrality();
      this->build_reverse_programs();
      for (const NonlinearElement& element: model.elements) {
         model.element_outputs.push_back(static_cast<std::uint32_t>(element.instruction_begin + element.instruction_count - 1));
         model.element_scales.push_back(element.scale);
      }
      for (int constraint = 0; constraint < model.number_constraints; ++constraint) {
         const FunctionBlock& block = model.constraint_block(constraint);
         model.constraint_constants.push_back(block.constant);
         model.constraint_element_begin.push_back(block.element_begin);
         model.constraint_element_end.push_back(block.element_end);
      }
      // workspace sizes
      model.total_tape_length = model.tapes.instructions.size();
      for (const FunctionBlock& block: model.functions) {
         model.maximum_function_tape_length = std::max(model.maximum_function_tape_length, block.tape_end - block.tape_begin);
      }
      for (const NonlinearElement& element: model.elements) {
         model.maximum_element_tape_length = std::max(model.maximum_element_tape_length, static_cast<std::size_t>(element.instruction_count));
         model.maximum_element_variables = std::max(model.maximum_element_variables, static_cast<std::size_t>(element.variable_count));
      }
      for (const DefinedTape& defined: model.defined_tapes) { // their Hessians use the same scratch
         model.maximum_element_tape_length = std::max(model.maximum_element_tape_length, defined.instruction_count);
         model.maximum_element_variables = std::max(model.maximum_element_variables, defined.variable_count);
      }
      for (int constraint = 0; constraint < model.number_constraints; ++constraint) {
         const FunctionBlock& block = model.constraint_block(constraint);
         if (block.has_quadratic_part()) model.quadratic_constraint_indices.push_back(constraint);
         if (block.has_elements()) model.nonlinear_constraint_indices.push_back(constraint);
      }
      // runs of nonlinear constraints whose tapes are adjacent: swept in one loop by the constraint/Jacobian evaluations
      for (std::size_t k = 0; k < model.nonlinear_constraint_indices.size(); ++k) {
         const FunctionBlock& block = model.constraint_block(model.nonlinear_constraint_indices[k]);
         if (!model.constraint_runs.empty() && model.constraint_runs.back().tape_end == block.tape_begin) {
            model.constraint_runs.back().tape_end = block.tape_end;
            model.constraint_runs.back().index_end = k + 1;
         }
         else model.constraint_runs.push_back({block.tape_begin, block.tape_end, k, k + 1});
      }
      // release the builder's scratch memory
      std::vector<SparseEntry>().swap(this->jacobian_entries);
      std::vector<SparseEntry>().swap(this->objective_gradient_entries);
      std::vector<SparseEntry>().swap(this->expression_linear_entries);
      return std::move(this->data);
   }

   /// Builds the CSR Jacobian pattern (union of the variables of the linear, quadratic and element parts of each
   /// constraint), the linear coefficients aligned with it, and the Jacobian positions of the quadratic and element
   /// variables.
   void ModelBuilder::assemble_constraint_rows() {
      ModelData& model = *this->data;
      const auto n = static_cast<std::size_t>(model.number_variables);
      const auto m = static_cast<std::size_t>(model.number_constraints);
      constexpr std::size_t unset = std::numeric_limits<std::size_t>::max();
      std::vector<std::size_t> position_of_variable(n, unset); // Jacobian position of each variable in the current row
      std::vector<int> row_variables;

      model.jacobian_row_starts.assign(m + 1, 0);
      model.jacobian_column_indices.reserve(this->jacobian_entries.size());
      std::size_t jacobian_cursor = 0, expression_cursor = 0;
      // skip the expression linear terms of the objectives (they come first after grouping)
      while (expression_cursor < this->expression_linear_entries.size() &&
            this->expression_linear_entries[expression_cursor].owner < model.number_objectives) {
         ++expression_cursor;
      }

      for (std::size_t constraint = 0; constraint < m; ++constraint) {
         const int function_index = model.number_objectives + static_cast<int>(constraint);
         FunctionBlock& block = model.functions[static_cast<std::size_t>(function_index)];
         row_variables.clear();
         auto touch = [&](int variable) {
            std::size_t& position = position_of_variable[static_cast<std::size_t>(variable)];
            if (position == unset) {
               position = 0;
               row_variables.push_back(variable);
            }
         };
         const std::size_t jacobian_begin = jacobian_cursor;
         while (jacobian_cursor < this->jacobian_entries.size() &&
               this->jacobian_entries[jacobian_cursor].owner == static_cast<int>(constraint)) {
            touch(this->jacobian_entries[jacobian_cursor++].variable);
         }
         const std::size_t expression_begin = expression_cursor;
         while (expression_cursor < this->expression_linear_entries.size() &&
               this->expression_linear_entries[expression_cursor].owner == function_index) {
            touch(this->expression_linear_entries[expression_cursor++].variable);
         }
         for (std::size_t q = block.quadratic_begin; q < block.quadratic_end; ++q) {
            touch(model.quadratic_rows[q]);
            touch(model.quadratic_columns[q]);
         }
         for (std::size_t e = block.element_begin; e < block.element_end; ++e) {
            const NonlinearElement& element = model.elements[e];
            for (std::size_t v = 0; v < element.variable_count; ++v) touch(model.tapes.element_variables[element.variable_begin + v]);
         }

         // row pattern, sorted by column
         std::sort(row_variables.begin(), row_variables.end());
         const std::size_t row_begin = model.jacobian_column_indices.size();
         for (int variable: row_variables) {
            position_of_variable[static_cast<std::size_t>(variable)] = model.jacobian_column_indices.size();
            model.jacobian_column_indices.push_back(variable);
         }
         model.jacobian_row_starts[constraint + 1] = model.jacobian_column_indices.size();
         block.linear_begin = row_begin;
         block.linear_end = model.jacobian_column_indices.size();

         // linear coefficients aligned with the pattern
         model.jacobian_linear_coefficients.resize(model.jacobian_column_indices.size(), 0.);
         for (std::size_t k = jacobian_begin; k < jacobian_cursor; ++k) {
            const SparseEntry& entry = this->jacobian_entries[k];
            model.jacobian_linear_coefficients[position_of_variable[static_cast<std::size_t>(entry.variable)]] += entry.coefficient;
         }
         for (std::size_t k = expression_begin; k < expression_cursor; ++k) {
            const SparseEntry& entry = this->expression_linear_entries[k];
            model.jacobian_linear_coefficients[position_of_variable[static_cast<std::size_t>(entry.variable)]] += entry.coefficient;
         }
         // positions of the nonlinear parts
         for (std::size_t q = block.quadratic_begin; q < block.quadratic_end; ++q) {
            model.quadratic_row_positions[q] = checked_position(position_of_variable[static_cast<std::size_t>(model.quadratic_rows[q])]);
            model.quadratic_column_positions[q] = checked_position(position_of_variable[static_cast<std::size_t>(model.quadratic_columns[q])]);
         }
         for (std::size_t e = block.element_begin; e < block.element_end; ++e) {
            const NonlinearElement& element = model.elements[e];
            for (std::size_t v = element.variable_begin; v < element.variable_begin + element.variable_count; ++v) {
               model.element_variable_positions[v] = checked_position(position_of_variable[static_cast<std::size_t>(model.tapes.element_variables[v])]);
            }
         }
         for (int variable: row_variables) position_of_variable[static_cast<std::size_t>(variable)] = unset;
      }


   }

   /// Objectives: merged sparse linear part; the gradient positions of the nonlinear parts are the variables.
   void ModelBuilder::assemble_objectives() {
      ModelData& model = *this->data;
      constexpr std::size_t unset = std::numeric_limits<std::size_t>::max();
      std::vector<std::size_t> position_of_variable(static_cast<std::size_t>(model.number_variables), unset);
      std::size_t gradient_cursor = 0, expression_cursor = 0;
      for (int objective = 0; objective < model.number_objectives; ++objective) {
         FunctionBlock& block = model.functions[static_cast<std::size_t>(objective)];
         block.linear_begin = model.objective_linear_variables.size();
         auto add = [&](const SparseEntry& entry) {
            std::size_t& position = position_of_variable[static_cast<std::size_t>(entry.variable)];
            if (position == unset) {
               position = model.objective_linear_variables.size();
               model.objective_linear_variables.push_back(entry.variable);
               model.objective_linear_coefficients.push_back(0.);
            }
            model.objective_linear_coefficients[position] += entry.coefficient;
         };
         while (gradient_cursor < this->objective_gradient_entries.size() &&
               this->objective_gradient_entries[gradient_cursor].owner == objective) {
            add(this->objective_gradient_entries[gradient_cursor++]);
         }
         while (expression_cursor < this->expression_linear_entries.size() &&
               this->expression_linear_entries[expression_cursor].owner == objective) {
            add(this->expression_linear_entries[expression_cursor++]);
         }
         block.linear_end = model.objective_linear_variables.size();
         for (std::size_t k = block.linear_begin; k < block.linear_end; ++k) {
            position_of_variable[static_cast<std::size_t>(model.objective_linear_variables[k])] = unset;
         }
         for (std::size_t q = block.quadratic_begin; q < block.quadratic_end; ++q) {
            model.quadratic_row_positions[q] = static_cast<std::uint32_t>(model.quadratic_rows[q]);
            model.quadratic_column_positions[q] = static_cast<std::uint32_t>(model.quadratic_columns[q]);
         }
         for (std::size_t e = block.element_begin; e < block.element_end; ++e) {
            const NonlinearElement& element = model.elements[e];
            for (std::size_t v = element.variable_begin; v < element.variable_begin + element.variable_count; ++v) {
               model.element_variable_positions[v] = static_cast<std::uint32_t>(model.tapes.element_variables[v]);
            }
         }
      }
   }

   /// Lower-triangular CSC pattern of the Lagrangian Hessian: union of the quadratic entries and the element variable
   /// pairs of the Hessian objective and of all constraints. Every contribution receives the index of its slot.
   void ModelBuilder::assemble_hessian_structure() {
      ModelData& model = *this->data;
      const auto n = static_cast<std::size_t>(model.number_variables);
      model.quadratic_hessian_slots.assign(model.quadratic_values.size(), no_hessian_slot);
      std::size_t number_element_slots = 0;
      for (NonlinearElement& element: model.elements) {
         element.hessian_slot_begin = static_cast<std::uint32_t>(number_element_slots);
         number_element_slots += element.variable_count * (element.variable_count + 1) / 2;
      }
      model.element_hessian_slots.assign(number_element_slots, no_hessian_slot);
      auto contributes = [&](std::size_t function_index) {
         return function_index >= static_cast<std::size_t>(model.number_objectives) ||
            static_cast<int>(function_index) == model.hessian_objective_index;
      };

      // Contributions bucketed by column: a quadratic entry is one record; the column b of an element is one segment
      // (its local variables are sorted, so the rows of column b are the variables b, b + 1, ... and its slots are
      // contiguous in the column-major packed storage). The records are thus proportional to the number of entries
      // of the quadratic forms plus the number of element variables, not to the size of the element Hessians.
      struct QuadraticRecord {
         int row;
         std::uint32_t entry;
      };
      struct SegmentRecord {
         std::uint32_t element;
         std::uint32_t local_column;
      };
      if (model.elements.size() + model.defined_tapes.size() >= no_hessian_slot) throw std::runtime_error("cppasl: too many elements");
      // shared defined variables read by the functions of the Lagrangian: their Hessians are contributions too
      // (a segment record with an index past the elements designates a defined tape)
      std::size_t number_defined_slots = 0;
      for (std::size_t f = 0; f < model.functions.size(); ++f) {
         if (!contributes(f)) continue;
         for (std::size_t e = model.functions[f].element_begin; e < model.functions[f].element_end; ++e) {
            const NonlinearElement& element = model.elements[e];
            if (!element.has_defined_leaves) continue;
            for (std::size_t i = element.instruction_begin; i < element.instruction_begin + element.instruction_count; ++i) {
               const TapeInstruction& instruction = model.tapes.instructions[i];
               if (instruction.operation == TapeOperation::DefinedValue) {
                  model.defined_tapes[this->defined_tape_index.at(static_cast<std::size_t>(instruction.second))].contributes_to_hessian = true;
               }
            }
         }
      }
      for (DefinedTape& defined: model.defined_tapes) {
         defined.hessian_slot_begin = number_defined_slots;
         if (defined.contributes_to_hessian) number_defined_slots += defined.variable_count * (defined.variable_count + 1) / 2;
      }
      model.defined_hessian_slots.assign(number_defined_slots, no_hessian_slot);
      std::vector<std::size_t> quadratic_starts(n + 1, 0), segment_starts(n + 1, 0);
      for (std::size_t f = 0; f < model.functions.size(); ++f) {
         if (!contributes(f)) continue;
         const FunctionBlock& block = model.functions[f];
         for (std::size_t q = block.quadratic_begin; q < block.quadratic_end; ++q) ++quadratic_starts[static_cast<std::size_t>(model.quadratic_columns[q]) + 1];
         for (std::size_t e = block.element_begin; e < block.element_end; ++e) {
            const NonlinearElement& element = model.elements[e];
            for (std::size_t v = 0; v < element.variable_count; ++v) {
               ++segment_starts[static_cast<std::size_t>(model.tapes.element_variables[element.variable_begin + v]) + 1];
            }
         }
      }
      for (const DefinedTape& defined: model.defined_tapes) {
         if (!defined.contributes_to_hessian) continue;
         for (std::size_t v = 0; v < defined.variable_count; ++v) {
            ++segment_starts[static_cast<std::size_t>(model.tapes.element_variables[defined.variable_begin + v]) + 1];
         }
      }
      for (std::size_t j = 0; j < n; ++j) {
         quadratic_starts[j + 1] += quadratic_starts[j];
         segment_starts[j + 1] += segment_starts[j];
      }
      const std::unique_ptr<QuadraticRecord[]> quadratic_records(new QuadraticRecord[quadratic_starts[n]]);
      const std::unique_ptr<SegmentRecord[]> segment_records(new SegmentRecord[segment_starts[n]]);
      {
         std::vector<std::size_t> quadratic_cursor(quadratic_starts.begin(), quadratic_starts.end() - 1);
         std::vector<std::size_t> segment_cursor(segment_starts.begin(), segment_starts.end() - 1);
         for (std::size_t f = 0; f < model.functions.size(); ++f) {
            if (!contributes(f)) continue;
            const FunctionBlock& block = model.functions[f];
            for (std::size_t q = block.quadratic_begin; q < block.quadratic_end; ++q) {
               quadratic_records[quadratic_cursor[static_cast<std::size_t>(model.quadratic_columns[q])]++] = {model.quadratic_rows[q],
                  static_cast<std::uint32_t>(q)};
            }
            for (std::size_t e = block.element_begin; e < block.element_end; ++e) {
               const NonlinearElement& element = model.elements[e];
               for (std::size_t v = 0; v < element.variable_count; ++v) {
                  const auto column = static_cast<std::size_t>(model.tapes.element_variables[element.variable_begin + v]);
                  segment_records[segment_cursor[column]++] = {static_cast<std::uint32_t>(e), static_cast<std::uint32_t>(v)};
               }
            }
         }
         for (std::size_t t = 0; t < model.defined_tapes.size(); ++t) {
            const DefinedTape& defined = model.defined_tapes[t];
            if (!defined.contributes_to_hessian) continue;
            for (std::size_t v = 0; v < defined.variable_count; ++v) {
               const auto column = static_cast<std::size_t>(model.tapes.element_variables[defined.variable_begin + v]);
               segment_records[segment_cursor[column]++] = {static_cast<std::uint32_t>(model.elements.size() + t), static_cast<std::uint32_t>(v)};
            }
         }
      }

      // for each column: distinct rows (marker array), sorted (or scanned when the column is dense), then slots
      model.hessian_column_starts.assign(n + 1, 0);
      model.hessian_row_indices.clear();
      // exact upper bound (only the touched pages are committed): no reallocation copies
      model.hessian_row_indices.reserve(quadratic_starts[n] + [&] {
         std::size_t bound = 0;
         for (std::size_t r = 0; r < segment_starts[n]; ++r) {
            const std::size_t index = segment_records[r].element, b = segment_records[r].local_column;
            bound += ((index < model.elements.size()) ? model.elements[index].variable_count
               : model.defined_tapes[index - model.elements.size()].variable_count) - b;
         }
         return bound;
      }());
      std::vector<std::uint32_t> local_slot_of_row(n, no_hessian_slot);
      std::vector<int> column_rows, sorted_rows;
      column_rows.reserve(n);
      sorted_rows.reserve(n);
      std::vector<std::uint32_t> final_slot_of_local;
      auto for_each_row = [&](std::size_t j, auto&& visit) { // visit(row, slot reference)
         for (std::size_t r = quadratic_starts[j]; r < quadratic_starts[j + 1]; ++r) {
            visit(quadratic_records[r].row, model.quadratic_hessian_slots[quadratic_records[r].entry]);
         }
         for (std::size_t r = segment_starts[j]; r < segment_starts[j + 1]; ++r) {
            const std::size_t index = segment_records[r].element, b = segment_records[r].local_column;
            std::size_t k, variable_begin;
            std::uint32_t* slots;
            if (index < model.elements.size()) {
               const NonlinearElement& element = model.elements[index];
               k = element.variable_count;
               variable_begin = element.variable_begin;
               slots = model.element_hessian_slots.data() + element.hessian_slot_begin;
            }
            else {
               const DefinedTape& defined = model.defined_tapes[index - model.elements.size()];
               k = defined.variable_count;
               variable_begin = defined.variable_begin;
               slots = model.defined_hessian_slots.data() + defined.hessian_slot_begin;
            }
            const int* variables = model.tapes.element_variables.data() + variable_begin;
            slots += packed_lower_index(b, b, k);
            for (std::size_t a = b; a < k; ++a) visit(variables[a], slots[a - b]);
         }
      };
      for (std::size_t j = 0; j < n; ++j) {
         column_rows.clear();
         for_each_row(j, [&](int row, std::uint32_t&) {
            std::uint32_t& local = local_slot_of_row[static_cast<std::size_t>(row)];
            if (local == no_hessian_slot) {
               local = static_cast<std::uint32_t>(column_rows.size());
               column_rows.push_back(row);
            }
         });
         sorted_rows.clear();
         if (16 * column_rows.size() >= n - j) {
            for (std::size_t row = j; row < n; ++row) if (local_slot_of_row[row] != no_hessian_slot) sorted_rows.push_back(static_cast<int>(row));
         }
         else {
            sorted_rows.assign(column_rows.begin(), column_rows.end());
            std::sort(sorted_rows.begin(), sorted_rows.end());
         }
         const std::size_t column_begin = model.hessian_row_indices.size();
         if (column_begin + sorted_rows.size() >= no_hessian_slot) {
            throw std::runtime_error("cppasl: more than 2^32-1 Hessian nonzeros are not supported");
         }
         final_slot_of_local.resize(column_rows.size());
         for (std::size_t k = 0; k < sorted_rows.size(); ++k) {
            final_slot_of_local[local_slot_of_row[static_cast<std::size_t>(sorted_rows[k])]] = static_cast<std::uint32_t>(column_begin + k);
         }
         model.hessian_row_indices.insert(model.hessian_row_indices.end(), sorted_rows.begin(), sorted_rows.end());
         for_each_row(j, [&](int row, std::uint32_t& slot) {
            slot = final_slot_of_local[local_slot_of_row[static_cast<std::size_t>(row)]];
         });
         for (int row: column_rows) local_slot_of_row[static_cast<std::size_t>(row)] = no_hessian_slot;
         model.hessian_column_starts[j + 1] = model.hessian_row_indices.size();
      }
   }

   /// Reverse programs of the functions (edges by decreasing source) and their variable leaves. Targets that are
   /// constants are dropped. The programs of consecutive functions are adjacent, so a run of constraints has one.
   void ModelBuilder::build_reverse_programs() {
      ModelData& model = *this->data;
      TapeStorage& tapes = model.tapes;
      const std::size_t number_instructions = tapes.instructions.size();
      model.number_partials = 2 * number_instructions + tapes.arguments.size() + 1;
      if (model.number_partials >= ReverseEdge::assign_flag) {
         throw std::runtime_error("cppasl: the tapes are too large (more than 2^31 partial derivatives)");
      }
      // upper bounds (reserving only maps virtual memory: untouched pages cost nothing)
      model.reverse_edges.reserve(number_instructions + tapes.arguments.size());
      model.variable_occurrences.reserve(number_instructions);
      // reached[i]: instruction i receives an adjoint (tape outputs, then targets of kept edges)
      std::vector<char> reached(number_instructions, 0);
      for (const NonlinearElement& element: model.elements) reached[element.instruction_begin + element.instruction_count - 1] = 1;
      const auto one = static_cast<std::uint32_t>(model.number_partials - 1);
      const auto operand_partials_begin = static_cast<std::uint32_t>(2 * number_instructions);
      auto is_constant = [&](std::int32_t instruction) {
         return tapes.instructions[static_cast<std::size_t>(instruction)].operation == TapeOperation::Constant;
      };
      // edges of the instructions [begin, end), by decreasing source; leaves are reported to `leaf`
      auto build_range = [&](std::size_t begin, std::size_t end, auto&& leaf) {
         for (std::size_t i = end; i-- > begin;) {
            const TapeInstruction& instruction = tapes.instructions[i];
            const auto source = static_cast<std::uint32_t>(i);
            if (!reached[i]) { // its adjoint is always 0 (e.g. inside a piecewise-constant operation)
               if (instruction.operation == TapeOperation::DefinedValue) { // an empty link for the Hessian kernels
                  tapes.instructions[i].constant = static_cast<double>(model.defined_leaf_links.size());
                  model.defined_leaf_links.push_back({0, 0});
               }
               continue;
            }
            auto add = [&](std::int32_t target, std::uint32_t partial) {
               if (is_constant(target)) return;
               char& target_reached = reached[static_cast<std::size_t>(target)];
               model.reverse_edges.push_back({source, static_cast<std::uint32_t>(target),
                  partial | (target_reached ? 0u : ReverseEdge::assign_flag)});
               target_reached = 1;
            };
            switch (derivative_category(instruction.operation)) {
               case DerivativeCategory::Leaf: leaf(source, instruction); break;
               case DerivativeCategory::Unary: add(instruction.first, 2 * source); break;
               case DerivativeCategory::Binary:
                  add(instruction.first, 2 * source);
                  add(instruction.second, 2 * source + 1);
                  break;
               case DerivativeCategory::Sum:
                  for (std::int32_t k = 0; k < instruction.second; ++k) add(tapes.arguments[static_cast<std::size_t>(instruction.first + k)], one);
                  break;
               case DerivativeCategory::Select:
                  for (std::int32_t k = 0; k < instruction.second; ++k) {
                     add(tapes.arguments[static_cast<std::size_t>(instruction.first + k)],
                        operand_partials_begin + static_cast<std::uint32_t>(instruction.first + k));
                  }
                  break;
               default: break; // piecewise-constant operations: zero derivatives
            }
         }
      };
      std::vector<std::size_t> defined_tape_of_variable(model.shared_defined_outputs.size(), 0);
      for (const auto& [defined, tape_index]: this->defined_tape_index) defined_tape_of_variable[defined] = tape_index;
      std::vector<std::int32_t> local_of_variable(static_cast<std::size_t>(model.number_variables), -1);

      for (FunctionBlock& block: model.functions) {
         block.edge_begin = model.reverse_edges.size();
         block.occurrence_begin = model.variable_occurrences.size();
         block.defined_occurrence_begin = model.defined_occurrences.size();
         for (std::size_t e = block.element_end; e-- > block.element_begin;) {
            const NonlinearElement& element = model.elements[e];
            const int* element_variables = tapes.element_variables.data() + element.variable_begin;
            for (std::size_t v = 0; v < element.variable_count; ++v) local_of_variable[static_cast<std::size_t>(element_variables[v])] = static_cast<std::int32_t>(v);
            build_range(element.instruction_begin, element.instruction_begin + element.instruction_count,
               [&](std::uint32_t source, const TapeInstruction& instruction) {
                  if (instruction.operation == TapeOperation::Variable) {
                     model.variable_occurrences.push_back({source, model.element_variable_positions[element.variable_begin +
                        static_cast<std::size_t>(instruction.second)]});
                  }
                  else if (instruction.operation == TapeOperation::DefinedValue) {
                     // the gradient of the defined variable lands on the positions of its variables in this function
                     const std::size_t tape_index = defined_tape_of_variable[static_cast<std::size_t>(instruction.second)];
                     const DefinedTape& defined = model.defined_tapes[tape_index];
                     tapes.instructions[source].constant = static_cast<double>(model.defined_leaf_links.size());
                     model.defined_leaf_links.push_back({model.defined_occurrence_positions.size(), static_cast<std::uint32_t>(defined.variable_count)});
                     model.defined_occurrences.push_back({source, static_cast<std::uint32_t>(defined.variable_count),
                        model.defined_occurrence_positions.size(), defined.gradient_begin});
                     for (std::size_t v = 0; v < defined.variable_count; ++v) {
                        const int variable = tapes.element_variables[defined.variable_begin + v];
                        const std::int32_t local = local_of_variable[static_cast<std::size_t>(variable)];
                        model.defined_occurrence_positions.push_back(model.element_variable_positions[element.variable_begin + static_cast<std::size_t>(local)]);
                        model.defined_pair_locals.push_back(static_cast<std::uint32_t>(local));
                        model.defined_pair_gradients.push_back(static_cast<std::uint32_t>(defined.gradient_begin + v));
                     }
                  }
               });
            for (std::size_t v = 0; v < element.variable_count; ++v) local_of_variable[static_cast<std::size_t>(element_variables[v])] = -1;
         }
         block.edge_end = model.reverse_edges.size();
         block.occurrence_end = model.variable_occurrences.size();
         block.defined_occurrence_end = model.defined_occurrences.size();
      }
      // program of the shared defined variables (their own storage): gradients w.r.t. their own variables
      const TapeStorage& defined_storage = model.defined_storage;
      const std::size_t number_defined_instructions = defined_storage.instructions.size();
      model.number_defined_partials = 2 * number_defined_instructions + defined_storage.arguments.size() + 1;
      std::vector<char> defined_reached(number_defined_instructions, 0);
      for (const DefinedTape& defined: model.defined_tapes) defined_reached[defined.instruction_begin + defined.instruction_count - 1] = 1;
      const auto defined_one = static_cast<std::uint32_t>(model.number_defined_partials - 1);
      for (std::size_t t = model.defined_tapes.size(); t-- > 0;) {
         const DefinedTape& defined = model.defined_tapes[t];
         for (std::size_t i = defined.instruction_begin + defined.instruction_count; i-- > defined.instruction_begin;) {
            const TapeInstruction& instruction = defined_storage.instructions[i];
            if (!defined_reached[i]) continue;
            const auto source = static_cast<std::uint32_t>(i);
            auto add = [&](std::int32_t target, std::uint32_t partial) {
               if (defined_storage.instructions[static_cast<std::size_t>(target)].operation == TapeOperation::Constant) return;
               char& target_reached = defined_reached[static_cast<std::size_t>(target)];
               model.defined_reverse_edges.push_back({source, static_cast<std::uint32_t>(target),
                  partial | (target_reached ? 0u : ReverseEdge::assign_flag)});
               target_reached = 1;
            };
            switch (derivative_category(instruction.operation)) {
               case DerivativeCategory::Leaf:
                  if (instruction.operation == TapeOperation::Variable) {
                     model.defined_variable_occurrences.push_back({source, static_cast<std::uint32_t>(defined.gradient_begin +
                        static_cast<std::size_t>(instruction.second))});
                  }
                  break;
               case DerivativeCategory::Unary: add(instruction.first, 2 * source); break;
               case DerivativeCategory::Binary:
                  add(instruction.first, 2 * source);
                  add(instruction.second, 2 * source + 1);
                  break;
               case DerivativeCategory::Sum:
                  for (std::int32_t k = 0; k < instruction.second; ++k) add(defined_storage.arguments[static_cast<std::size_t>(instruction.first + k)], defined_one);
                  break;
               case DerivativeCategory::Select:
                  for (std::int32_t k = 0; k < instruction.second; ++k) {
                     add(defined_storage.arguments[static_cast<std::size_t>(instruction.first + k)],
                        static_cast<std::uint32_t>(2 * number_defined_instructions) + static_cast<std::uint32_t>(instruction.first + k));
                  }
                  break;
               default: break;
            }
         }
      }
   }

   /// AMPL orders the variables: nonlinear in constraints and objectives, nonlinear in constraints only, nonlinear in
   /// objectives only (the integer ones last in each group), linear arcs, other linear, binary, integer.
   void ModelBuilder::compute_variable_integrality() {
      ModelData& model = *this->data;
      const NlHeader& h = model.header;
      model.variable_is_integer.assign(static_cast<std::size_t>(model.number_variables), 0);
      auto mark_integer = [&](int begin, int end) {
         for (int i = std::max(begin, 0); i < std::min(end, model.number_variables); ++i) {
            model.variable_is_integer[static_cast<std::size_t>(i)] = 1;
         }
      };
      const int in_both = h.number_nonlinear_variables_in_both;
      const int constraints_only = h.number_nonlinear_variables_in_constraints - in_both;
      const int objectives_only = h.number_nonlinear_variables_in_objectives - in_both;
      int group_end = in_both;
      mark_integer(group_end - h.number_nonlinear_integer_variables_in_both, group_end);
      group_end += std::max(constraints_only, 0);
      mark_integer(group_end - h.number_nonlinear_integer_variables_in_constraints, group_end);
      group_end += std::max(objectives_only, 0);
      mark_integer(group_end - h.number_nonlinear_integer_variables_in_objectives, group_end);
      mark_integer(model.number_variables - h.number_binary_variables - h.number_integer_variables, model.number_variables);
   }

} // namespace cppasl::detail
