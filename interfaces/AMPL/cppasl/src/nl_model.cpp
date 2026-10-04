#include <algorithm>
#include <charconv>
#include <fstream>
#include <string>
#include <cstring>
#include <limits>
#include <stdexcept>
#include "cppasl/nl_model.hpp"
#include "model_data.hpp"

namespace cppasl {

   using detail::FunctionBlock;
   using detail::InstructionPartials;
   using detail::ModelData;
   using detail::NonlinearElement;

   namespace {
      /// 1/2 x^T H x for the quadratic part of a function
      inline double quadratic_value(const ModelData& model, const FunctionBlock& block, const double* x) {
         const int* rows = model.quadratic_rows.data();
         const int* columns = model.quadratic_columns.data();
         const double* values = model.quadratic_values.data();
         double diagonal_sum = 0., off_diagonal_sum = 0.;
         for (std::size_t q = block.quadratic_begin; q < block.quadratic_off_diagonal_begin; ++q) {
            const double xi = x[rows[q]];
            diagonal_sum += values[q] * xi * xi;
         }
         for (std::size_t q = block.quadratic_off_diagonal_begin; q < block.quadratic_end; ++q) {
            off_diagonal_sum += values[q] * x[rows[q]] * x[columns[q]];
         }
         return 0.5 * diagonal_sum + off_diagonal_sum;
      }

      /// gradient[positions] += weight * H x for the quadratic part of a function (positions: variables for
      /// objectives, Jacobian nonzeros for constraints)
      inline void add_quadratic_gradient(const ModelData& model, const FunctionBlock& block, const double* x, double weight,
            double* gradient) {
         const int* rows = model.quadratic_rows.data();
         const int* columns = model.quadratic_columns.data();
         const double* values = model.quadratic_values.data();
         const std::uint32_t* row_positions = model.quadratic_row_positions.data();
         const std::uint32_t* column_positions = model.quadratic_column_positions.data();
         for (std::size_t q = block.quadratic_begin; q < block.quadratic_off_diagonal_begin; ++q) {
            gradient[row_positions[q]] += weight * values[q] * x[rows[q]];
         }
         for (std::size_t q = block.quadratic_off_diagonal_begin; q < block.quadratic_end; ++q) {
            const double h = weight * values[q];
            gradient[row_positions[q]] += h * x[columns[q]];
            gradient[column_positions[q]] += h * x[rows[q]];
         }
      }

      /// result += weight * H v (full symmetric product with the lower triangle, indexed by variables)
      inline void add_quadratic_product(const ModelData& model, const FunctionBlock& block, double weight, const double* v,
            double* result) {
         const int* rows = model.quadratic_rows.data();
         const int* columns = model.quadratic_columns.data();
         const double* values = model.quadratic_values.data();
         for (std::size_t q = block.quadratic_begin; q < block.quadratic_off_diagonal_begin; ++q) {
            result[rows[q]] += weight * values[q] * v[rows[q]];
         }
         for (std::size_t q = block.quadratic_off_diagonal_begin; q < block.quadratic_end; ++q) {
            const double h = weight * values[q];
            result[rows[q]] += h * v[columns[q]];
            result[columns[q]] += h * v[rows[q]];
         }
      }
   } // namespace

   /// Evaluation kernels with access to the workspace internals
   class EvaluationKernels {
   public:
      /// records x in the workspace; the tape values of all functions become stale if x changed
      static void set_point(EvaluationWorkspace& workspace, const double* x) {
         const std::size_t bytes = workspace.current_point.size() * sizeof(double);
         if (std::memcmp(workspace.current_point.data(), x, bytes) != 0) {
            std::memcpy(workspace.current_point.data(), x, bytes);
            ++workspace.point_version;
         }
      }

      /// forward sweep of all the shared defined variables (once per point, before any function reads them)
      static void evaluate_defined_variables(const ModelData& model, EvaluationWorkspace& workspace, const double* x) {
         if (model.defined_tapes.empty() || workspace.defined_point_version == workspace.point_version) return;
         detail::evaluate_tape_values_and_partials(model.defined_tapes_view(), x, workspace.defined_values.data(),
            workspace.defined_partials.data(), workspace.defined_partials.data() + 2 * model.defined_storage.instructions.size());
         workspace.defined_point_version = workspace.point_version;
      }

      /// gradients of the shared defined variables w.r.t. their own variables (once per point): one reverse program
      static void evaluate_defined_gradients(const ModelData& model, EvaluationWorkspace& workspace) {
         if (workspace.defined_gradient_version == workspace.point_version) return;
         double* adjoints = workspace.defined_adjoints.data();
         double* gradients = workspace.defined_gradients.data();
         for (const detail::DefinedTape& defined: model.defined_tapes) adjoints[defined.instruction_begin + defined.instruction_count - 1] = 1.;
         detail::propagate_adjoints(model.defined_reverse_edges.data(), model.defined_reverse_edges.data() + model.defined_reverse_edges.size(),
            adjoints, workspace.defined_partials.data());
         std::fill(gradients, gradients + model.number_defined_gradients, 0.);
         for (const detail::VariableOccurrence& occurrence: model.defined_variable_occurrences) {
            gradients[occurrence.position] += adjoints[occurrence.instruction];
         }
         workspace.defined_gradient_version = workspace.point_version;
      }

      /// forward sweep of all the element tapes of a function, with first partials (cached per point)
      static void evaluate_elements(const ModelData& model, EvaluationWorkspace& workspace, std::size_t function_index,
            const double* x) {
         if (workspace.function_point_version[function_index] == workspace.point_version) return;
         EvaluationKernels::evaluate_defined_variables(model, workspace, x);
         const FunctionBlock& block = model.functions[function_index];
         detail::evaluate_tape_values_and_partials(model.tape_of(block), x, workspace.tape_values.data(),
            workspace.tape_partials.data(), EvaluationKernels::operand_partials(model, workspace), workspace.defined_values.data());
         workspace.function_point_version[function_index] = workspace.point_version;
      }

      /// forward sweep of a run of constraints (cached per point)
      static void evaluate_run(const ModelData& model, EvaluationWorkspace& workspace, const detail::ConstraintRun& run,
            const double* x) {
         const auto offset = static_cast<std::size_t>(model.number_objectives);
         bool is_current = true;
         for (std::size_t k = run.index_begin; k < run.index_end && is_current; ++k) {
            is_current = workspace.function_point_version[offset + static_cast<std::size_t>(model.nonlinear_constraint_indices[k])] == workspace.point_version;
         }
         if (is_current) return;
         EvaluationKernels::evaluate_defined_variables(model, workspace, x);
         detail::evaluate_tape_values_and_partials(model.tape_of(run), x, workspace.tape_values.data(),
            workspace.tape_partials.data(), EvaluationKernels::operand_partials(model, workspace), workspace.defined_values.data());
         for (std::size_t k = run.index_begin; k < run.index_end; ++k) {
            workspace.function_point_version[offset + static_cast<std::size_t>(model.nonlinear_constraint_indices[k])] = workspace.point_version;
         }
      }

      static double* operand_partials(const ModelData& model, EvaluationWorkspace& workspace) {
         return workspace.tape_partials.data() + 2 * model.total_tape_length;
      }

      /// The gradients of the elements of the functions [first, last] (whose tapes and reverse programs are adjacent),
      /// weighted by weight * scale, added to their positions: one branch-free reverse program.
      static void add_gradients(const ModelData& model, EvaluationWorkspace& workspace, const FunctionBlock& first,
            const FunctionBlock& last, double weight, double* gradient) {
         double* adjoints = workspace.adjoints.data(); // no zeroing: the first edge into an instruction assigns
         const std::uint32_t* outputs = model.element_outputs.data();
         const double* scales = model.element_scales.data();
         for (std::size_t e = first.element_begin; e < last.element_end; ++e) adjoints[outputs[e]] = weight * scales[e];
         detail::propagate_adjoints(model.reverse_edges.data() + first.edge_begin, model.reverse_edges.data() + last.edge_end,
            adjoints, workspace.tape_partials.data());
         const detail::VariableOccurrence* occurrence = model.variable_occurrences.data() + first.occurrence_begin;
         const detail::VariableOccurrence* end = model.variable_occurrences.data() + last.occurrence_end;
         for (; occurrence != end; ++occurrence) gradient[occurrence->position] += adjoints[occurrence->instruction];
         // shared defined variables: adjoint * gradient of the defined variable
         if (first.defined_occurrence_begin < last.defined_occurrence_end) {
            EvaluationKernels::evaluate_defined_gradients(model, workspace);
            for (std::size_t k = first.defined_occurrence_begin; k < last.defined_occurrence_end; ++k) {
               const detail::DefinedOccurrence& occurrence = model.defined_occurrences[k];
               const double adjoint = adjoints[occurrence.instruction];
               if (adjoint == 0.) continue;
               const double* defined_gradient = workspace.defined_gradients.data() + occurrence.gradient_begin;
               const std::uint32_t* positions = model.defined_occurrence_positions.data() + occurrence.position_begin;
               for (std::size_t v = 0; v < occurrence.variable_count; ++v) gradient[positions[v]] += adjoint * defined_gradient[v];
            }
         }
      }


      /// Jacobian entries of a run of constraints
      static void add_run_gradients(const ModelData& model, EvaluationWorkspace& workspace, const detail::ConstraintRun& run,
            double* jacobian_values) {
         add_gradients(model, workspace, model.constraint_block(model.nonlinear_constraint_indices[run.index_begin]),
            model.constraint_block(model.nonlinear_constraint_indices[run.index_end - 1]), 1., jacobian_values);
      }

      static double element_sum(const ModelData& model, const EvaluationWorkspace& workspace, std::size_t element_begin,
            std::size_t element_end) {
         const std::uint32_t* outputs = model.element_outputs.data();
         const double* scales = model.element_scales.data();
         const double* values = workspace.tape_values.data();
         double sum = 0.;
         for (std::size_t e = element_begin; e < element_end; ++e) sum += scales[e] * values[outputs[e]];
         return sum;
      }

      static InstructionPartials* partials(EvaluationWorkspace& workspace) {
         return reinterpret_cast<InstructionPartials*>(workspace.partials.data());
      }

      /// gradient[positions] += weight * sum_e scale_e grad(element_e) for all elements of a function (values evaluated)
      static void add_element_gradients(const ModelData& model, EvaluationWorkspace& workspace, const FunctionBlock& block,
            double weight, double* gradient) {
         add_gradients(model, workspace, block, block, weight, gradient);
      }

      /// hessian_values[slots] += weight * scale * Hessian(element): group formula phi'' g g^T + phi' sum H_k for group
      /// elements, vector forward-over-reverse sweeps (up to 8 directions at once) otherwise
      static void add_element_hessians(const ModelData& model, EvaluationWorkspace& workspace, const FunctionBlock& block,
            double weight, double* hessian_values) {
         InstructionPartials* partials = EvaluationKernels::partials(workspace);
         double* packed = workspace.element_hessian.data();
         const double* values = workspace.tape_values.data();
         detail::DefinedHessianLinks links{model.defined_leaf_links.data(), model.defined_pair_locals.data(),
            model.defined_pair_gradients.data(), workspace.defined_gradients.data(), workspace.defined_weights.data()};
         for (std::size_t e = block.element_begin; e < block.element_end; ++e) {
            const NonlinearElement& element = model.elements[e];
            const detail::TapeView tape = model.tape_of(element);
            if (element.has_defined_leaves) EvaluationKernels::evaluate_defined_gradients(model, workspace);
            const std::size_t size = element.variable_count * (element.variable_count + 1) / 2;
            const double factor = weight * element.scale;
            std::fill(packed, packed + size, 0.);
            bool done = false;
            if (element.group_sum != NonlinearElement::no_group) {
               done = detail::evaluate_group_hessian(tape, element.group_sum,
                  model.tapes.group_operand_ranges.data() + element.group_range_begin, element.group_range_end - element.group_range_begin,
                  values, partials, workspace.adjoints.data(), workspace.tangents.data(), workspace.second_order_adjoints.data(),
                  workspace.direction_of_local.data(), workspace.local_of_direction.data(), workspace.local_direction.data(),
                  element.variable_count, factor, packed);
            }
            if (!done) {
               detail::evaluate_tape_hessian(tape, values, partials, workspace.adjoints.data(), workspace.tangents.data(),
                  workspace.second_order_adjoints.data(), element.variable_count, factor, packed,
                  element.has_defined_leaves ? &links : nullptr);
            }
            const std::uint32_t* slots = model.element_hessian_slots.data() + element.hessian_slot_begin;
            for (std::size_t k = 0; k < size; ++k) hessian_values[slots[k]] += packed[k];
         }
      }

      /// Prepares the Hessian-vector products through the shared defined variables: their tangents along the vector
      /// (from their gradients) and zero seeds. Returns the links passed to the element kernels.
      static detail::DefinedVariableLinks prepare_defined_links(const ModelData& model, EvaluationWorkspace& workspace,
            const double* x, const double* vector) {
         EvaluationKernels::evaluate_defined_variables(model, workspace, x);
         EvaluationKernels::evaluate_defined_gradients(model, workspace);
         const std::size_t length = model.defined_storage.instructions.size();
         std::fill(workspace.defined_adjoints.begin(), workspace.defined_adjoints.end(), 0.);
         std::fill(workspace.defined_second_order_adjoints.begin(), workspace.defined_second_order_adjoints.begin() + static_cast<std::ptrdiff_t>(length), 0.);
         for (const detail::DefinedTape& defined: model.defined_tapes) {
            const int* variables = model.tapes.element_variables.data() + defined.variable_begin;
            const double* gradient = workspace.defined_gradients.data() + defined.gradient_begin;
            double tangent = 0.;
            for (std::size_t v = 0; v < defined.variable_count; ++v) tangent += gradient[v] * vector[variables[v]];
            workspace.defined_direction_tangents[defined.instruction_begin + defined.instruction_count - 1] = tangent;
         }
         return {workspace.defined_direction_tangents.data(), workspace.defined_adjoints.data(), workspace.defined_second_order_adjoints.data()};
      }

      /// hessian_values[slots] += sum over the shared defined variables of weight_v * Hessian(v) (after the element kernels)
      static void add_defined_hessians(const ModelData& model, EvaluationWorkspace& workspace, double* hessian_values) {
         InstructionPartials* partials = EvaluationKernels::partials(workspace);
         double* packed = workspace.element_hessian.data();
         for (const detail::DefinedTape& defined: model.defined_tapes) {
            const double weight = workspace.defined_weights[defined.instruction_begin + defined.instruction_count - 1];
            if (weight == 0. || !defined.contributes_to_hessian) continue;
            const std::size_t size = defined.variable_count * (defined.variable_count + 1) / 2;
            std::fill(packed, packed + size, 0.);
            detail::evaluate_tape_hessian(model.tape_of(defined), workspace.defined_values.data(), partials, workspace.adjoints.data(),
               workspace.tangents.data(), workspace.second_order_adjoints.data(), defined.variable_count, weight, packed);
            const std::uint32_t* slots = model.defined_hessian_slots.data() + defined.hessian_slot_begin;
            for (std::size_t k = 0; k < size; ++k) hessian_values[slots[k]] += packed[k];
         }
      }

      /// result += Hessian of the seeded shared defined variables * vector (after the element kernels)
      static void add_defined_hessian_vector_product(const ModelData& model, EvaluationWorkspace& workspace,
            const double* vector, double* result) {
         detail::evaluate_seeded_hessian_vector_product(model.defined_tapes_view(), workspace.defined_values.data(),
            reinterpret_cast<InstructionPartials*>(workspace.defined_second_partials.data()), workspace.defined_adjoints.data(),
            workspace.defined_tangents.data(), workspace.defined_second_order_adjoints.data(), vector, result);
      }

      /// result[variables] += weight * scale * Hessian(element) * vector
      static void add_element_hessian_vector_products(const ModelData& model, EvaluationWorkspace& workspace,
            const FunctionBlock& block, double weight, const double* vector, double* result,
            const detail::DefinedVariableLinks* links) {
         InstructionPartials* partials = EvaluationKernels::partials(workspace);
         for (std::size_t e = block.element_begin; e < block.element_end; ++e) {
            const NonlinearElement& element = model.elements[e];
            const int* variables = model.tapes.element_variables.data() + element.variable_begin;
            bool zero_direction = true;
            for (std::size_t v = 0; v < element.variable_count; ++v) {
               workspace.local_direction[v] = vector[variables[v]];
               zero_direction = zero_direction && (workspace.local_direction[v] == 0.);
            }
            if (element.has_defined_leaves) {
               // through the shared defined variables: the DefinedValue leaves are linked
               detail::evaluate_tape_hessian_vector_product(model.tape_of(element), workspace.tape_values.data(),
                  partials, workspace.adjoints.data(), workspace.tangents.data(), workspace.second_order_adjoints.data(),
                  workspace.local_direction.data(), weight * element.scale, result, links);
               continue;
            }
            if (zero_direction) continue;
            detail::evaluate_tape_hessian_vector_product(model.tape_of(element), workspace.tape_values.data(), partials, workspace.adjoints.data(),
               workspace.tangents.data(), workspace.second_order_adjoints.data(), workspace.local_direction.data(),
               weight * element.scale, result);
         }
      }
   };

   // ------------------------------------------------------------------------------------------------------------------
   // construction and accessors
   // ------------------------------------------------------------------------------------------------------------------

   NlModel::NlModel(std::shared_ptr<const detail::ModelData> data): data(std::move(data)) {}

   EvaluationWorkspace::EvaluationWorkspace(const NlModel& model): data(model.data) {
      const ModelData& d = *this->data;
      this->tape_values.assign(d.total_tape_length, 0.);
      this->defined_gradients.assign(d.number_defined_gradients, 0.);
      this->defined_values.assign(d.defined_storage.instructions.size(), 0.);
      this->defined_adjoints.assign(d.defined_storage.instructions.size(), 0.);
      this->defined_partials.assign(d.number_defined_partials, 0.);
      const std::size_t number_defined_instructions = d.defined_storage.instructions.size();
      this->defined_tangents.assign(number_defined_instructions, 0.);
      this->defined_direction_tangents.assign(number_defined_instructions, 0.);
      this->defined_second_order_adjoints.assign(number_defined_instructions, 0.);
      this->defined_second_partials.assign(5 * number_defined_instructions, 0.);
      this->defined_weights.assign(number_defined_instructions, 0.);
      if (d.number_defined_partials > 0) this->defined_partials.back() = 1.;
      this->tape_partials.assign(d.number_partials, 0.);
      if (d.number_partials > 0) this->tape_partials.back() = 1.; // the partials of sums
      this->function_point_version.assign(d.functions.size(), 0);
      this->current_point.assign(static_cast<std::size_t>(d.number_variables), std::numeric_limits<double>::quiet_NaN());
      this->partials.assign(5 * d.maximum_element_tape_length, 0.);
      this->adjoints.assign(d.total_tape_length, 0.); // one-sweep gradients of whole runs of constraints
      this->tangents.assign(detail::maximum_hessian_directions * d.maximum_element_tape_length, 0.);
      this->second_order_adjoints.assign(detail::maximum_hessian_directions * d.maximum_element_tape_length, 0.);
      this->local_direction.assign(d.maximum_element_variables, 0.);
      this->element_hessian.assign(d.maximum_element_variables * (d.maximum_element_variables + 1) / 2, 0.);
      this->direction_of_local.assign(d.maximum_element_variables, -1);
      this->local_of_direction.assign(d.maximum_element_variables, 0);
   }

   const NlHeader& NlModel::header() const { return this->data->header; }
   int NlModel::number_variables() const { return this->data->number_variables; }
   int NlModel::number_constraints() const { return this->data->number_constraints; }
   int NlModel::number_objectives() const { return this->data->number_objectives; }
   const std::vector<double>& NlModel::variable_lower_bounds() const { return this->data->variable_lower_bounds; }
   const std::vector<double>& NlModel::variable_upper_bounds() const { return this->data->variable_upper_bounds; }
   const std::vector<double>& NlModel::constraint_lower_bounds() const { return this->data->constraint_lower_bounds; }
   const std::vector<double>& NlModel::constraint_upper_bounds() const { return this->data->constraint_upper_bounds; }
   const std::vector<double>& NlModel::initial_primal_point() const { return this->data->initial_primal_point; }
   const std::vector<double>& NlModel::initial_dual_point() const { return this->data->initial_dual_point; }
   const std::vector<int>& NlModel::complementary_variables() const { return this->data->complementary_variables; }
   bool NlModel::is_integer_variable(int variable_index) const {
      return this->data->variable_is_integer[static_cast<std::size_t>(variable_index)] != 0;
   }
   bool NlModel::is_maximization(int objective_index) const {
      return this->data->objective_is_maximization[static_cast<std::size_t>(objective_index)] != 0;
   }
   const std::vector<Suffix>& NlModel::suffixes() const { return this->data->suffixes; }

   namespace {
      FunctionStructure structure_of(const FunctionBlock& block) {
         if (block.has_elements()) return FunctionStructure::Nonlinear;
         if (block.has_quadratic_part()) return FunctionStructure::Quadratic;
         return FunctionStructure::Linear;
      }
      QuadraticFormView quadratic_form_of(const ModelData& model, const FunctionBlock& block) {
         QuadraticFormView view;
         view.number_nonzeros = block.quadratic_end - block.quadratic_begin;
         view.number_diagonal_nonzeros = block.quadratic_off_diagonal_begin - block.quadratic_begin;
         view.row_indices = model.quadratic_rows.data() + block.quadratic_begin;
         view.column_indices = model.quadratic_columns.data() + block.quadratic_begin;
         view.values = model.quadratic_values.data() + block.quadratic_begin;
         return view;
      }
   } // namespace

   FunctionStructure NlModel::objective_structure(int objective_index) const {
      return structure_of(this->data->objective_block(objective_index));
   }
   FunctionStructure NlModel::constraint_structure(int constraint_index) const {
      return structure_of(this->data->constraint_block(constraint_index));
   }
   QuadraticFormView NlModel::objective_quadratic_form(int objective_index) const {
      return quadratic_form_of(*this->data, this->data->objective_block(objective_index));
   }
   QuadraticFormView NlModel::constraint_quadratic_form(int constraint_index) const {
      return quadratic_form_of(*this->data, this->data->constraint_block(constraint_index));
   }

   std::size_t NlModel::number_jacobian_nonzeros() const { return this->data->jacobian_column_indices.size(); }
   const std::vector<std::size_t>& NlModel::jacobian_row_starts() const { return this->data->jacobian_row_starts; }
   const std::vector<int>& NlModel::jacobian_row_indices() const {
      const ModelData& model = *this->data;
      std::call_once(model.jacobian_row_indices_built, [&] {
         model.jacobian_row_indices.resize(model.jacobian_column_indices.size());
         for (int i = 0; i < model.number_constraints; ++i) {
            std::fill(model.jacobian_row_indices.begin() + static_cast<std::ptrdiff_t>(model.jacobian_row_starts[static_cast<std::size_t>(i)]),
               model.jacobian_row_indices.begin() + static_cast<std::ptrdiff_t>(model.jacobian_row_starts[static_cast<std::size_t>(i) + 1]), i);
         }
      });
      return model.jacobian_row_indices;
   }
   const std::vector<int>& NlModel::jacobian_column_indices() const { return this->data->jacobian_column_indices; }
   int NlModel::hessian_objective_index() const { return this->data->hessian_objective_index; }
   std::size_t NlModel::number_hessian_nonzeros() const { return this->data->hessian_row_indices.size(); }
   const std::vector<std::size_t>& NlModel::hessian_column_starts() const { return this->data->hessian_column_starts; }
   const std::vector<int>& NlModel::hessian_row_indices() const { return this->data->hessian_row_indices; }
   const std::vector<int>& NlModel::hessian_column_indices() const {
      const ModelData& model = *this->data;
      std::call_once(model.hessian_column_indices_built, [&] {
         model.hessian_column_indices.resize(model.hessian_row_indices.size());
         for (int j = 0; j < model.number_variables; ++j) {
            std::fill(model.hessian_column_indices.begin() + static_cast<std::ptrdiff_t>(model.hessian_column_starts[static_cast<std::size_t>(j)]),
               model.hessian_column_indices.begin() + static_cast<std::ptrdiff_t>(model.hessian_column_starts[static_cast<std::size_t>(j) + 1]), j);
         }
      });
      return model.hessian_column_indices;
   }

   // ------------------------------------------------------------------------------------------------------------------
   // evaluations
   // ------------------------------------------------------------------------------------------------------------------

   double NlModel::evaluate_objective(EvaluationWorkspace& workspace, const double* x, int objective_index) const {
      const ModelData& model = *this->data;
      const FunctionBlock& block = model.objective_block(objective_index);
      double value = block.constant;
      for (std::size_t k = block.linear_begin; k < block.linear_end; ++k) {
         value += model.objective_linear_coefficients[k] * x[model.objective_linear_variables[k]];
      }
      if (block.has_quadratic_part()) value += quadratic_value(model, block, x);
      if (block.has_elements()) {
         EvaluationKernels::set_point(workspace, x);
         EvaluationKernels::evaluate_elements(model, workspace, static_cast<std::size_t>(objective_index), x);
         value += EvaluationKernels::element_sum(model, workspace, block.element_begin, block.element_end);
      }
      return value;
   }

   void NlModel::evaluate_objective_gradient(EvaluationWorkspace& workspace, const double* x, double* gradient,
         int objective_index) const {
      const ModelData& model = *this->data;
      const FunctionBlock& block = model.objective_block(objective_index);
      std::fill(gradient, gradient + model.number_variables, 0.);
      for (std::size_t k = block.linear_begin; k < block.linear_end; ++k) {
         gradient[model.objective_linear_variables[k]] += model.objective_linear_coefficients[k];
      }
      if (block.has_quadratic_part()) add_quadratic_gradient(model, block, x, 1., gradient);
      if (block.has_elements()) {
         EvaluationKernels::set_point(workspace, x);
         EvaluationKernels::evaluate_elements(model, workspace, static_cast<std::size_t>(objective_index), x);
         EvaluationKernels::add_element_gradients(model, workspace, block, 1., gradient);
      }
   }

   void NlModel::evaluate_constraints(EvaluationWorkspace& workspace, const double* x, double* constraints) const {
      const ModelData& model = *this->data;
      const std::size_t* row_starts = model.jacobian_row_starts.data();
      const int* columns = model.jacobian_column_indices.data();
      const double* coefficients = model.jacobian_linear_coefficients.data();
      // constant + linear part: sparse matrix-vector product with the aligned coefficients
      for (int i = 0; i < model.number_constraints; ++i) {
         double value = model.constraint_constants[static_cast<std::size_t>(i)];
         for (std::size_t k = row_starts[i]; k < row_starts[i + 1]; ++k) value += coefficients[k] * x[columns[k]];
         constraints[i] = value;
      }
      for (int i: model.quadratic_constraint_indices) constraints[i] += quadratic_value(model, model.constraint_block(i), x);
      if (!model.nonlinear_constraint_indices.empty()) {
         EvaluationKernels::set_point(workspace, x);
         for (const detail::ConstraintRun& run: model.constraint_runs) {
            EvaluationKernels::evaluate_run(model, workspace, run, x);
            for (std::size_t k = run.index_begin; k < run.index_end; ++k) {
               const int i = model.nonlinear_constraint_indices[k];
               const auto ui = static_cast<std::size_t>(i);
               constraints[i] += EvaluationKernels::element_sum(model, workspace, model.constraint_element_begin[ui],
                  model.constraint_element_end[ui]);
            }
         }
      }
   }

   void NlModel::evaluate_jacobian(EvaluationWorkspace& workspace, const double* x, double* jacobian_values) const {
      const ModelData& model = *this->data;
      if (!model.jacobian_linear_coefficients.empty()) {
         std::memcpy(jacobian_values, model.jacobian_linear_coefficients.data(), model.jacobian_linear_coefficients.size() * sizeof(double));
      }
      for (int i: model.quadratic_constraint_indices) add_quadratic_gradient(model, model.constraint_block(i), x, 1., jacobian_values);
      if (!model.nonlinear_constraint_indices.empty()) {
         EvaluationKernels::set_point(workspace, x);
         for (const detail::ConstraintRun& run: model.constraint_runs) {
            EvaluationKernels::evaluate_run(model, workspace, run, x);
            EvaluationKernels::add_run_gradients(model, workspace, run, jacobian_values);
         }
      }
   }

   void NlModel::evaluate_lagrangian_hessian(EvaluationWorkspace& workspace, const double* x, double objective_multiplier,
         const double* constraint_multipliers, double* hessian_values) const {
      const ModelData& model = *this->data;
      std::fill(hessian_values, hessian_values + model.hessian_row_indices.size(), 0.);
      auto add_function = [&](std::size_t function_index, double weight) {
         const FunctionBlock& block = model.functions[function_index];
         const double* values = model.quadratic_values.data();
         const std::uint32_t* slots = model.quadratic_hessian_slots.data();
         for (std::size_t q = block.quadratic_begin; q < block.quadratic_end; ++q) hessian_values[slots[q]] += weight * values[q];
         if (block.has_elements()) {
            EvaluationKernels::evaluate_elements(model, workspace, function_index, x);
            EvaluationKernels::add_element_hessians(model, workspace, block, weight, hessian_values);
         }
      };
      EvaluationKernels::set_point(workspace, x);
      std::fill(workspace.defined_weights.begin(), workspace.defined_weights.end(), 0.);
      if (model.hessian_objective_index >= 0 && objective_multiplier != 0.) {
         add_function(static_cast<std::size_t>(model.hessian_objective_index), objective_multiplier);
      }
      const auto offset = static_cast<std::size_t>(model.number_objectives);
      for (int i: model.quadratic_constraint_indices) {
         const double y = constraint_multipliers[i];
         if (y == 0.) continue;
         const FunctionBlock& block = model.constraint_block(i);
         const double* values = model.quadratic_values.data();
         const std::uint32_t* slots = model.quadratic_hessian_slots.data();
         for (std::size_t q = block.quadratic_begin; q < block.quadratic_end; ++q) hessian_values[slots[q]] += y * values[q];
      }
      for (int i: model.nonlinear_constraint_indices) {
         const double y = constraint_multipliers[i];
         if (y == 0.) continue;
         const std::size_t function_index = offset + static_cast<std::size_t>(i);
         EvaluationKernels::evaluate_elements(model, workspace, function_index, x);
         EvaluationKernels::add_element_hessians(model, workspace, model.functions[function_index], y, hessian_values);
      }
      // curvature of the shared defined variables, weighted by the adjoints collected from all the functions
      if (!model.defined_tapes.empty()) EvaluationKernels::add_defined_hessians(model, workspace, hessian_values);
   }

   void NlModel::evaluate_lagrangian_hessian_vector_product(EvaluationWorkspace& workspace, const double* x,
         double objective_multiplier, const double* constraint_multipliers, const double* vector, double* result) const {
      const ModelData& model = *this->data;
      std::fill(result, result + model.number_variables, 0.);
      EvaluationKernels::set_point(workspace, x);
      // the shared defined variables collect seeds from the elements, then contribute once
      detail::DefinedVariableLinks links{};
      if (!model.defined_tapes.empty()) links = EvaluationKernels::prepare_defined_links(model, workspace, x, vector);
      auto add_function = [&](std::size_t function_index, double weight, bool quadratic, bool elements) {
         const FunctionBlock& block = model.functions[function_index];
         if (quadratic) add_quadratic_product(model, block, weight, vector, result);
         if (elements && block.has_elements()) {
            EvaluationKernels::evaluate_elements(model, workspace, function_index, x);
            EvaluationKernels::add_element_hessian_vector_products(model, workspace, block, weight, vector, result, &links);
         }
      };
      if (model.hessian_objective_index >= 0 && objective_multiplier != 0.) {
         add_function(static_cast<std::size_t>(model.hessian_objective_index), objective_multiplier, true, true);
      }
      const auto offset = static_cast<std::size_t>(model.number_objectives);
      for (int i: model.quadratic_constraint_indices) {
         if (constraint_multipliers[i] != 0.) add_function(offset + static_cast<std::size_t>(i), constraint_multipliers[i], true, false);
      }
      for (int i: model.nonlinear_constraint_indices) {
         if (constraint_multipliers[i] != 0.) add_function(offset + static_cast<std::size_t>(i), constraint_multipliers[i], false, true);
      }
      if (!model.defined_tapes.empty()) EvaluationKernels::add_defined_hessian_vector_product(model, workspace, vector, result);
   }

   // ------------------------------------------------------------------------------------------------------------------
   // linear operators
   // ------------------------------------------------------------------------------------------------------------------

   void NlModel::multiply_jacobian(const double* jacobian_values, const double* vector, double* result) const {
      const ModelData& model = *this->data;
      const std::size_t* row_starts = model.jacobian_row_starts.data();
      const int* columns = model.jacobian_column_indices.data();
      for (int i = 0; i < model.number_constraints; ++i) {
         double sum = 0.;
         for (std::size_t k = row_starts[i]; k < row_starts[i + 1]; ++k) sum += jacobian_values[k] * vector[columns[k]];
         result[i] = sum;
      }
   }

   void NlModel::multiply_jacobian_transpose(const double* jacobian_values, const double* vector, double* result) const {
      const ModelData& model = *this->data;
      const std::size_t* row_starts = model.jacobian_row_starts.data();
      const int* columns = model.jacobian_column_indices.data();
      std::fill(result, result + model.number_variables, 0.);
      for (int i = 0; i < model.number_constraints; ++i) {
         const double vi = vector[i];
         if (vi == 0.) continue;
         for (std::size_t k = row_starts[i]; k < row_starts[i + 1]; ++k) result[columns[k]] += jacobian_values[k] * vi;
      }
   }

   void NlModel::multiply_hessian(const double* hessian_values, const double* vector, double* result) const {
      const ModelData& model = *this->data;
      const std::size_t* column_starts = model.hessian_column_starts.data();
      const int* rows = model.hessian_row_indices.data();
      std::fill(result, result + model.number_variables, 0.);
      for (int j = 0; j < model.number_variables; ++j) {
         const double vj = vector[j];
         double column_sum = 0.;
         for (std::size_t k = column_starts[j]; k < column_starts[j + 1]; ++k) {
            const int i = rows[k];
            result[i] += hessian_values[k] * vj;
            if (i != j) column_sum += hessian_values[k] * vector[i];
         }
         result[j] += column_sum;
      }
   }

   // ------------------------------------------------------------------------------------------------------------------
   // solution file
   // ------------------------------------------------------------------------------------------------------------------

   void NlModel::write_solution(const std::string& file_name, const std::string& message, const double* x,
         const double* y, int solve_result_code, const std::vector<Suffix>& output_suffixes) const {
      const ModelData& model = *this->data;
      const NlHeader& header = model.header;
      std::string contents;
      auto number = [&](double value) {
         char text[32];
         const auto result = std::to_chars(text, text + sizeof(text), value);
         contents.append(text, result.ptr);
         contents += '\n';
      };
      auto integer = [&](long long value) { contents += std::to_string(value); contents += '\n'; };
      // message lines (an empty line would end the message: it is written as a blank)
      std::size_t begin = 0;
      std::string trimmed = message;
      while (!trimmed.empty() && trimmed.back() == '\n') trimmed.pop_back();
      while (begin < trimmed.size()) {
         std::size_t end = trimmed.find('\n', begin);
         if (end == std::string::npos) end = trimmed.size();
         contents += (end == begin) ? std::string(" ") : trimmed.substr(begin, end - begin);
         contents += '\n';
         begin = end + 1;
      }
      contents += '\n';
      const int number_options = header.ampl_options[0];
      if (number_options > 0) {
         contents += "Options\n";
         integer(number_options + (header.ampl_options[2] == 3 ? 2 : 0));
         for (int k = 1; k <= number_options; ++k) integer(header.ampl_options[k]);
         integer(model.number_constraints);
         integer(y != nullptr ? model.number_constraints : 0);
         integer(model.number_variables);
         integer(x != nullptr ? model.number_variables : 0);
         if (header.ampl_options[2] == 3) number(header.variable_bound_tolerance);
      }
      if (y != nullptr) for (int i = 0; i < model.number_constraints; ++i) number(y[i]);
      if (x != nullptr) for (int j = 0; j < model.number_variables; ++j) number(x[j]);
      if (header.flags & 1) contents += "objno 0 " + std::to_string(solve_result_code) + "\n";
      // suffix sections (ASL's write_sol): "suffix kind n namelen tablen tablines", name, then "index value" lines
      for (const Suffix& suffix: output_suffixes) {
         std::size_t number_nonzeros = 0;
         for (const double value: suffix.values) number_nonzeros += (value != 0.);
         if (number_nonzeros == 0) continue;
         const int kind = static_cast<int>(suffix.target) | (suffix.is_real ? 4 : 0);
         contents += "suffix " + std::to_string(kind) + " " + std::to_string(number_nonzeros) + " " +
            std::to_string(suffix.name.size() + 1) + " 0 0\n" + suffix.name + "\n";
         for (std::size_t index = 0; index < suffix.values.size(); ++index) {
            if (suffix.values[index] == 0.) continue;
            contents += std::to_string(index) + " ";
            if (suffix.is_real) number(suffix.values[index]);
            else integer(static_cast<long long>(suffix.values[index]));
         }
      }
      std::ofstream file(file_name, std::ios::binary);
      if (!file) throw std::runtime_error("cppasl: cannot write " + file_name);
      file << contents;
   }

} // namespace cppasl
