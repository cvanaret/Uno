// Copyright (c) 2025 Stefano Lovato and Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <vector>
#include "MatlabModel.hpp"
#include "../unomex/unomex_function.hpp"
#include "../unomex/unomex_conversion.hpp"
#include "../unomex/unomex_validation.hpp"
#include "../unomex/unomex_utils.hpp"
#include "mex.h"

namespace uno {

    MatlabModel::MatlabModel(const MatlabUserModel& user_model):
            Model("Matlab model", static_cast<size_t>(user_model.number_variables), static_cast<size_t>(user_model.number_constraints),
            static_cast<double>(user_model.optimization_sense), user_model.lagrangian_sign_convention),
            user_model(user_model),
            equality_constraints_collection(this->equality_constraints),
            inequality_constraints_collection(this->inequality_constraints) {
        this->find_fixed_variables(this->fixed_variables);
        this->partition_constraints(this->equality_constraints, this->inequality_constraints);
    }

    ProblemType MatlabModel::get_problem_type() const {
        return this->user_model.problem_type;
    }

    bool MatlabModel::has_jacobian_operator() const {
        return this->user_model.jacobian_operator != nullptr;
    }

    bool MatlabModel::has_jacobian_transposed_operator() const {
        return this->user_model.jacobian_transposed_operator != nullptr;
    }

    bool MatlabModel::has_hessian_operator() const {
        return this->user_model.lagrangian_hessian_operator != nullptr;
    }

    bool MatlabModel::has_hessian_matrix() const {
        return this->user_model.lagrangian_hessian != nullptr;
    }

    double MatlabModel::evaluate_objective(const Vector<double>& x) const {
        double objective_value{0.};
        // objective_value = objective_function(x);
        if (this->user_model.objective_function) {
            const std::vector<mxArray*> inputs({vector_to_mxArray(x, this->number_variables)});
            std::vector<mxArray*> outputs(1);
            MxArrayVectorGuard input_guard(inputs);
            MxArrayVectorGuard output_guard(outputs);
            try {
                call_matlab_function(user_model.objective_function, inputs, outputs);
            } catch (const MatlabFunctionError& err) {
                mexWarnMsgIdAndTxt("uno:warning", "%s\n", err.what());
                throw FunctionEvaluationError();
            }
            if (std::string errmsg; !validate_double_scalar_output(outputs[0], errmsg)) {
                mexWarnMsgIdAndTxt("uno:warning", "Error in objective function.\n%s\n", errmsg.c_str());
                throw FunctionEvaluationError();
            }
            objective_value = this->optimization_sense * mxArray_to_scalar<double>(outputs[0]);
            ++this->number_model_evaluations.objective;
        }
        return objective_value;
    }

    void MatlabModel::evaluate_constraints(const Vector<double>& x, Vector<double>& constraints) const {
        // constraint = constraint_function(x);
        if (this->user_model.constraint_functions) {
            std::vector<mxArray*> inputs({vector_to_mxArray(x, this->number_variables)});
            std::vector<mxArray*> outputs(1);
            MxArrayVectorGuard input_guard(inputs);
            MxArrayVectorGuard output_guard(outputs);
            try {
                call_matlab_function(user_model.constraint_functions, inputs, outputs);
            } catch (const MatlabFunctionError& err) {
                mexWarnMsgIdAndTxt("uno:warning", "%s", err.what());
                throw FunctionEvaluationError();
            }
            if (std::string errmsg; !validate_double_vector_output(outputs[0], this->number_constraints, errmsg)) {
                mexWarnMsgIdAndTxt("uno:warning", "Error in constraint function.\n%s\n", errmsg.c_str());
                throw FunctionEvaluationError();
            }
            mxArray_to_vector(outputs[0], constraints);
            ++this->number_model_evaluations.constraints;
        }
    }

    void MatlabModel::evaluate_objective_gradient(const Vector<double>& x, Vector<double>& gradient) const {
        // gradient = objective_gradient(x);
        if (this->user_model.objective_gradient) {
            std::vector<mxArray*> inputs({vector_to_mxArray(x, this->number_variables)});
            std::vector<mxArray*> outputs(1);
            MxArrayVectorGuard input_guard(inputs);
            MxArrayVectorGuard output_guard(outputs);
            try {
                call_matlab_function(user_model.objective_gradient, inputs, outputs);
            } catch (const MatlabFunctionError& err) {
                mexWarnMsgIdAndTxt("uno:warning", "%s", err.what());
                throw GradientEvaluationError();
            }
            if (std::string errmsg; !validate_double_vector_output(outputs[0], this->number_variables, errmsg)) {
                mexWarnMsgIdAndTxt("uno:warning", "Error in objective gradient.\n%s\n", errmsg.c_str());
                throw FunctionEvaluationError();
            }
            mxArray_to_vector(outputs[0], gradient);
            for (size_t variable_index: Range(this->number_variables)) {
                gradient[variable_index] *= this->optimization_sense;
            }
            ++this->number_model_evaluations.objective_gradient;
        }
    }

    void MatlabModel::compute_constraint_jacobian_sparsity(int* row_indices, int* column_indices, int solver_indexing,
            MatrixOrder /*matrix_order*/) const {
        // copy the indices of the user sparsity patterns to the Uno vectors
        for (size_t nonzero_index: Range(static_cast<size_t>(this->user_model.number_jacobian_nonzeros))) {
            row_indices[nonzero_index] = this->user_model.jacobian_row_indices[nonzero_index];
            column_indices[nonzero_index] = this->user_model.jacobian_column_indices[nonzero_index];
        }
        // TODO matrix_order

        // handle the solver indexing
        if (this->user_model.base_indexing != solver_indexing) {
            const int indexing_difference = solver_indexing - this->user_model.base_indexing;
            for (size_t nonzero_index: Range(static_cast<size_t>(this->user_model.number_jacobian_nonzeros))) {
                row_indices[nonzero_index] += indexing_difference;
                column_indices[nonzero_index] += indexing_difference;
            }
        }
    }

    void MatlabModel::compute_hessian_sparsity(int* row_indices, int* column_indices, int solver_indexing) const {
        // copy the indices of the user sparsity patterns to the Uno vectors
        const size_t number_hessian_nonzeros = this->number_hessian_nonzeros();
        for (size_t nonzero_index: Range(number_hessian_nonzeros)) {
            row_indices[nonzero_index] = this->user_model.hessian_row_indices[nonzero_index];
            column_indices[nonzero_index] = this->user_model.hessian_column_indices[nonzero_index];
        }

        // handle the solver indexing
        if (this->user_model.base_indexing != solver_indexing) {
            const int indexing_difference = solver_indexing - this->user_model.base_indexing;
            for (size_t nonzero_index: Range(number_hessian_nonzeros)) {
                row_indices[nonzero_index] += indexing_difference;
                column_indices[nonzero_index] += indexing_difference;
            }
        }
    }
    
    void MatlabModel::evaluate_constraint_jacobian(const Vector<double>& x, double* jacobian_values) const {
        // jacobian_values = constraint_jacobian(x);
        if (this->user_model.constraint_jacobian) {
            std::vector<mxArray*> inputs({vector_to_mxArray(x, this->number_variables)});
            std::vector<mxArray*> outputs(1);
            MxArrayVectorGuard input_guard(inputs);
            MxArrayVectorGuard output_guard(outputs);
            try {
                call_matlab_function(user_model.constraint_jacobian, inputs, outputs);
            } catch (const MatlabFunctionError& err) {
                mexWarnMsgIdAndTxt("uno:warning", "%s", err.what());
                throw GradientEvaluationError();
            }
            if (std::string errmsg; !validate_double_vector_output(outputs[0], this->number_jacobian_nonzeros(), errmsg)) {
                mexWarnMsgIdAndTxt("uno:warning", "Error in constraint Jacobian.\n%s\n", errmsg.c_str());
                throw FunctionEvaluationError();
            }
            mxArray_to_pointer(outputs[0], jacobian_values);
            ++this->number_model_evaluations.jacobian;
        }
    }

    void MatlabModel::evaluate_lagrangian_hessian(const Vector<double>& x, double objective_multiplier, const Vector<double>& multipliers,
         double* hessian_values) const {
        // hessian_values = lagrangian_hessian(x, rho, y);
        if (this->user_model.lagrangian_hessian) {
            objective_multiplier *= this->optimization_sense;
            // if the model has a different sign convention for the Lagrangian than Uno, flip the signs of the multipliers
            if (this->user_model.lagrangian_sign_convention == UNO_MULTIPLIER_POSITIVE) {
                const_cast<Vector<double>&>(multipliers).scale(-1.);
            }
            // eval matlab function
            std::vector<mxArray*> inputs({vector_to_mxArray(x, this->number_variables), scalar_to_mxArray(objective_multiplier),vector_to_mxArray(multipliers)});
            std::vector<mxArray*> outputs(1);
            MxArrayVectorGuard input_guard(inputs);
            MxArrayVectorGuard output_guard(outputs);
            // flip the signs of the multipliers back
            if (this->user_model.lagrangian_sign_convention == UNO_MULTIPLIER_POSITIVE) {
                const_cast<Vector<double>&>(multipliers).scale(-1.);
            }
            try {
                call_matlab_function(user_model.lagrangian_hessian, inputs, outputs);
            } catch (const MatlabFunctionError& err) {
                mexWarnMsgIdAndTxt("uno:warning", "%s", err.what());
                throw HessianEvaluationError();
            }
            if (std::string errmsg; !validate_double_vector_output(outputs[0], this->number_hessian_nonzeros(), errmsg)) {
                mexWarnMsgIdAndTxt("uno:warning", "Error in Lagrangian Hessian.\n%s\n", errmsg.c_str());
                throw FunctionEvaluationError();
            }
            mxArray_to_pointer(outputs[0], hessian_values);
            ++this->number_model_evaluations.hessian;
        }
        else {
            throw std::runtime_error("evaluate_lagrangian_hessian not implemented");
        }
    }

    void MatlabModel::compute_jacobian_vector_product(const double* x, const double* vector, double* result) const {
        // result = jacobian_operator(x, v);
        if (this->user_model.jacobian_operator) {
            mxArray* x_arr = pointer_to_mxArray(x, this->number_variables);
            mxArray* vector_arr = pointer_to_mxArray(vector, this->number_variables);
            std::vector<mxArray*> inputs({x_arr, vector_arr});
            std::vector<mxArray*> outputs(1);
            MxArrayVectorGuard input_guard(inputs);
            MxArrayVectorGuard output_guard(outputs);
            try {
                call_matlab_function(user_model.jacobian_operator, inputs, outputs);
            } catch (const MatlabFunctionError& err) {
                mexWarnMsgIdAndTxt("uno:warning", "%s", err.what());
                throw GradientEvaluationError();
            }
            if (std::string errmsg; !validate_double_vector_output(outputs[0], this->number_constraints, errmsg)) {
                mexWarnMsgIdAndTxt("uno:warning", "Error in Jacobian operator.\n%s\n", errmsg.c_str());
                throw FunctionEvaluationError();
            }
            mxArray_to_pointer(outputs[0], result);
        }
        else {
            throw std::runtime_error("compute_jacobian_vector_product not implemented");
        }
    }

    void MatlabModel::compute_jacobian_transposed_vector_product(const double* x, const double* vector, double* result) const {
        // result = jacobian_transposed_operator(x, v);
        if ((this->user_model.jacobian_operator != nullptr) &&
            !mxIsEmpty(this->user_model.jacobian_operator)) {
            mxArray* x_arr = pointer_to_mxArray(x, this->number_variables);
            mxArray* vector_arr = pointer_to_mxArray(vector, this->number_constraints);
            std::vector<mxArray*> inputs({x_arr, vector_arr});
            std::vector<mxArray*> outputs(1);
            MxArrayVectorGuard input_guard(inputs);
            MxArrayVectorGuard output_guard(outputs);
            try {
                call_matlab_function(user_model.jacobian_transposed_operator, inputs, outputs);
            } catch (const MatlabFunctionError& err) {
                mexWarnMsgIdAndTxt("uno:warning", "%s", err.what());
                throw GradientEvaluationError();
            }
            if (std::string errmsg; !validate_double_vector_output(outputs[0], this->number_variables, errmsg)) {
                mexWarnMsgIdAndTxt("uno:warning", "Error in Jacobian transposed operator.\n%s\n", errmsg.c_str());
                throw FunctionEvaluationError();
            }
            mxArray_to_pointer(outputs[0], result);
        }
        else {
            throw std::runtime_error("compute_jacobian_transposed_vector_product not implemented");
        }
    }

    void MatlabModel::compute_hessian_vector_product(const double* x, const double* vector, double objective_multiplier, const Vector<double>& multipliers,
         double* result) const {
        // result = lagrangian_hessian_operator(x, rho, y, v);
        if (this->user_model.lagrangian_hessian_operator) {
            objective_multiplier *= this->optimization_sense;
            // if the model has a different sign convention for the Lagrangian than Uno, flip the signs of the multipliers
            if (this->user_model.lagrangian_sign_convention == UNO_MULTIPLIER_POSITIVE) {
                const_cast<Vector<double>&>(multipliers).scale(-1.);
            }
            mxArray* x_arr = pointer_to_mxArray(x, this->number_variables);
            mxArray* vector_arr = pointer_to_mxArray(vector, this->number_variables);
            std::vector<mxArray*> inputs({x_arr, scalar_to_mxArray(objective_multiplier), vector_to_mxArray(multipliers), vector_arr});
            std::vector<mxArray*> outputs(1);
            MxArrayVectorGuard input_guard(inputs);
            MxArrayVectorGuard output_guard(outputs);
            // flip the signs of the multipliers back
            if (this->user_model.lagrangian_sign_convention == UNO_MULTIPLIER_POSITIVE) {
                const_cast<Vector<double>&>(multipliers).scale(-1.);
            }
            try {
                call_matlab_function(user_model.lagrangian_hessian_operator, inputs, outputs);
            } catch (const MatlabFunctionError& err) {
                mexWarnMsgIdAndTxt("uno:warning", "%s", err.what());
                throw HessianEvaluationError();
            }
            if (std::string errmsg; !validate_double_vector_output(outputs[0], this->number_variables, errmsg)) {
                mexWarnMsgIdAndTxt("uno:warning", "Error in Lagrangian Hessian operator.\n%s\n", errmsg.c_str());
                throw FunctionEvaluationError();
            }
            mxArray_to_pointer(outputs[0], result);
        } 
        else {
            throw std::runtime_error("compute_hessian_vector_product not implemented");
        }
    }

    double MatlabModel::variable_lower_bound(size_t variable_index) const {
        Vector<double> variables_lower_bounds = mxArray_to_vector<double>(this->user_model.variables_lower_bounds);
        return variables_lower_bounds[variable_index];
    }

    double MatlabModel::variable_upper_bound(size_t variable_index) const {
        Vector<double> variables_upper_bounds = mxArray_to_vector<double>(this->user_model.variables_upper_bounds);
        return variables_upper_bounds[variable_index];
    }

    const SparseVector<size_t>& MatlabModel::get_slacks() const {
        return this->slacks;
    }

    const Vector<size_t>& MatlabModel::get_fixed_variables() const {
        return this->fixed_variables;
    }

    double MatlabModel::constraint_lower_bound(size_t constraint_index) const {
        Vector<double> constraints_lower_bounds = mxArray_to_vector<double>(this->user_model.constraints_lower_bounds);
        return constraints_lower_bounds[constraint_index];
    }

    double MatlabModel::constraint_upper_bound(size_t constraint_index) const {
        Vector<double> constraints_upper_bounds = mxArray_to_vector<double>(this->user_model.constraints_upper_bounds);
        return constraints_upper_bounds[constraint_index];
    }

    const Collection<size_t>& MatlabModel::get_equality_constraints() const {
        return this->equality_constraints_collection;
    }

    const Collection<size_t>& MatlabModel::get_inequality_constraints() const {
        return this->inequality_constraints_collection;
    }

    const Collection<size_t>& MatlabModel::get_linear_constraints() const {
        return this->linear_constraints;
    }

    void MatlabModel::initial_primal_point(Vector<double>& x) const {
        x.fill(0);
        if (this->user_model.initial_primal_iterate) {
            mxArray_to_vector(this->user_model.initial_primal_iterate, x);
        }
    }

    void MatlabModel::initial_dual_point(Vector<double>& multipliers) const {
        multipliers.fill(0);
        if (this->user_model.initial_dual_iterate) {
            mxArray_to_vector(this->user_model.initial_dual_iterate, multipliers);
        }
        if (this->user_model.lagrangian_sign_convention == UNO_MULTIPLIER_POSITIVE) {
            multipliers.scale(-1.);
        }
    }
    
    void MatlabModel::postprocess_solution(Iterate& /* iterate */) const {
        // do nothing
    }

    size_t MatlabModel::number_jacobian_nonzeros() const {
        return static_cast<size_t>(this->user_model.number_jacobian_nonzeros);
    }

    size_t MatlabModel::number_hessian_nonzeros() const {
        if (this->user_model.number_hessian_nonzeros.has_value()) {
            return static_cast<size_t>(*this->user_model.number_hessian_nonzeros);
        }
        else {
            throw std::runtime_error("The number of Hessian nonzeros is not available in MatlabModel");
        }
    }

    size_t MatlabModel::number_model_objective_evaluations() const {
      return this->number_model_evaluations.objective;
   }

   size_t MatlabModel::number_model_constraints_evaluations() const {
      return this->number_model_evaluations.constraints;
   }

   size_t MatlabModel::number_model_objective_gradient_evaluations() const {
      return this->number_model_evaluations.objective_gradient;
   }

   size_t MatlabModel::number_model_jacobian_evaluations() const {
      return this->number_model_evaluations.jacobian;
   }

   size_t MatlabModel::number_model_hessian_evaluations() const {
      return this->number_model_evaluations.hessian;
   }

   void MatlabModel::reset_number_evaluations() const {
      this->number_model_evaluations.reset();
   }
    
} // namespace


