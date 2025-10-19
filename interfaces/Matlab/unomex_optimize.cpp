// Copyright (c) 2025 Stefano Lovato and Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <stdint.h>
#include <algorithm>
#include "../UserModel.hpp"
#include "Uno.hpp"
#include "linear_algebra/SparseVector.hpp"
#include "linear_algebra/Vector.hpp"
#include "model/Model.hpp"
#include "options/DefaultOptions.hpp"
#include "options/Presets.hpp"
#include "optimization/EvaluationErrors.hpp"
#include "optimization/Iterate.hpp"
#include "symbolic/CollectionAdapter.hpp"
#include "symbolic/Range.hpp"
#include "tools/Infinity.hpp"
#include "tools/Logger.hpp"
#include "tools/NumberModelEvaluations.hpp"
#include "tools/UserCallbacks.hpp"
#include "cpp_classes/MxStruct.hpp"
#include "cpp_classes/ErrorString.hpp"
#include "unomex/unomex_function.hpp"
#include "unomex/unomex_conversion.hpp"
#include "unomex/unomex_validation.hpp"
#include "unomex/unomex_utils.hpp"
#include "mex.h"

using namespace uno;

using MatlabUserModel = UserModel<handle_t, handle_t, handle_t,handle_t, handle_t,handle_t,
   handle_t, handle_t, mxArray*, void*>; // user data are handled in matlab wrapper

// UnoModel contains an instance of UserModel and complies with the Model interface
class UnoModel: public Model {
public:
    explicit UnoModel(const MatlabUserModel& user_model):
            Model("Matlab model", static_cast<size_t>(user_model.number_variables), static_cast<size_t>(user_model.number_constraints),
            static_cast<double>(user_model.optimization_sense)),
            user_model(user_model),
            equality_constraints_collection(this->equality_constraints),
            inequality_constraints_collection(this->inequality_constraints) {
        this->find_fixed_variables(this->fixed_variables);
        this->partition_constraints(this->equality_constraints, this->inequality_constraints);
    }

    // availability of linear operators
    [[nodiscard]] bool has_jacobian_operator() const override {
        return this->user_model.jacobian_operator != nullptr;
    }

    [[nodiscard]] bool has_jacobian_transposed_operator() const override {
        return this->user_model.jacobian_transposed_operator != nullptr;
    }

    [[nodiscard]] bool has_hessian_operator() const override {
        return this->user_model.lagrangian_hessian_operator != nullptr;
    }

    [[nodiscard]] bool has_hessian_matrix() const override {
        return this->user_model.lagrangian_hessian != nullptr;
    }

    // function evaluations
    [[nodiscard]] double evaluate_objective(const Vector<double>& x) const override {
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

    void evaluate_constraints(const Vector<double>& x, Vector<double>& constraints) const override {
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

    // dense objective gradient
    void evaluate_objective_gradient(const Vector<double>& x, Vector<double>& gradient) const override {
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

    // sparsity patterns of Jacobian and Hessian
    void compute_constraint_jacobian_sparsity(int* row_indices, int* column_indices, int solver_indexing,
            MatrixOrder /*matrix_order*/) const override {
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

    void compute_hessian_sparsity(int* row_indices, int* column_indices, int solver_indexing) const override {
        // copy the indices of the user sparsity patterns to the Uno vectors
        for (size_t nonzero_index: Range(static_cast<size_t>(this->user_model.number_hessian_nonzeros))) {
            row_indices[nonzero_index] = this->user_model.hessian_row_indices[nonzero_index];
            column_indices[nonzero_index] = this->user_model.hessian_column_indices[nonzero_index];
        }

        // handle the solver indexing
        if (this->user_model.base_indexing != solver_indexing) {
            const int indexing_difference = solver_indexing - this->user_model.base_indexing;
            for (size_t nonzero_index: Range(static_cast<size_t>(this->user_model.number_hessian_nonzeros))) {
                row_indices[nonzero_index] += indexing_difference;
                column_indices[nonzero_index] += indexing_difference;
            }
        }
    }
    
    // numerical evaluations of Jacobian and Hessian
    void evaluate_constraint_jacobian(const Vector<double>& x, double* jacobian_values) const override {
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

    void evaluate_lagrangian_hessian(const Vector<double>& x, double objective_multiplier, const Vector<double>& multipliers,
         double* hessian_values) const override {
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

    void compute_jacobian_vector_product(const double* x, const double* vector, double* result) const override {
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

    void compute_jacobian_transposed_vector_product(const double* x, const double* vector, double* result) const override {
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

    void compute_hessian_vector_product(const double* x, const double* vector, double objective_multiplier, const Vector<double>& multipliers,
         double* result) const override {
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

    [[nodiscard]] double variable_lower_bound(size_t variable_index) const override {
        Vector<double> variables_lower_bounds = mxArray_to_vector<double>(this->user_model.variables_lower_bounds);
        return variables_lower_bounds[variable_index];
    }

    [[nodiscard]] double variable_upper_bound(size_t variable_index) const override {
        Vector<double> variables_upper_bounds = mxArray_to_vector<double>(this->user_model.variables_upper_bounds);
        return variables_upper_bounds[variable_index];
    }

    [[nodiscard]] const SparseVector<size_t>& get_slacks() const override {
        return this->slacks;
    }

    [[nodiscard]] const Vector<size_t>& get_fixed_variables() const override {
        return this->fixed_variables;
    }

    [[nodiscard]] double constraint_lower_bound(size_t constraint_index) const override {
        Vector<double> constraints_lower_bounds = mxArray_to_vector<double>(this->user_model.constraints_lower_bounds);
        return constraints_lower_bounds[constraint_index];
    }

    [[nodiscard]] double constraint_upper_bound(size_t constraint_index) const override {
        Vector<double> constraints_upper_bounds = mxArray_to_vector<double>(this->user_model.constraints_upper_bounds);
        return constraints_upper_bounds[constraint_index];
    }

    [[nodiscard]] const Collection<size_t>& get_equality_constraints() const override {
        return this->equality_constraints_collection;
    }

    [[nodiscard]] const Collection<size_t>& get_inequality_constraints() const override {
        return this->inequality_constraints_collection;
    }

    [[nodiscard]] const Collection<size_t>& get_linear_constraints() const override {
        return this->linear_constraints;
    }

    void initial_primal_point(Vector<double>& x) const override {
        x.fill(0);
        if (this->user_model.initial_primal_iterate) {
            mxArray_to_vector(this->user_model.initial_primal_iterate, x);
        }
    }

    void initial_dual_point(Vector<double>& multipliers) const override {
        multipliers.fill(0);
        if (this->user_model.initial_dual_iterate) {
            mxArray_to_vector(this->user_model.initial_dual_iterate, multipliers);
        }
        if (this->user_model.lagrangian_sign_convention == UNO_MULTIPLIER_POSITIVE) {
            multipliers.scale(-1.);
        }
    }
    
    void postprocess_solution(Iterate& iterate) const override {
        // flip the signs of the multipliers, depending on what the sign convention of the Lagrangian is, and whether
        // we maximize
        iterate.multipliers.constraints *= -this->user_model.lagrangian_sign_convention * this->optimization_sense;
        iterate.multipliers.lower_bounds *= -this->user_model.lagrangian_sign_convention * this->optimization_sense;
        iterate.multipliers.upper_bounds *= -this->user_model.lagrangian_sign_convention * this->optimization_sense;
        iterate.evaluations.objective *= this->optimization_sense;
    }

    [[nodiscard]] size_t number_jacobian_nonzeros() const override {
        return static_cast<size_t>(this->user_model.number_jacobian_nonzeros);
    }

    [[nodiscard]] size_t number_hessian_nonzeros() const override {
        return static_cast<size_t>(this->user_model.number_hessian_nonzeros);
    }

    [[nodiscard]] size_t number_model_objective_evaluations() const override {
      return this->number_model_evaluations.objective;
   }

   [[nodiscard]] size_t number_model_constraints_evaluations() const override {
      return this->number_model_evaluations.constraints;
   }

   [[nodiscard]] size_t number_model_objective_gradient_evaluations() const override {
      return this->number_model_evaluations.objective_gradient;
   }

   [[nodiscard]] size_t number_model_jacobian_evaluations() const override {
      return this->number_model_evaluations.jacobian;
   }

   [[nodiscard]] size_t number_model_hessian_evaluations() const override {
      return this->number_model_evaluations.hessian;
   }

   void reset_number_evaluations() const override {
      this->number_model_evaluations.reset();
   }

protected:
   const MatlabUserModel& user_model;
   mutable NumberModelEvaluations number_model_evaluations;
   const SparseVector<size_t> slacks{};
   Vector<size_t> fixed_variables{};
   const ForwardRange linear_constraints{0};
   std::vector<size_t> equality_constraints;
   CollectionAdapter<std::vector<size_t>> equality_constraints_collection;
   std::vector<size_t> inequality_constraints;
   CollectionAdapter<std::vector<size_t>> inequality_constraints_collection;

};

class MatlabStreamCallback : public UserStreamCallback {
public: 
    MatlabStreamCallback(handle_t logger_stream_callback) :
    UserStreamCallback(), logger_stream_callback(logger_stream_callback) { }
    ~MatlabStreamCallback() override { }

int32_t operator()(const char* buf, int32_t len) const override {
    if (this->logger_stream_callback) {
        // logger_stream_callback(str)
        std::string str(buf, len);
        std::vector<mxArray*> inputs = { string_to_mxArray(str) };
        std::vector<mxArray*> outputs(0);
        MxArrayVectorGuard input_guard(inputs);
        // MxArrayVectorGuard output_guard(outputs);
        try {
            call_matlab_function(this->logger_stream_callback, inputs, outputs);
        }
        catch (const MatlabFunctionError& err) {
            mexWarnMsgIdAndTxt("uno:warning", "%s", err.what());
        }
        // no output to validate
        mexEvalString("pause(0);"); // force output to appear immediately 
        return len;
    } 
    else {
        // call mexPrintf
        mexPrintf("%.*s", static_cast<int>(len), buf);
        mexEvalString("pause(0);"); // force output to appear immediately 
        return len;
    }
}

private:
    handle_t logger_stream_callback;
};

#ifdef HAS_MXLIBUT
// These functions are to allow interrupt from MATLAB using CTRL+C (see 
// https://undocumentedmatlab.com/articles/mex-ctrl-c-interrupt) 
extern "C" bool utIsInterruptPending();
extern "C" bool utSetInterruptPending(bool);
#else
// libut is not available, use dummy
constexpr bool utIsInterruptPending() { return false; }
constexpr bool utSetInterruptPending(bool) { return false; }
#endif

// Matlab user callbacks
class MatlabUserCallbacks : public UserCallbacks {
public:
    MatlabUserCallbacks(handle_t notify_acceptable_iterate_callback, handle_t user_termination_callback) : UserCallbacks(),
        notify_acceptable_iterate_callback(notify_acceptable_iterate_callback),
        user_termination_callback(user_termination_callback) { }

    void notify_acceptable_iterate(const Vector<double>& primals, const Multipliers& multipliers, double objective_multiplier, double primal_feasibility, double stationarity, double complementarity) override {
        // notify_acceptable_iterate_callback(x, yl, yb, y, rho, feas, stat, compl);
        if (this->notify_acceptable_iterate_callback) {
            std::vector<mxArray*> inputs = {vector_to_mxArray(primals), vector_to_mxArray(multipliers.lower_bounds), vector_to_mxArray(multipliers.upper_bounds), vector_to_mxArray(multipliers.constraints), scalar_to_mxArray(objective_multiplier), scalar_to_mxArray(primal_feasibility), scalar_to_mxArray(stationarity), scalar_to_mxArray(complementarity)};
            std::vector<mxArray*> outputs(0);
            MxArrayVectorGuard input_guard(inputs);
            // MxArrayVectorGuard output_guard(outputs);
            try {
                call_matlab_function(this->notify_acceptable_iterate_callback, inputs, outputs);
            }
            catch (const MatlabFunctionError& err) {
                mexWarnMsgIdAndTxt("uno:warning", "%s", err.what());
            }
            // no output to validate
        }
    }
    bool user_termination(const Vector<double>& primals, const Multipliers& multipliers, double objective_multiplier, double primal_feasibility, double stationarity, double complementarity) override {
        // handle matlab CTRL+C event
        if (utIsInterruptPending()) {
            utSetInterruptPending(false);
            return true;
        }
        // terminate = user_termination_callback(x, yl, yb, y, rho, feas, stat, compl);
        if (this->user_termination_callback) {
            std::vector<mxArray*> inputs = {vector_to_mxArray(primals), vector_to_mxArray(multipliers.lower_bounds), vector_to_mxArray(multipliers.upper_bounds), vector_to_mxArray(multipliers.constraints), scalar_to_mxArray(objective_multiplier), scalar_to_mxArray(primal_feasibility), scalar_to_mxArray(stationarity), scalar_to_mxArray(complementarity)};
            std::vector<mxArray*> outputs(1);
            MxArrayVectorGuard input_guard(inputs);
            MxArrayVectorGuard output_guard(outputs);
            try {
                call_matlab_function(this->user_termination_callback, inputs, outputs);
            }  
            catch (const MatlabFunctionError& err) {
                mexWarnMsgIdAndTxt("uno:warning", "%s", err.what());
                return false;
            }
            if (std::string errmsg; !validate_bool_scalar_output(outputs[0], errmsg)) {
                mexWarnMsgIdAndTxt("uno:warning", "Error in user termination callback.\n%s\n", errmsg.c_str());
                return false;
            }
            return mxArray_to_scalar<bool>(outputs[0]);
        } else {
            return false; // never terminate
        }
    }

private:
    handle_t notify_acceptable_iterate_callback;
    handle_t user_termination_callback;
};

// gateway function: result = uno_optimize(model[, options, callbacks])
void mexFunction( int /* nlhs */, mxArray* plhs[], int nrhs, const mxArray* prhs[]) {
    std::string errmsg; // declare once for all
    // validate argument list
    if (nrhs < 1) {
        errmsg = ErrorString::format_error(ErrorType::NARGIN_NOTENOUGH, 1);
        mexErrMsgIdAndTxt("uno:error", errmsg.c_str());
    }
    if (nrhs > 3) {
        errmsg = ErrorString::format_error(ErrorType::NARGIN_TOOMANY);
        mexErrMsgIdAndTxt("uno:error", errmsg.c_str());
    }

    // model (mandatory)
    if (!validate_struct_input(prhs[0], 1, errmsg)) {
        mexErrMsgIdAndTxt("uno:error", errmsg.c_str());
    }
    const MxStruct model = mxArray_to_mxStruct(prhs[0]);

    // options (optional)
    MxStruct options;
    if (nrhs > 1) {
        if (!validate_struct_input(prhs[1], 2, errmsg)) {
            mexErrMsgIdAndTxt("uno:error", errmsg.c_str());
        }
        options = mxArray_to_mxStruct(prhs[1]);
    }

    // callbacks (optional)
    MxStruct callbacks;
    if (nrhs > 2) {
        if (!validate_struct_input(prhs[2], 3, errmsg)) {
            mexErrMsgIdAndTxt("uno:error", errmsg.c_str());
        }
        callbacks = mxArray_to_mxStruct(prhs[2]);
    }

    // problem type (mandatory)
    mxArray* problem_type = model["problem_type"];
    if (!validate_char_field(problem_type, {UNO_PROBLEM_LINEAR, UNO_PROBLEM_QUADRATIC, UNO_PROBLEM_NONLINEAR}, "problem_type", errmsg)) {
        mexErrMsgIdAndTxt("uno:error", "Invalid model. %s", errmsg.c_str());
    }

    // base indexing (mandatory)
    mxArray* base_indexing = model["base_indexing"];
    if (!validate_positive_integer_field(base_indexing, "base_indexing", errmsg)) {
        mexErrMsgIdAndTxt("uno:error", "Invalid model. %s", errmsg.c_str());
    }

    // variables (mandatory)
    mxArray* number_variables = model["number_variables"];
    mxArray* variables_lower_bounds = model["variables_lower_bounds"];
    mxArray* variables_upper_bounds = model["variables_upper_bounds"];
    if (!validate_positive_integer_field(number_variables, "number_variables", errmsg)) {
        mexErrMsgIdAndTxt("uno:error", "Invalid model. %s", errmsg.c_str());
    }
    if (!validate_double_vector_field(variables_lower_bounds, static_cast<size_t>(mxArray_to_scalar<double>(number_variables)), "variables_lower_bounds", errmsg)) {
        mexErrMsgIdAndTxt("uno:error", "Invalid model. %s", errmsg.c_str());
    }
    if (!validate_double_vector_field(variables_upper_bounds, static_cast<size_t>(mxArray_to_scalar<double>(number_variables)), "variables_lower_bounds", errmsg)) {
        mexErrMsgIdAndTxt("uno:error", "Invalid model. %s", errmsg.c_str());
    }

    // objective (mandatory)
    mxArray* objective_function = model["objective_function"];
    mxArray* objective_gradient = model["objective_gradient"];
    mxArray* optimization_sense = model["optimization_sense"];
    // objective = objective_function(x)
    if (!validate_matlab_handle(objective_function, "objective_function", errmsg)) {
        mexErrMsgIdAndTxt("uno:error", "Invalid model objective function. %s", errmsg.c_str());
    }
    // gradient = objective_gradient(x)
    if (!validate_matlab_handle(objective_gradient, "objective_gradient", errmsg)) {
        mexErrMsgIdAndTxt("uno:error", "Invalid model objective gradient. %s", errmsg.c_str());
    }
    if (!validate_unitary_field(optimization_sense, "optimization_sense", errmsg)) {
        mexErrMsgIdAndTxt("uno:error", "Invalid model. %s", errmsg.c_str());
    }

    // constraints (mandatory)
    mxArray* number_constraints = model["number_constraints"];
    mxArray* constraints_lower_bounds = model["constraints_lower_bounds"];
    mxArray* constraints_upper_bounds = model["constraints_upper_bounds"];
    mxArray* constraint_function = model["constraint_function"];
    if (!validate_positive_integer_field(number_constraints, "number_constraints", errmsg)) {
        mexErrMsgIdAndTxt("uno:error", "Invalid model. %s", errmsg.c_str());
    }
    if (!validate_double_vector_field(constraints_lower_bounds, static_cast<size_t>(mxArray_to_scalar<double>(number_constraints)), "constraints_lower_bounds", errmsg)) {
        mexErrMsgIdAndTxt("uno:error", "Invalid model. %s", errmsg.c_str());
    }
    if (!validate_double_vector_field(constraints_upper_bounds, static_cast<size_t>(mxArray_to_scalar<double>(number_constraints)), "constraints_upper_bounds", errmsg)) {
        mexErrMsgIdAndTxt("uno:error", "Invalid model. %s", errmsg.c_str());
    }
    // constraints = constrain_function(x)
    if (!validate_matlab_handle(constraint_function, "constraint_function", errmsg)) {
        mexErrMsgIdAndTxt("uno:error", "Invalid model constraint function. %s", errmsg.c_str());
    }

    // jacobian (mandatory)
    mxArray* number_jacobian_nonzeros = model["number_jacobian_nonzeros"];
    mxArray* jacobian_row_indices = model["jacobian_row_indices"];
    mxArray* jacobian_column_indices = model["jacobian_column_indices"];
    mxArray* constraint_jacobian = model["constraint_jacobian"];
    if (!validate_positive_integer_field(number_jacobian_nonzeros, "number_jacobian_nonzeros", errmsg)) {
        mexErrMsgIdAndTxt("uno:error", "Invalid model. %s", errmsg.c_str());
    }
    if (!validate_integer_vector_field(jacobian_row_indices, static_cast<size_t>(mxArray_to_scalar<double>(number_jacobian_nonzeros)), "jacobian_row_indices", errmsg)) {
        mexErrMsgIdAndTxt("uno:error", "Invalid model. %s", errmsg.c_str());
    }
    if (!validate_integer_vector_field(jacobian_column_indices, static_cast<size_t>(mxArray_to_scalar<double>(number_jacobian_nonzeros)), "jacobian_column_indices", errmsg)) {
        mexErrMsgIdAndTxt("uno:error", "Invalid model. %s", errmsg.c_str());
    }
    // jacobian = constraint_jacobian(x)
    if (!validate_matlab_handle(constraint_jacobian, "constraint_jacobian", errmsg)) {
        mexErrMsgIdAndTxt("uno:error", "Invalid model constraint Jacobian. %s", errmsg.c_str());
    }

    // jacobian operators (optional)
    mxArray* jacobian_operator = model["jacobian_operator"];
    mxArray* jacobian_transposed_operator = model["jacobian_transposed_operator"];
    if (jacobian_operator != nullptr) {
        // result = jacobian_operator(x,v)
        if (!validate_matlab_handle(jacobian_operator, "jacobian_operator", errmsg)) {
            mexErrMsgIdAndTxt("uno:error", "Invalid model Jacobian operator. %s", errmsg.c_str());
        }
    }
    if (jacobian_transposed_operator != nullptr) {
        // jacobian_transposed_operator(x,v)
        if (!validate_matlab_handle(jacobian_transposed_operator, "jacobian_transposed_operator", errmsg)) {
            mexErrMsgIdAndTxt("uno:error", "Invalid model Jacobian transposed operator. %s", errmsg.c_str());
        }
    }

    // lagrangian sign convention (mandatory)
    mxArray* lagrangian_sign_convention = model["lagrangian_sign_convention"];
    if (!validate_unitary_field(lagrangian_sign_convention, "lagrangian_sign_convention", errmsg)) {
        mexErrMsgIdAndTxt("uno:error", "Invalid model. %s", errmsg.c_str());
    }

    // hessian (optional)
    mxArray* hessian_triangular_part = model["hessian_triangular_part"];
    mxArray* number_hessian_nonzeros = model["number_hessian_nonzeros"];
    mxArray* hessian_row_indices = model["hessian_row_indices"];
    mxArray* hessian_column_indices = model["hessian_column_indices"];
    mxArray* lagrangian_hessian = model["lagrangian_hessian"];
    if (lagrangian_hessian != nullptr) {
        if (!validate_char_field(hessian_triangular_part, {UNO_LOWER_TRIANGLE, UNO_UPPER_TRIANGLE}, "hessian_triangular_part", errmsg)) {
            mexErrMsgIdAndTxt("uno:error", "Invalid model. %s", errmsg.c_str());
        }
        if (!validate_positive_integer_field(number_hessian_nonzeros, "number_hessian_nonzeros", errmsg)) {
            mexErrMsgIdAndTxt("uno:error", "Invalid model. %s", errmsg.c_str());
        }
        if (!validate_integer_vector_field(hessian_row_indices, static_cast<size_t>(mxArray_to_scalar<double>(number_hessian_nonzeros)), "hessian_row_indices", errmsg)) {
            mexErrMsgIdAndTxt("uno:error", "Invalid model. %s", errmsg.c_str());
        }
        if (!validate_integer_vector_field(hessian_column_indices, static_cast<size_t>(mxArray_to_scalar<double>(number_hessian_nonzeros)), "hessian_column_indices", errmsg)) {
            mexErrMsgIdAndTxt("uno:error", "Invalid model. %s", errmsg.c_str());
        }
        // hessian = lagrangian_hessian(x,rho,y)
        if (!validate_matlab_handle(lagrangian_hessian, "lagrangian_hessian", errmsg)) {
            mexErrMsgIdAndTxt("uno:error", "Invalid model Lagrangian Hessian. %s", errmsg.c_str());
        }
    }

    // hessian operators (optional)
    mxArray* lagrangian_hessian_operator = model["lagrangian_hessian_operator"];
    if (lagrangian_hessian_operator != nullptr) {
        // lagrangian_hessian_operator(x,rho,y,v)
        if (!validate_matlab_handle(lagrangian_hessian_operator, "lagrangian_hessian_operator", errmsg)) {
            mexErrMsgIdAndTxt("uno:error", "Invalid model Lagrangian Hessian operator. %s", errmsg.c_str());
        }
    }

    // initial iterates (optional)
    mxArray* initial_primal_iterate = model["initial_primal_iterate"];
    mxArray* initial_dual_iterate = model["initial_dual_iterate"];
    if (initial_primal_iterate != nullptr) {
        if (!validate_double_vector_field(initial_primal_iterate, static_cast<size_t>(mxArray_to_scalar<double>(number_variables)), "initial_primal_iterate", errmsg)) {
            mexErrMsgIdAndTxt("uno:error", "Invalid model. %s", errmsg.c_str());
        }
    }
    if (initial_dual_iterate != nullptr) {
        if (!validate_double_vector_field(initial_dual_iterate, static_cast<size_t>(mxArray_to_scalar<double>(number_constraints)), "initial_dual_iterate", errmsg)) {
            mexErrMsgIdAndTxt("uno:error", "Invalid model. %s", errmsg.c_str());
        }
    }

    // callbacks (optional)
    mxArray* logger_stream_callback = callbacks["logger_stream_callback"];
    mxArray* notify_acceptable_iterate_callback = callbacks["notify_acceptable_iterate_callback"];
    mxArray* user_termination_callback = callbacks["user_termination_callback"];
    if (logger_stream_callback != nullptr) {
        // logger_stream_callback(str)
        if (!validate_matlab_handle(logger_stream_callback, "logger_stream_callback", errmsg)) {
            mexErrMsgIdAndTxt("uno:error", "Invalid logger stream callback. %s", errmsg.c_str());
        }
    }
    if (notify_acceptable_iterate_callback != nullptr) {
        // notify_acceptable_iterate_callback(x, yl, yb, y, rho, feas, stat, compl)
        if (!validate_matlab_handle(notify_acceptable_iterate_callback, "notify_acceptable_iterate_callback", errmsg)) {
            mexErrMsgIdAndTxt("uno:error", "Invalid notify acceptable iterate callback. %s", errmsg.c_str());
        }
    }
    if (user_termination_callback != nullptr) {
        // terminate = user_termination_callback(x, yl, yb, y, rho, feas, stat, compl)
        if (!validate_matlab_handle(user_termination_callback, "user_termination_callback", errmsg)) {
            mexErrMsgIdAndTxt("uno:error", "Invalid user termination callback. %s", errmsg.c_str());
        }
    }
    
    // set the logger stream
    MatlabStreamCallback matlab_stream_callback(logger_stream_callback);
    UserOStream matlab_ostream(&matlab_stream_callback);
    Logger::set_stream(matlab_ostream);
    
    // create user callbacks
    MatlabUserCallbacks user_callbacks(notify_acceptable_iterate_callback, user_termination_callback);

    // create user model
    MatlabUserModel user_model(mxArray_to_scalar<char>(problem_type), 
        static_cast<int32_t>(mxArray_to_scalar<double>(number_variables)), 
        static_cast<int32_t>(mxArray_to_scalar<double>(base_indexing)));
    user_model.number_constraints = static_cast<int32_t>(mxArray_to_scalar<double>(number_constraints));
    user_model.variables_lower_bounds = variables_lower_bounds;
    user_model.variables_upper_bounds = variables_upper_bounds;
    // objective
    user_model.objective_function = objective_function;
    user_model.objective_gradient = objective_gradient;
    user_model.optimization_sense = static_cast<int32_t>(mxArray_to_scalar<double>(optimization_sense));
    // constraint
    user_model.constraint_functions = constraint_function;
    user_model.constraints_lower_bounds = constraints_lower_bounds;
    user_model.constraints_upper_bounds = constraints_upper_bounds;
    user_model.constraint_jacobian = constraint_jacobian;
    user_model.number_jacobian_nonzeros = static_cast<int32_t>(mxArray_to_scalar<double>(number_jacobian_nonzeros));
    Vector jacobian_row_indices_vector = convert_vector_type<int32_t>( mxArray_to_vector<double>(jacobian_row_indices) );
    Vector jacobian_column_indices_vector = convert_vector_type<int32_t>( mxArray_to_vector<double>(jacobian_column_indices) );
    user_model.jacobian_row_indices = std::vector(jacobian_row_indices_vector.begin(), jacobian_row_indices_vector.end());
    user_model.jacobian_column_indices = std::vector(jacobian_column_indices_vector.begin(), jacobian_column_indices_vector.end());
    // jacobian operators
    user_model.jacobian_operator = jacobian_operator;
    user_model.jacobian_transposed_operator = jacobian_transposed_operator;
    // hessian
    user_model.lagrangian_sign_convention = static_cast<int32_t>(mxArray_to_scalar<double>(lagrangian_sign_convention));
    user_model.hessian_triangular_part = mxArray_to_scalar<char>(hessian_triangular_part);
    user_model.lagrangian_hessian = lagrangian_hessian;
    user_model.number_hessian_nonzeros = static_cast<int32_t>(mxArray_to_scalar<double>(number_hessian_nonzeros));
    Vector hessian_row_indices_vector = convert_vector_type<int32_t>( mxArray_to_vector<double>(hessian_row_indices) );
    Vector hessian_column_indices_vector = convert_vector_type<int32_t>( mxArray_to_vector<double>(hessian_column_indices) );
    user_model.hessian_row_indices = std::vector(hessian_row_indices_vector.begin(), hessian_row_indices_vector.end());
    user_model.hessian_column_indices = std::vector(hessian_column_indices_vector.begin(), hessian_column_indices_vector.end());
    // force lower triangular hessian
    if (user_model.lagrangian_hessian) {
        if (user_model.hessian_triangular_part == 'U') {
            std::swap(user_model.hessian_row_indices, user_model.hessian_column_indices);
            user_model.hessian_triangular_part = 'L';
        }
    }
    // hessian operator
    user_model.lagrangian_hessian_operator = lagrangian_hessian_operator;
    // initial iterates
    user_model.initial_primal_iterate = initial_primal_iterate;
    user_model.initial_dual_iterate = initial_dual_iterate;

    // create Uno model
    const UnoModel uno_model(user_model);
    
    // create Uno solver
    Uno uno_solver;

    // create Uno options
    Options uno_options;
    // default options
    DefaultOptions::load(uno_options);
    // add default preset
    const Options preset_options = Presets::get_preset_options(std::nullopt);
    uno_options.overwrite_with(preset_options);
    // overwrite with user options
    const Options user_options = mxStruct_to_options(options);
    uno_options.overwrite_with(user_options);   

    // solve
    Logger::set_logger(uno_options.get_string("logger"));
    MxStruct result;
    try {
        Result uno_result = uno_solver.solve(uno_model, uno_options, user_callbacks);
        result = result_to_mxStruct(uno_result);
    }
    catch (const std::exception& e) {
        errmsg = ErrorString::format_error(ErrorType::UNO, e.what());
        mexErrMsgIdAndTxt("uno:error", errmsg.c_str());
    }

    // output
    plhs[0] = mxStruct_to_mxArray(result);

    // flush the logger
    Logger::flush();
}
