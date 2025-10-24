// Copyright (c) 2025 Stefano Lovato and Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <stdint.h>
#include <algorithm>
#include "Uno.hpp"
#include "options/DefaultOptions.hpp"
#include "options/Presets.hpp"
#include "tools/Logger.hpp"
#include "cpp_classes/MatlabModel.hpp"
#include "cpp_classes/MatlabUserCallbacks.hpp"
#include "cpp_classes/MatlabStreamCallback.hpp"
#include "cpp_classes/MxStruct.hpp"
#include "cpp_classes/ErrorString.hpp"
#include "unomex/unomex_function.hpp"
#include "unomex/unomex_conversion.hpp"
#include "unomex/unomex_validation.hpp"
#include "unomex/unomex_utils.hpp"
#include "mex.h"

using namespace uno;

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
    if (!validate_string_field(problem_type, {UNO_PROBLEM_LINEAR, UNO_PROBLEM_QUADRATIC, UNO_PROBLEM_NONLINEAR}, "problem_type", errmsg)) {
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
    if (!validate_matlab_handle_field(objective_function, "objective_function", errmsg)) {
        mexErrMsgIdAndTxt("uno:error", "Invalid model objective function. %s", errmsg.c_str());
    }
    // gradient = objective_gradient(x)
    if (!validate_matlab_handle_field(objective_gradient, "objective_gradient", errmsg)) {
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
    if (!validate_matlab_handle_field(constraint_function, "constraint_function", errmsg)) {
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
    if (!validate_matlab_handle_field(constraint_jacobian, "constraint_jacobian", errmsg)) {
        mexErrMsgIdAndTxt("uno:error", "Invalid model constraint Jacobian. %s", errmsg.c_str());
    }

    // jacobian operators (optional)
    mxArray* jacobian_operator = model["jacobian_operator"];
    mxArray* jacobian_transposed_operator = model["jacobian_transposed_operator"];
    if (jacobian_operator != nullptr) {
        // result = jacobian_operator(x,v)
        if (!validate_matlab_handle_field(jacobian_operator, "jacobian_operator", errmsg)) {
            mexErrMsgIdAndTxt("uno:error", "Invalid model Jacobian operator. %s", errmsg.c_str());
        }
    }
    if (jacobian_transposed_operator != nullptr) {
        // jacobian_transposed_operator(x,v)
        if (!validate_matlab_handle_field(jacobian_transposed_operator, "jacobian_transposed_operator", errmsg)) {
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
        if (!validate_matlab_handle_field(lagrangian_hessian, "lagrangian_hessian", errmsg)) {
            mexErrMsgIdAndTxt("uno:error", "Invalid model Lagrangian Hessian. %s", errmsg.c_str());
        }
    }

    // hessian operators (optional)
    mxArray* lagrangian_hessian_operator = model["lagrangian_hessian_operator"];
    if (lagrangian_hessian_operator != nullptr) {
        // lagrangian_hessian_operator(x,rho,y,v)
        if (!validate_matlab_handle_field(lagrangian_hessian_operator, "lagrangian_hessian_operator", errmsg)) {
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
        if (!validate_matlab_handle_field(logger_stream_callback, "logger_stream_callback", errmsg)) {
            mexErrMsgIdAndTxt("uno:error", "Invalid logger stream callback. %s", errmsg.c_str());
        }
    }
    if (notify_acceptable_iterate_callback != nullptr) {
        // notify_acceptable_iterate_callback(x, yl, yb, y, rho, feas, stat, compl)
        if (!validate_matlab_handle_field(notify_acceptable_iterate_callback, "notify_acceptable_iterate_callback", errmsg)) {
            mexErrMsgIdAndTxt("uno:error", "Invalid notify acceptable iterate callback. %s", errmsg.c_str());
        }
    }
    if (user_termination_callback != nullptr) {
        // terminate = user_termination_callback(x, yl, yb, y, rho, feas, stat, compl)
        if (!validate_matlab_handle_field(user_termination_callback, "user_termination_callback", errmsg)) {
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
    MatlabUserModel user_model(mxArray_to_string(problem_type).c_str(), 
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
    const MatlabModel uno_model(user_model);
    
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
