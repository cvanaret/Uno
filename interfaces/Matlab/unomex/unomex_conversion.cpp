// Copyright (c) 2025 Stefano Lovato and Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include "tools/Logger.hpp"
#include "../cpp_classes/ErrorString.hpp"
#include "unomex_utils.hpp"
#include "unomex_conversion.hpp"

namespace uno {

    std::string mxArray_to_string(const mxArray* arr) {
        std::string str;
        if (!isvalid(arr)) {
            return str;
        }
        else if (mxIsChar(arr)) {
            char* cstr = mxArrayToString(arr);
            str = std::string(cstr);
            mxFree(cstr);
            return str;
        }
        else if (mxIsClass(arr,"string")) {
            // convert matlab string to char and call mxArray_to_string
            // TODO handle non-scalar matlab string
            mxArray* arr_str = mxDuplicateArray(arr);
            mxArray* arr_char;
            mexCallMATLAB(1, &arr_char, 1, &arr_str, "char");
            str = mxArray_to_string(arr_char);
            mxDestroyArray(arr_str);
            mxDestroyArray(arr_char);
            return str;
        }
        else {
            return str;
        }
    }

    mxArray* string_to_mxArray(const std::string str) {
        return mxCreateString(str.c_str());
    }

    MxStruct mxArray_to_mxStruct(const mxArray* s) {
        if (!s) return MxStruct();
        MxStruct result;
        int nFields = mxGetNumberOfFields(s);
        for (int i = 0; i < nFields; ++i) {
            const char* fieldName = mxGetFieldNameByNumber(s, i);
            mxArray* fieldValue = mxGetFieldByNumber(s, 0, i);
            result.insert(fieldName, fieldValue);
        }
        return result;
    }

    mxArray* mxStruct_to_mxArray(const MxStruct& mxStruct) {
        Vector<const char*> fieldNames;
        for (const auto& [key, value]: mxStruct) {
            fieldNames.push_back(key.c_str());
        }
        mxArray* s = mxCreateStructMatrix(1, 1, static_cast<int>(fieldNames.size()), fieldNames.data());
        for (const auto& [key, value]: mxStruct) {
            mxSetField(s, 0, key.c_str(), value);
        }
        return s;
    }

    MxStruct options_to_mxStruct(const Options& uno_options) {
        MxStruct mxStruct;
        // integer options (matlab double)
        std::map<std::string, int32_t> integer_options = uno_options.get_integer_options();
        for (const auto& [option_name, option_value]: integer_options) {
            mxStruct.insert(option_name, scalar_to_mxArray<double>(static_cast<double>(option_value)));
        }
        // double options (matlab double)
        std::map<std::string, double> double_options = uno_options.get_double_options();
        for (const auto& [option_name, option_value]: double_options) {
            mxStruct.insert(option_name, scalar_to_mxArray<double>(option_value));
        }
        // bool options (matlab logical)
        std::map<std::string, bool> bool_options = uno_options.get_bool_options();
        for (const auto& [option_name, option_value]: bool_options) {
            mxStruct.insert(option_name, scalar_to_mxArray<bool>(option_value));
        }
        // string options (matlab char)
        std::map<std::string, std::string> string_options = uno_options.get_string_options();
        for (const auto& [option_name, option_value]: string_options) {
            mxStruct.insert(option_name, string_to_mxArray(option_value));
        }        
        return mxStruct;
    }

    Options mxStruct_to_options(const MxStruct& options) {
        Options uno_options;
        for (const auto& [field_name, field_value]: options) {
            try {
                // set option with check between type and mxArray class
                const OptionType option_type = uno_options.get_option_type(field_name);
                // OptionType::DOUBLE accepts double
                if (option_type == OptionType::DOUBLE) {
                    if (isa<double>(field_value)) {
                        uno_options.set_double(field_name, mxArray_to_scalar<double>(field_value));
                    }
                    else {
                        INFO << ErrorString::format_error(ErrorType::OPTION_TYPE, field_name.c_str()) << std::endl;
                    }
                }
                // OptionType::INTEGER accepts double
                else if (option_type == OptionType::INTEGER) {
                    if (isa<double>(field_value)) {
                        uno_options.set_integer(field_name, static_cast<int32_t>(mxArray_to_scalar<double>(field_value))); // cast to int
                    }
                    else {
                        INFO << ErrorString::format_error(ErrorType::OPTION_TYPE, field_name.c_str()) << std::endl;
                    }
                }
                // OptionType::BOOL accepts logical
                else if (option_type == OptionType::BOOL) {
                    if (isa<bool>(field_value)) { 
                        uno_options.set_bool(field_name, mxArray_to_scalar<bool>(field_value)); 
                    }
                    else {
                        INFO << ErrorString::format_error(ErrorType::OPTION_TYPE, field_name.c_str()) << std::endl;
                    }
                }
                // OptionType::STRING accepts char
                else if (option_type == OptionType::STRING) {
                    if (isa<char>(field_value) || isa<std::string>(field_value)) { 
                        uno_options.set_string(field_name, mxArray_to_string(field_value)); 
                    }
                    else {
                        INFO << ErrorString::format_error(ErrorType::OPTION_TYPE, field_name.c_str()) << std::endl;
                    }
                }
            } catch (const std::out_of_range&) {
                // set the option with type depending on mxArray class (only double, bool, char)
                if (isa<double>(field_value)) {
                    uno_options.set_double(field_name, mxArray_to_scalar<double>(field_value));
                }
                else if (isa<bool>(field_value)) {
                    uno_options.set_bool(field_name, mxArray_to_scalar<bool>(field_value));
                }
                else if (isa<char>(field_value) || isa<std::string>(field_value)) {
                    uno_options.set_string(field_name, mxArray_to_string(field_value));
                }
                else {
                    INFO << ErrorString::format_error(ErrorType::OPTION_TYPE, field_name.c_str()) << std::endl;
                }
            }
        }
        return uno_options;
    }

    MxStruct result_to_mxStruct(const Result& uno_result) {
        MxStruct result;
        result.insert("optimization_status", scalar_to_mxArray(static_cast<double>(uno_result.optimization_status)));
        result.insert("solution_status", scalar_to_mxArray(static_cast<double>(uno_result.solution_status)));
        result.insert("solution_objective", scalar_to_mxArray(uno_result.solution_objective));
        result.insert("solution_primal_feasibility", scalar_to_mxArray(uno_result.solution_primal_feasibility));
        result.insert("solution_stationarity", scalar_to_mxArray(uno_result.solution_stationarity));
        result.insert("solution_complementarity", scalar_to_mxArray(uno_result.solution_complementarity));
        result.insert("primal_solution", vector_to_mxArray(uno_result.primal_solution));
        result.insert("constraint_dual_solution", vector_to_mxArray(uno_result.constraint_dual_solution));
        result.insert("lower_bound_dual_solution", vector_to_mxArray(uno_result.lower_bound_dual_solution));
        result.insert("upper_bound_dual_solution", vector_to_mxArray(uno_result.upper_bound_dual_solution));
        result.insert("number_iterations", scalar_to_mxArray(static_cast<double>(uno_result.number_iterations)));
        result.insert("cpu_time", scalar_to_mxArray(uno_result.cpu_time));
        result.insert("number_objective_evaluations", scalar_to_mxArray(static_cast<double>(uno_result.number_objective_evaluations)));
        result.insert("number_constraint_evaluations", scalar_to_mxArray(static_cast<double>(uno_result.number_constraint_evaluations)));
        result.insert("number_objective_gradient_evaluations", scalar_to_mxArray(static_cast<double>(uno_result.number_objective_gradient_evaluations)));
        result.insert("number_jacobian_evaluations", scalar_to_mxArray(static_cast<double>(uno_result.number_jacobian_evaluations)));
        result.insert("number_hessian_evaluations", scalar_to_mxArray(static_cast<double>(uno_result.number_hessian_evaluations)));
        result.insert("number_subproblems_solved", scalar_to_mxArray(static_cast<double>(uno_result.number_subproblems_solved)));
        return result;
    }

} // namespace