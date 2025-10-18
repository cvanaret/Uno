// Copyright (c) 2025 Stefano Lovato and Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include "../cpp_classes/ErrorString.hpp"
#include "unomex_utils.hpp"
#include "unomex_conversion.hpp"
#include "unomex_validation.hpp"

namespace uno {

    // validate Matlab function handle with no calls
    bool validate_matlab_handle_no_call(handle_t handle, int32_t expected_nin, int32_t expected_nout,
        const std::string field_name, std::string& errmsg) {
        if (!isvalid(handle)) {
            errmsg = ErrorString::format_error(ErrorType::MISSING_FIELD, field_name.c_str());
            return false;
        }
        if (!isa<handle_t>(handle)) {
            errmsg = ErrorString::format_error(ErrorType::INVALID_HANDLE);
            return false;
        }
        // nargin
        const int32_t nin = nargin(handle);
        if (nin > expected_nin) {
            errmsg = ErrorString::format_error(ErrorType::NARGIN_TOOMANY);
            return false;
        } 
        else if (nin < expected_nin) {
            errmsg = ErrorString::format_error(ErrorType::NARGIN_NOTENOUGH, expected_nin-nin);
            return false;
        }
        // nargout
        const int32_t nout = nargout(handle);
        if (nout > expected_nout) {
            errmsg = ErrorString::format_error(ErrorType::NARGOUT_TOOMANY);
            return false;
        } 
        else if (nout < expected_nout) {
            errmsg = ErrorString::format_error(ErrorType::NARGOUT_NOTENOUGH, expected_nout-nout);
            return false;
        }
        return true;
    }

    // validate Matlab function handle
    bool validate_matlab_handle(handle_t handle,
        const std::vector<int32_t> input_dims, const std::vector<int32_t> output_dims, 
        const std::string field_name, std::string& errmsg) {
        const int32_t nin = static_cast<int32_t>(input_dims.size());
        const int32_t nout = static_cast<int32_t>(output_dims.size());
        if (!validate_matlab_handle_no_call(handle, nin, nout, field_name, errmsg)) {
            return false;
        }
        // test the function call
        std::vector<mxArray*> inputs(nin);
        std::vector<mxArray*> outputs(nout);
        for (size_t input_index: Range(nin)) {
            const int32_t nrows = input_dims[input_index];
            const Vector<double> input(nrows, 0.);
            inputs[input_index] = vector_to_mxArray(input);
        }
        try {
            call_matlab_function(handle, inputs, outputs);
            destroy_mxArray_vector(inputs);
        } catch (const MatlabFunctionError& err) {
            errmsg = err.what(); 
            return false;
        }
        // validate the outputs
        for (size_t output_index: Range(nout)) {
            const mxArray* output = outputs[output_index];
            const int32_t nrows = output_dims[output_index];
            if (!isa<double>(output)) {
                errmsg = ErrorString::format_error(ErrorType::OUTPUT_TYPE, output_index+1);
                destroy_mxArray_vector(outputs);
                return false;
            }
            if (!isvector(output, nrows)) {
                errmsg = ErrorString::format_error(ErrorType::OUTPUT_SIZE, output_index+1);
                destroy_mxArray_vector(outputs);
                return false;
            }
        }
        destroy_mxArray_vector(outputs);
        return true;
    }

    bool validate_struct_input(const mxArray* arr, const int32_t position, std::string& errmsg) {
        if (!isvalid(arr) || !isa<MxStruct>(arr)) {
            errmsg = ErrorString::format_error(ErrorType::INPUT_TYPE_STRUCT, position);
            return false;
        }
        if (!isscalar(arr)) {
            errmsg = ErrorString::format_error(ErrorType::INPUT_SCALAR, position);
            return false;
        }
        return true;
    }

    bool validate_char(const mxArray* arr, const std::vector<char> chars, const std::string field_name, std::string& errmsg) {
        if (!isvalid(arr)) {
            errmsg = ErrorString::format_error(ErrorType::MISSING_FIELD, field_name.c_str());
            return false;
        }
        if (!isa<char>(arr) || !isscalar(arr)) {
            errmsg = ErrorString::format_error(ErrorType::FIELD_CHAR, field_name.c_str());
            return false;
        }
        if (chars.size()>0) {
            char c = mxArray_to_scalar<char>(arr);
            if (std::find(chars.begin(), chars.end(), c) == chars.end()) {
                errmsg = "Field '" + field_name + "' must be one of { ";
                for (size_t i: Range(chars.size())) {
                    if (i > 0) {
                        errmsg += ", ";
                    }
                    errmsg += "'" + std::string({chars[i]}) + "'";
                }
                errmsg += " }.";
                return false;
            }
        }
        return true;
    }

    bool validate_positive_integer(const mxArray* arr, const std::string field_name, std::string& errmsg) {
        if (!isvalid(arr)) {
            errmsg = ErrorString::format_error(ErrorType::MISSING_FIELD, field_name.c_str());
            return false;
        }
        if (!isa<double>(arr) || !isscalar(arr) || !ispositive(arr)  || !isinteger(arr) ) {
            errmsg = ErrorString::format_error(ErrorType::FIELD_POSITIVE_INTEGER, field_name.c_str());
            return false;
        }
        return true;
    }

    bool validate_unitary(const mxArray* arr, const std::string field_name, std::string& errmsg) {
        if (!isvalid(arr)) {
            errmsg = ErrorString::format_error(ErrorType::MISSING_FIELD, field_name.c_str());
            return false;
        }
        if (!isa<double>(arr) || !isscalar(arr) || !isinteger(arr) || !isunitary(arr) ) {
            errmsg = ErrorString::format_error(ErrorType::FIELD_UNITARY, field_name.c_str());
            return false;
        }
        return true;
    }

    bool validate_double_vector(const mxArray* arr, const int32_t len, const std::string field_name, std::string& errmsg) {
        if (!isvalid(arr)) {
            errmsg = ErrorString::format_error(ErrorType::MISSING_FIELD, field_name.c_str());
            return false;
        }
        if (!isa<double>(arr) || !isvector(arr, len) ) {
            errmsg = ErrorString::format_error(ErrorType::FIELD_VECTOR, field_name.c_str(), len);
            return false;
        }
        return true;
    }

    bool validate_integer_vector(const mxArray* arr, const int32_t len, const std::string field_name, std::string& errmsg) {
        if (!isvalid(arr)) {
            errmsg = ErrorString::format_error(ErrorType::MISSING_FIELD, field_name.c_str());
            return false;
        }
        if (!isa<double>(arr) || !isvector(arr, len) ) {
            errmsg = ErrorString::format_error(ErrorType::FIELD_VECTOR_INT, field_name.c_str(), len);
            return false;
        }
        const Vector vector = mxArray_to_vector<double>(arr);
        double intpart;
        for (auto value : vector) {
            if (std::modf(value, &intpart) != 0.0) {
                errmsg = ErrorString::format_error(ErrorType::FIELD_VECTOR_INT, field_name.c_str(), len);
                return false;
            }
        }
        return true;
    }
    
} // namespace