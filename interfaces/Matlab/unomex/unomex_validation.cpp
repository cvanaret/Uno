// Copyright (c) 2025 Stefano Lovato and Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include "../cpp_classes/ErrorString.hpp"
#include "unomex_utils.hpp"
#include "unomex_conversion.hpp"
#include "unomex_validation.hpp"

namespace uno {

    // validate Matlab function handle
    bool validate_matlab_handle_field(handle_t handle, const std::string field_name, std::string& errmsg) {
        // Matlab throws error when wrong nargin/nargour etc... just check the validity of the handle
        if (!isvalid(handle)) {
            errmsg = ErrorString::format_error(ErrorType::MISSING_FIELD, field_name.c_str());
            return false;
        }
        if (!isa<handle_t>(handle)) {
            errmsg = ErrorString::format_error(ErrorType::INVALID_HANDLE);
            return false;
        }
        return true;
    }

    bool validate_struct_input(const mxArray* arr, const size_t position, std::string& errmsg) {
        if (!isvalid(arr) || !isa<MxStruct>(arr)) {
            errmsg = ErrorString::format_error(ErrorType::INPUT_STRUCT, position);
            return false;
        }
        if (!isscalar(arr)) {
            errmsg = ErrorString::format_error(ErrorType::INPUT_SCALAR, position);
            return false;
        }
        return true;
    }

    bool validate_char_field(const mxArray* arr, const std::vector<char> chars, const std::string field_name, std::string& errmsg) {
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

    bool validate_string_field(const mxArray* arr, const std::vector<std::string> strings, const std::string field_name, std::string& errmsg) {
        if (!isvalid(arr)) {
            errmsg = ErrorString::format_error(ErrorType::MISSING_FIELD, field_name.c_str());
            return false;
        }
        if (!isa<char>(arr) && !isa<std::string>(arr)) {
            errmsg = ErrorString::format_error(ErrorType::FIELD_STRING, field_name.c_str());
            return false;
        }
        if (strings.size()>0) {
            std::string str = mxArray_to_string(arr);
            if (std::find(strings.begin(), strings.end(), str) == strings.end()) {
                errmsg = "Field '" + field_name + "' must be one of { ";
                for (size_t i: Range(strings.size())) {
                    if (i > 0) {
                        errmsg += ", ";
                    }
                    errmsg += "'" + strings[i] + "'";
                }
                errmsg += " }.";
                return false;
            }
        }
    }

    bool validate_positive_integer_field(const mxArray* arr, const std::string field_name, std::string& errmsg) {
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

    bool validate_unitary_field(const mxArray* arr, const std::string field_name, std::string& errmsg) {
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

    bool validate_double_vector_field(const mxArray* arr, const size_t len, const std::string field_name, std::string& errmsg) {
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

    bool validate_integer_vector_field(const mxArray* arr, const size_t len, const std::string field_name, std::string& errmsg) {
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

    bool validate_double_scalar_output(const mxArray* arr, std::string& errmsg) {
        if (!isa<double>(arr) || !isscalar(arr)) {
            errmsg = ErrorString::format_error(ErrorType::OUTPUT_SCALAR, 1);
            return false;
        }
        return true;
    }

    bool validate_double_vector_output(const mxArray* arr, size_t n, std::string& errmsg) {
        if (!isa<double>(arr) || !isvector(arr, n)) {
            errmsg = ErrorString::format_error(ErrorType::OUTPUT_VECTOR, 1, n);
            return false;
        }
        return true;
    }

    bool validate_bool_scalar_output(const mxArray* arr, std::string& errmsg) {
        if (!isa<bool>(arr) || !isscalar(arr)) {
            errmsg = ErrorString::format_error(ErrorType::OUTPUT_BOOL_SCALAR, 1);
            return false;
        }
        return true;
    }

} // namespace