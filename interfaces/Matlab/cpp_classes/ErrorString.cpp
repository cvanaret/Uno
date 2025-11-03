// Copyright (c) 2025 Stefano Lovato and Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include "ErrorString.hpp"
#include <cstdarg>

namespace uno {

    const std::map<ErrorType, std::string> ErrorString::error_strings = {
        {ErrorType::NARGIN_NOTENOUGH, "Invalid argument list. Function requires %d more input(s)."},
        {ErrorType::NARGIN_TOOMANY, "Too many input arguments."},
        {ErrorType::INPUT_INVALID, "Invalid input argument at position %d."},
        {ErrorType::INPUT_STRING, "Invalid argument at position %d. Value must be char or string."},
        {ErrorType::INPUT_STRUCT, "Invalid argument at position %d. Value must be struct."},  
        {ErrorType::INPUT_SCALAR, "Invalid input argument at position %d. Value must be scalar."},
        {ErrorType::OUTPUT_SCALAR, "Invalid output argument at position %d. Value must be a numeric scalar."},
        {ErrorType::OUTPUT_VECTOR, "Invalid output argument at position %d. Value must be a numeric vector of %d element(s)."},
        {ErrorType::OUTPUT_BOOL_SCALAR, "Invalid output argument at position %d. Value must be a logical scalar."},
        {ErrorType::INVALID_HANDLE, "Value must be a function handle."},
        {ErrorType::EVALUATION, "Error in function evaluation."},
        {ErrorType::OPTION_TYPE, "Invalid type of option '%s'."},
        {ErrorType::MISSING_FIELD, "Missing field '%s'."},
        {ErrorType::FIELD_CHAR, "Type of field '%s' must be char."},
        {ErrorType::FIELD_STRING, "Type of field '%s' must be char array or string."},
        {ErrorType::FIELD_POSITIVE_INTEGER, "Field '%s' must be a positive integer value."},
        {ErrorType::FIELD_UNITARY, "Field '%s' must be either +1 or -1."},
        {ErrorType::FIELD_VECTOR, "Field '%s' must be a vector of %d element(s)."},
        {ErrorType::FIELD_VECTOR_INT, "Field '%s' must be a vector of %d integer element(s)."},
        {ErrorType::UNO, "Error in UNO. %s"}
    };

    std::string ErrorString::format_error(ErrorType error_type, ...) {
        auto it = error_strings.find(error_type);
        if (it == error_strings.end()) {
            return "Unknown error.";
        }
        // format the string (avoid using std::format which is C++20)
        const std::string& error_fmt = it->second;
        va_list args1;
        va_start(args1, error_type);
        va_list args2;
        va_copy(args2, args1);
        int size = std::vsnprintf(nullptr, 0, error_fmt.c_str(), args2);
        va_end(args2);
        std::string error_msg(size, '\0');
        std::vsnprintf(error_msg.data(), error_msg.size() + 1, error_fmt.c_str(), args1);
        va_end(args1);
        return error_msg;
    }

} // namespace 