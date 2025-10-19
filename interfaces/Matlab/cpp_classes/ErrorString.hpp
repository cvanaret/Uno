// Copyright (c) 2025 Stefano Lovato and Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_ERRORSTRING_H
#define UNO_ERRORSTRING_H

#include <string>
#include <map>

namespace uno {

    enum class ErrorType {
        NARGIN_NOTENOUGH,
        NARGIN_TOOMANY,
        INPUT_INVALID,
        INPUT_STRING,
        INPUT_STRUCT,
        INPUT_SCALAR,
        OUTPUT_SCALAR,
        OUTPUT_VECTOR,
        OUTPUT_BOOL_SCALAR,
        INVALID_HANDLE,
        EVALUATION,
        OPTION_TYPE,
        MISSING_FIELD,
        FIELD_CHAR,
        FIELD_POSITIVE_INTEGER,
        FIELD_UNITARY,
        FIELD_VECTOR,
        FIELD_VECTOR_INT,
        UNO
    };

    // map the error type to error strings
    class ErrorString {
    public:
        static std::string format_error(ErrorType error_type, ...);
    private:
        static const std::map<ErrorType, std::string> error_strings;
    };

}; // namespace

#endif // UNO_ERRORSTRING_H