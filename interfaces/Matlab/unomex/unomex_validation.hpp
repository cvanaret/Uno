// Copyright (c) 2025 Stefano Lovato and Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNOMEX_VALIDATION_H
#define UNOMEX_VALIDATION_H

#include <string>
#include <vector>
#include "unomex_function.hpp"
#include "mex.h"

namespace uno {

    bool validate_matlab_handle_no_call(handle_t handle, int32_t expected_nin, int32_t expected_nout,
        const std::string field_name, std::string& errmsg);
    bool validate_matlab_handle(handle_t handle,
        const std::vector<int32_t> input_dims, const std::vector<int32_t> output_dims, 
        const std::string field_name, std::string& errmsg);

    bool validate_struct_input(const mxArray* arr, const int32_t position, std::string& errmsg);

    bool validate_char(const mxArray* arr, const std::vector<char> chars, const std::string field_name, std::string& errmsg);
    bool validate_positive_integer(const mxArray* arr, const std::string field_name, std::string& errmsg);
    bool validate_unitary(const mxArray* arr, const std::string field_name, std::string& errmsg);
    
    bool validate_double_vector(const mxArray* arr, const int32_t len, const std::string field_name, std::string& errmsg);
    bool validate_integer_vector(const mxArray* arr, const int32_t len, const std::string field_name, std::string& errmsg);
    
} // namespace



#endif // UNOMEX_VALIDATION_H