// Copyright (c) 2025 Stefano Lovato and Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNOMEX_VALIDATION_H
#define UNOMEX_VALIDATION_H

#include <string>
#include <vector>
#include "unomex_function.hpp"
#include "mex.h"

namespace uno {

    bool validate_matlab_handle(handle_t handle, const std::string field_name, std::string& errmsg);

    bool validate_struct_input(const mxArray* arr, const size_t position, std::string& errmsg);

    bool validate_char_field(const mxArray* arr, const std::vector<char> chars, const std::string field_name, std::string& errmsg);
    bool validate_positive_integer_field(const mxArray* arr, const std::string field_name, std::string& errmsg);
    bool validate_unitary_field(const mxArray* arr, const std::string field_name, std::string& errmsg);
    
    bool validate_double_vector_field(const mxArray* arr, const size_t len, const std::string field_name, std::string& errmsg);
    bool validate_integer_vector_field(const mxArray* arr, const size_t len, const std::string field_name, std::string& errmsg);

    bool validate_double_scalar_output(const mxArray* arr,  std::string& errmsg);
    bool validate_bool_scalar_output(const mxArray* arr, std::string& errmsg);
    bool validate_double_vector_output(const mxArray* arr, size_t n, std::string& errmsg);
    bool validate_gradient(const mxArray* arr, size_t n, std::string& errmsg);
    bool validate_constraint_jacobian(const mxArray* arr, size_t n, std::string& errmsg);
    
} // namespace



#endif // UNOMEX_VALIDATION_H