// Copyright (c) 2025 Stefano Lovato and Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNOMEX_FUNCTION_H
#define UNOMEX_FUNCTION_H

#include <stdint.h>
#include <string>
#include <vector>
#include "optimization/EvaluationErrors.hpp"
#include "mex.h"

namespace uno {

    // matlab function handle type
    typedef mxArray* handle_t;

    struct MatlabFunctionError : EvaluationError {
        std::string msg;
        explicit MatlabFunctionError(const std::string& str);
        const char* what() const noexcept override;
    };

    void call_matlab_function(handle_t handle, const std::vector<mxArray*>& inputs, std::vector<mxArray*>& outputs);
    
    int32_t nargin(handle_t handle);
    int32_t nargout(handle_t handle);

}; // namespace


#endif // UNOMEX_FUNCTION_H