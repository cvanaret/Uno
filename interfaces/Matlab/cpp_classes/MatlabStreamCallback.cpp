// Copyright (c) 2025 Stefano Lovato and Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <vector>
#include "MatlabStreamCallback.hpp"
#include "../unomex/unomex_function.hpp"
#include "../unomex/unomex_conversion.hpp"
#include "../unomex/unomex_validation.hpp"
#include "../unomex/unomex_utils.hpp"
#include "mex.h"

namespace uno {
    
    MatlabStreamCallback::MatlabStreamCallback(handle_t logger_stream_callback) :
    UserStreamCallback(), logger_stream_callback(logger_stream_callback) { }

    int32_t MatlabStreamCallback::operator()(const char* buf, int32_t len) const {
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
    
} // namespace
