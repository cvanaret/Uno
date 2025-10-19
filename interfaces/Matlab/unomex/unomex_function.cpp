// Copyright (c) 2025 Stefano Lovato and Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include "unomex_function.hpp"
#include "unomex_conversion.hpp"
#include "../cpp_classes\ErrorString.hpp"

namespace uno {

    MatlabFunctionError::MatlabFunctionError(const std::string& str) 
        : msg(ErrorString::format_error(ErrorType::EVALUATION) + "\n" + str) {}

    const char* MatlabFunctionError::what() const noexcept {
        return this->msg.c_str();
    }

    void call_matlab_function(handle_t handle, const std::vector<mxArray*>& inputs, std::vector<mxArray*>& outputs) {
        std::vector<mxArray*> inputs_all = inputs;
        inputs_all.insert(inputs_all.begin(), handle);
        const int nrhs = static_cast<int>(inputs_all.size());
        const int nlhs = static_cast<int>(outputs.size());
        mxArray** prhs = inputs_all.data();
        mxArray** plhs = outputs.data();
        mxArray* err = mexCallMATLABWithTrap(nlhs, plhs, nrhs, prhs, "feval");
        if (err) {
            mxArray* errmsg;
            mexCallMATLAB(1, &errmsg, 1, &err, "getReport");
            const std::string strerr = mxArrayToString(errmsg);
            mxDestroyArray(err);
            mxDestroyArray(errmsg);
            throw MatlabFunctionError(strerr);
        }
    }

    int32_t nargin(handle_t handle) {
        mxArray* out;
        mexCallMATLAB(1, &out, 1, &handle, "nargin");
        const double n = mxArray_to_scalar<double>(out);
        mxDestroyArray(out);
        return static_cast<int32_t>(n);
    }

    int32_t nargout(handle_t handle) {
        mxArray* out;
        mexCallMATLAB(1, &out, 1, &handle, "nargout");
        const double n = mxArray_to_scalar<double>(out);
        mxDestroyArray(out);
        return static_cast<int32_t>(n);
    }
}