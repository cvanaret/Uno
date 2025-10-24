// Copyright (c) 2025 Stefano Lovato and Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <vector>
#include "MatlabUserCallbacks.hpp"
#include "../unomex/unomex_function.hpp"
#include "../unomex/unomex_conversion.hpp"
#include "../unomex/unomex_validation.hpp"
#include "../unomex/unomex_utils.hpp"
#include "mex.h"

#ifdef HAS_MXLIBUT
// These functions are to allow interrupt from MATLAB using CTRL+C (see 
// https://undocumentedmatlab.com/articles/mex-ctrl-c-interrupt) 
extern "C" bool utIsInterruptPending();
extern "C" bool utSetInterruptPending(bool);
#else
// libut is not available, use dummy
constexpr bool utIsInterruptPending() { return false; }
constexpr bool utSetInterruptPending(bool) { return false; }
#endif

namespace uno {

    MatlabUserCallbacks::MatlabUserCallbacks(handle_t notify_acceptable_iterate_callback, handle_t user_termination_callback) : UserCallbacks(),
        notify_acceptable_iterate_callback(notify_acceptable_iterate_callback),
        user_termination_callback(user_termination_callback) { }

    void MatlabUserCallbacks::notify_acceptable_iterate(const Vector<double>& primals, const Multipliers& multipliers, double objective_multiplier, double primal_feasibility, double stationarity, double complementarity) {
        // notify_acceptable_iterate_callback(x, yl, yb, y, rho, feas, stat, compl);
        if (this->notify_acceptable_iterate_callback) {
            std::vector<mxArray*> inputs = {vector_to_mxArray(primals), vector_to_mxArray(multipliers.lower_bounds), vector_to_mxArray(multipliers.upper_bounds), vector_to_mxArray(multipliers.constraints), scalar_to_mxArray(objective_multiplier), scalar_to_mxArray(primal_feasibility), scalar_to_mxArray(stationarity), scalar_to_mxArray(complementarity)};
            std::vector<mxArray*> outputs(0);
            MxArrayVectorGuard input_guard(inputs);
            // MxArrayVectorGuard output_guard(outputs);
            try {
                call_matlab_function(this->notify_acceptable_iterate_callback, inputs, outputs);
            }
            catch (const MatlabFunctionError& err) {
                mexWarnMsgIdAndTxt("uno:warning", "%s", err.what());
            }
            // no output to validate
        }
    }

    bool MatlabUserCallbacks::user_termination(const Vector<double>& primals, const Multipliers& multipliers, double objective_multiplier, double primal_feasibility, double stationarity, double complementarity) {
        // handle matlab CTRL+C event
        if (utIsInterruptPending()) {
            utSetInterruptPending(false);
            return true;
        }
        // terminate = user_termination_callback(x, yl, yb, y, rho, feas, stat, compl);
        if (this->user_termination_callback) {
            std::vector<mxArray*> inputs = {vector_to_mxArray(primals), vector_to_mxArray(multipliers.lower_bounds), vector_to_mxArray(multipliers.upper_bounds), vector_to_mxArray(multipliers.constraints), scalar_to_mxArray(objective_multiplier), scalar_to_mxArray(primal_feasibility), scalar_to_mxArray(stationarity), scalar_to_mxArray(complementarity)};
            std::vector<mxArray*> outputs(1);
            MxArrayVectorGuard input_guard(inputs);
            MxArrayVectorGuard output_guard(outputs);
            try {
                call_matlab_function(this->user_termination_callback, inputs, outputs);
            }  
            catch (const MatlabFunctionError& err) {
                mexWarnMsgIdAndTxt("uno:warning", "%s", err.what());
                return false;
            }
            if (std::string errmsg; !validate_bool_scalar_output(outputs[0], errmsg)) {
                mexWarnMsgIdAndTxt("uno:warning", "Error in user termination callback.\n%s\n", errmsg.c_str());
                return false;
            }
            return mxArray_to_scalar<bool>(outputs[0]);
        } else {
            return false; // never terminate
        }
    }
    
} // namespace