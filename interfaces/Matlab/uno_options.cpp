// Copyright (c) 2025 Stefano Lovato and Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include "Uno.hpp"
#include "options/DefaultOptions.hpp"
#include "options/Presets.hpp"
#include "cpp_classes/ErrorString.hpp"
#include "unomex/unomex_conversion.hpp"
#include "unomex/unomex_utils.hpp"
#include "mex.h"

using namespace uno;

// gateway function: options = uno_options([preset]);
void mexFunction( int /* nlhs */, mxArray* plhs[], int nrhs, const mxArray* prhs[]) {
    // validate arguments
    if (nrhs > 1) {
        const std::string errmsg = ErrorString::format_error(ErrorType::NARGIN_TOOMANY);
        mexErrMsgIdAndTxt("uno:error", errmsg.c_str());
    }
    // preset (optional)
    std::string preset;
    if (nrhs == 1 && !isempty(prhs[0])) {
        if (!isa<char>(prhs[0]) && !isa<std::string>(prhs[0])) {
            const std::string errmsg = ErrorString::format_error(ErrorType::INPUT_STRING, 1);
            mexErrMsgIdAndTxt("uno:error", errmsg.c_str());
        }
        preset = mxArray_to_string(prhs[0]);
    }

    // create Uno options
    Options uno_options;

    // default options
    DefaultOptions::load(uno_options);

    // add default preset
    const Options preset_options = Presets::get_preset_options(std::nullopt);
    uno_options.overwrite_with(preset_options);

    // set preset
    if (!preset.empty()) {
        try {
            Presets::set(uno_options, preset);
        }
        catch (const std::runtime_error& err) {
            const std::string errmsg = ErrorString::format_error(ErrorType::INPUT_INVALID, 1);
            mexErrMsgIdAndTxt("uno:error", "%s %s.", errmsg.c_str(), err.what());
        }
    }

    // options
    MxStruct options = options_to_mxStruct(uno_options);

    // output
    plhs[0] = mxStruct_to_mxArray(options);
}