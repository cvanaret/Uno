// Copyright (c) 2025 Stefano Lovato and Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <string>
#include <cmath>
#include "unomex_utils.hpp"
#include "unomex_function.hpp"
#include "../cpp_classes/MxStruct.hpp"

namespace uno {

    MxArrayVectorGuard::MxArrayVectorGuard(const std::vector<mxArray*>& vector)
    : vector(vector) {}

    MxArrayVectorGuard::~MxArrayVectorGuard() {
        for (auto* arr : vector) {
            mxDestroyArray(arr);
        }
    }

    template<> mxClassID get_mxClassID<double>() { return mxDOUBLE_CLASS; }
    template<> mxClassID get_mxClassID<float>() { return mxSINGLE_CLASS; }
    template<> mxClassID get_mxClassID<int8_t>() { return mxINT8_CLASS; }
    template<> mxClassID get_mxClassID<uint8_t>() { return mxUINT8_CLASS; }
    template<> mxClassID get_mxClassID<int16_t>() { return mxINT16_CLASS; }
    template<> mxClassID get_mxClassID<uint16_t>() { return mxUINT16_CLASS; }
    template<> mxClassID get_mxClassID<int32_t>() { return mxINT32_CLASS; }
    template<> mxClassID get_mxClassID<uint32_t>() { return mxUINT32_CLASS; }
    template<> mxClassID get_mxClassID<int64_t>() { return mxINT64_CLASS; }
    template<> mxClassID get_mxClassID<uint64_t>() { return mxUINT64_CLASS; }
    template<> mxClassID get_mxClassID<bool>() { return mxLOGICAL_CLASS; }
    template<> mxClassID get_mxClassID<char>() { return mxCHAR_CLASS; }
    template<> mxClassID get_mxClassID<void>() { return mxVOID_CLASS; }
    template<> mxClassID get_mxClassID<handle_t>() { return mxFUNCTION_CLASS; }
    template<> mxClassID get_mxClassID<MxStruct>() { return mxSTRUCT_CLASS; }
    template<> mxClassID get_mxClassID<std::string>() { return static_cast<mxClassID>(mxSTRING_CLASS); }

    bool isvalid(const mxArray* arr) {
        return (arr != nullptr);
    }
    
    bool isempty(const mxArray* arr) {
        return mxIsEmpty(arr);
    }

    bool has_size(const mxArray* arr, const size_t nrows, const size_t ncolumns) {
        if (mxGetNumberOfDimensions(arr)>2) {
            return false;;
        }
        return (static_cast<size_t>(mxGetM(arr))==nrows) && (static_cast<size_t>(mxGetN(arr))==ncolumns);
    }

    bool ispositive(const mxArray* arr) {
        const double value = mxGetScalar(arr);
        return value>=0;
    }

    bool isinteger(const mxArray* arr) {
        const double value = mxGetScalar(arr);
        if (!std::isfinite(value)) {
            return false;
        }
        double intpart;
        return std::modf(value, &intpart) == 0.0;
    }

    bool isunitary(const mxArray* arr) {
        const double value = mxGetScalar(arr);
        return (value==1.0) || (value==-1.0);
    }

    bool isscalar(const mxArray* arr) {
        return mxGetNumberOfElements(arr)==1;
    }

    bool isvector(const mxArray* arr, const size_t len) {
        return has_size(arr, len, 1) || has_size(arr, 1, len);
    }

}; // namespace