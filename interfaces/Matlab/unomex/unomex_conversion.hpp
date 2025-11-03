// Copyright (c) 2025 Stefano Lovato and Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNOMEX_CONVERSION_H
#define UNOMEX_CONVERSION_H

#include "linear_algebra/Vector.hpp"
#include "options/Options.hpp"
#include "optimization/Result.hpp"
#include "../cpp_classes/MxStruct.hpp"
#include "unomex_utils.hpp"
#include "mex.h"

namespace uno {

    template<typename T>
    void mxArray_to_pointer(const mxArray* arr, T* pointer) {
        const size_t n = mxGetNumberOfElements(arr);
        const T* ptr = static_cast<T*>(mxGetData(arr));
        std::copy(ptr, ptr+n, pointer);
    }

    template<typename T>
    mxArray* pointer_to_mxArray(const T* pointer, const size_t n) {
        const mxClassID classid = get_mxClassID<T>();
        mxArray* arr = mxCreateNumericMatrix(n, 1, classid, mxREAL);
        T* ptr = static_cast<T*>(mxGetData(arr));
        std::copy(pointer, pointer+n, ptr);
        return arr;
    }

    template<typename T>
    T mxArray_to_scalar(const mxArray* arr) {
        return *(static_cast<T*>(mxGetData(arr)));
    }

    template<typename T>
    mxArray* scalar_to_mxArray(const T value) {
        const mxClassID classid = get_mxClassID<T>();
        mxArray* arr = mxCreateNumericMatrix(1, 1, classid, mxREAL);
        T* ptr = static_cast<T*>(mxGetData(arr));
        *ptr = value;
        return arr;
    }

    template<typename T>
    void mxArray_to_vector(const mxArray* arr, Vector<T>& vec) {
        const size_t n = mxGetNumberOfElements(arr); // assume vec has the correct size
        const T* ptr = mxGetPr(arr);
        std::copy(ptr, ptr+n, vec.data());
    }
    template<typename T>
    Vector<T> mxArray_to_vector(const mxArray* arr) {
        const size_t n = mxGetNumberOfElements(arr);
        Vector<T> vec(n);
        mxArray_to_vector(arr, vec);
        return vec;
    }

    template<typename T>
    mxArray* vector_to_mxArray(const Vector<T>& vec) {
        const mwSize n = static_cast<mwSize>(vec.size());
        const mxClassID classid = get_mxClassID<T>();
        mxArray* arr = mxCreateNumericMatrix(n, 1, classid, mxREAL);
        T* ptr = static_cast<T*>(mxGetData(arr));
        std::copy(vec.begin(), vec.end(), ptr);
        return arr;
    }
    template<typename T>
    mxArray* vector_to_mxArray(const Vector<T>& vec, const size_t n) {
        const mxClassID classid = get_mxClassID<T>();
        mxArray* arr = mxCreateNumericMatrix(static_cast<mwSize>(n), 1, classid, mxREAL);
        T* ptr = static_cast<T*>(mxGetData(arr));
        std::copy(vec.begin(), vec.begin()+n, ptr);
        return arr;
    }

    std::string mxArray_to_string(const mxArray* arr);
    mxArray* string_to_mxArray(const std::string str);

    MxStruct mxArray_to_mxStruct(const mxArray* s);
    mxArray* mxStruct_to_mxArray(const MxStruct& mxStruct);
    MxStruct options_to_mxStruct(const Options& uno_options);
    Options mxStruct_to_options(const MxStruct& options);
    MxStruct result_to_mxStruct(const Result& uno_result);

} // namespace

#endif // UNOMEX_CONVERSION_H