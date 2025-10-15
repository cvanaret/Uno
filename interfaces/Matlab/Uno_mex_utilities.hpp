// Copyright (c) 2025 Stefano Lovato and Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_MEX_UTILITIES_H
#define UNO_MEX_UTILITIES_H

#include <stdint.h>
#include <map>
#include <string>
#include <algorithm>
#include "mex.h"
#include "Uno.hpp"
#include "linear_algebra/SparseVector.hpp"
#include "linear_algebra/Vector.hpp"
#include "options/Options.hpp"
#include "tools/Logger.hpp"
#include "tools/Logger.hpp"

namespace uno {
    // map C++ type to MATLAB class ID
    template<typename T>
    mxClassID getMxClassID();
    template<> mxClassID getMxClassID<double>()  { return mxDOUBLE_CLASS;  }
    template<> mxClassID getMxClassID<float>()   { return mxSINGLE_CLASS;  }
    template<> mxClassID getMxClassID<int8_t>()  { return mxINT8_CLASS;    }
    template<> mxClassID getMxClassID<uint8_t>() { return mxUINT8_CLASS;   }
    template<> mxClassID getMxClassID<int16_t>() { return mxINT16_CLASS;   }
    template<> mxClassID getMxClassID<uint16_t>(){ return mxUINT16_CLASS;  }
    template<> mxClassID getMxClassID<int32_t>() { return mxINT32_CLASS;   }
    template<> mxClassID getMxClassID<uint32_t>(){ return mxUINT32_CLASS;  }
    template<> mxClassID getMxClassID<int64_t>() { return mxINT64_CLASS;   }
    template<> mxClassID getMxClassID<uint64_t>(){ return mxUINT64_CLASS;  }
    template<> mxClassID getMxClassID<bool>()    { return mxLOGICAL_CLASS; }
    template<> mxClassID getMxClassID<char>()    { return mxCHAR_CLASS;    }

    // std::map to store matlab struct
    class MxStruct {
    public:
        mxArray* operator[](const std::string& key) const {
            auto it = data.find(key);
            return (it != data.end()) ? it->second : nullptr;
        }
        void insert(const std::string& key, mxArray* value) {
            data[key] = value;
        }
        bool contains(const std::string& key) const {
            return data.find(key) != data.end();
        }
        auto begin() const { return data.begin(); }
        auto end() const { return data.end(); }
    private:
        std::map<std::string, mxArray*> data;
    };

    // call matlab function with error trapping
    int32_t call_matlab_function(mxArray* handle, const std::vector<mxArray*>& inputs, std::vector<mxArray*>& outputs) {
        std::vector<mxArray*> inputs_all = inputs;
        inputs_all.insert(inputs_all.begin(), handle);
        const int nrhs = static_cast<int>(inputs_all.size());
        const int nlhs = static_cast<int>(outputs.size());
        mxArray** prhs = inputs_all.data();
        mxArray** plhs = outputs.data();
        mxArray* err = mexCallMATLABWithTrap(nlhs, plhs, nrhs, prhs, "feval");
        if (err) {
            // display the error message but do not terminate the mex function
            mxArray* errmsg;
            mexCallMATLAB(1, &errmsg, 1, &err, "getReport"); // get err message from MException 
            mexWarnMsgIdAndTxt("Uno:MatlabError", "%s", mxArrayToString(errmsg)); // trow warning but do not terminate
            mxDestroyArray(err);
            mxDestroyArray(errmsg);
            return 1;
        } else {
            return 0;
        }
    }

    // convert Vector type
    template <typename OutType, typename InType>
    Vector<OutType> convert_vector_type(const Vector<InType>& input) {
        Vector<OutType> output(input.size());
        std::transform(input.begin(), input.end(), output.begin(),
                    [](const InType& val) { return static_cast<OutType>(val); });
        return output;
    }

    // destroy mxArray vector
    void destroy_mxArray_vector(std::vector<mxArray*>& vector) {
        for (auto* arr : vector) {
            mxDestroyArray(arr);
        }
    }

    // copy mxArray to raw pointer
    template<typename T>
    void mxArray_to_pointer(const mxArray* array, T* pointer) {
        size_t n = mxGetNumberOfElements(array);
        T* ptr = static_cast<T*>(mxGetData(array));
        std::copy(ptr, ptr+n, pointer);
    }

    // copy raw pointer to mxArray
    template<typename T>
    mxArray* pointer_to_mxArray(const T* pointer, size_t n) {
        mxArray* array = mxCreateNumericMatrix(n, 1, getMxClassID<T>(), mxREAL);
        T* ptr = static_cast<T*>(mxGetData(array));
        std::copy(pointer, pointer+n, ptr);
        return array;
    }

    // convert mxArray to scalar
    template<typename T>
    T mxArray_to_scalar(const mxArray* array) {
        return *(static_cast<T*>(mxGetData(array)));
    }

    // convert scalar to mxArray
    template<typename T>
    mxArray* scalar_to_mxArray(const T value) {
        mxArray* arr = mxCreateNumericMatrix(1, 1, getMxClassID<T>(), mxREAL);
        T* ptr = static_cast<T*>(mxGetData(arr));
        *ptr = value;
        return arr;
    }

    // convert mxArray to string
    std::string mxArray_to_string(const mxArray* array) {
        std::string str;
        if (!array) {
            return str;
        }
        else if (mxIsChar(array)) {
            char* cstr = mxArrayToString(array);
            str = std::string(cstr);
            mxFree(cstr);
            return str;
        }
        else if (mxIsClass(array,"string")) {
            // convert matlab string to char and call mxArray_to_string
            // TODO handle non-scalar matlab string
            mxArray* array_str = mxDuplicateArray(array);
            mxArray* array_char = nullptr;
            int return_value = mexCallMATLAB(1, &array_char, 1, &array_str, "char");
            if (return_value == 0) {
                str = mxArray_to_string(array_char);
            }
            mxDestroyArray(array_str);
            mxDestroyArray(array_char);
            return str;
        }
        else {
            return str;
        }
    }

    // convert string to mxArray
    mxArray* string_to_mxArray(const std::string str) {
        return mxCreateString(str.c_str());
    }

    // convert mxArray to Uno::Vector
    template<typename T>
    void mxArray_to_vector(const mxArray* array, Vector<T>& vec) {
        size_t n = mxGetNumberOfElements(array); // assume vec has the correct size
        T* ptr = mxGetPr(array);
        std::copy(ptr, ptr+n, vec.data());
    }
    template<typename T>
    Vector<T> mxArray_to_vector(const mxArray* array) {
        size_t n = mxGetNumberOfElements(array);
        Vector<T> vec(n);
        mxArray_to_vector(array, vec);
        return vec;
    }

    // convert Uno::Vector to mxArray
    template<typename T>
    mxArray* vector_to_mxArray(const Vector<T>& vec) {
        mwSize n = static_cast<mwSize>(vec.size());
        mxArray* arr = mxCreateNumericMatrix(n, 1, getMxClassID<T>(), mxREAL);
        T* ptr = static_cast<T*>(mxGetData(arr));
        std::copy(vec.begin(), vec.end(), ptr);
        return arr;
    }
    
    // convert MxStruct to mxArray
    MxStruct mxArray_to_mxStruct(const mxArray* s) {
        if (!s) return MxStruct();
        MxStruct result;
        int nFields = mxGetNumberOfFields(s);
        for (int i = 0; i < nFields; ++i) {
            const char* fieldName = mxGetFieldNameByNumber(s, i);
            mxArray* fieldValue = mxGetFieldByNumber(s, 0, i);
            result.insert(fieldName, fieldValue);
        }
        return result;
    }

    // convert MxStruct to mxArray
    mxArray* mxStruct_to_mxArray(const MxStruct& mxStruct) {
        Vector<const char*> fieldNames;
        for (const auto& [key, value]: mxStruct) {
            fieldNames.push_back(key.c_str());
        }
        mxArray* s = mxCreateStructMatrix(1, 1, static_cast<int>(fieldNames.size()), fieldNames.data());
        for (const auto& [key, value]: mxStruct) {
            mxSetField(s, 0, key.c_str(), value);
        }
        return s;
    }

    // convert Uno::options to MxStruct
    MxStruct options_to_mxStruct(const Options& uno_options) {
        MxStruct mxStruct;
        // integer options (matlab double)
        std::map<std::string, int32_t> integer_options = uno_options.get_integer_options();
        for (const auto& [option_name, option_value]: integer_options) {
            mxStruct.insert(option_name, scalar_to_mxArray<double>(static_cast<double>(option_value)));
        }
        // double options (matlab double)
        std::map<std::string, double> double_options = uno_options.get_double_options();
        for (const auto& [option_name, option_value]: double_options) {
            mxStruct.insert(option_name, scalar_to_mxArray<double>(option_value));
        }
        // bool options (matlab logical)
        std::map<std::string, bool> bool_options = uno_options.get_bool_options();
        for (const auto& [option_name, option_value]: bool_options) {
            mxStruct.insert(option_name, scalar_to_mxArray<bool>(option_value));
        }
        // string options (matlab char)
        std::map<std::string, std::string> string_options = uno_options.get_string_options();
        for (const auto& [option_name, option_value]: string_options) {
            mxStruct.insert(option_name, string_to_mxArray(option_value));
        }        
        return mxStruct;
    }

    // convert MxStruct to Uno::Options
    Options mxStruct_to_options(const MxStruct& options) {
        Options uno_options;
        for (const auto& [field_name, field_value]: options) {
            try {
                // set option with check between type and mxArray class
                OptionType option_type = uno_options.get_option_type(field_name);
                // OptionType::DOUBLE accepts double
                if (option_type == OptionType::DOUBLE) {
                    if (mxGetClassID(field_value) == getMxClassID<double>()) {
                        uno_options.set_double(field_name, mxArray_to_scalar<double>(field_value));
                    }
                    else {
                        INFO << "Incorrect type of option '" << field_name << "'." << std::endl;
                    }
                }
                // OptionType::INTEGER accepts double
                else if (option_type == OptionType::INTEGER) {
                    if (mxGetClassID(field_value) == getMxClassID<double>()) {
                        uno_options.set_integer(field_name, static_cast<int32_t>(mxArray_to_scalar<double>(field_value))); // cast to int
                    }
                    else {
                        INFO << "Incorrect type of option '" << field_name << "'." << std::endl;
                    }
                }
                // OptionType::BOOL accepts logical
                else if (option_type == OptionType::BOOL) {
                    if (mxGetClassID(field_value) == getMxClassID<bool>()) { 
                        uno_options.set_bool(field_name, mxArray_to_scalar<bool>(field_value)); 
                    }
                    else {
                        INFO << "Incorrect type of option '" << field_name << "'." << std::endl;
                    }
                }
                // OptionType::STRING accepts char
                else if (option_type == OptionType::STRING) {
                    if (mxGetClassID(field_value) == getMxClassID<char>()) { 
                        uno_options.set_string(field_name, mxArray_to_string(field_value)); 
                    }
                    else {
                        INFO << "Incorrect type of option '" << field_name << "'." << std::endl;
                    }
                }
            } catch (const std::out_of_range&) {
                // set the option with type depending on mxArray class (only double, bool, char)
                if (mxGetClassID(field_value) == getMxClassID<double>()) {
                    uno_options.set_double(field_name, mxArray_to_scalar<double>(field_value));
                }
                else if (mxGetClassID(field_value) == getMxClassID<bool>()) {
                    uno_options.set_bool(field_name, mxArray_to_scalar<bool>(field_value));
                }
                else if (mxGetClassID(field_value) == getMxClassID<char>()) {
                    uno_options.set_string(field_name, mxArray_to_string(field_value));
                }
                else {
                    INFO << "Type of option '" << field_name << "' is not allowed. Allowed types are: double, int32, logical, char" << std::endl;
                }
            }
        }
        return uno_options;
    }
    
    // convert Uno::Result to MxStruct
    MxStruct result_to_mxStruct(const Result& uno_result) {
        MxStruct result;
        result.insert("optimization_status", scalar_to_mxArray(static_cast<double>(uno_result.optimization_status)));
        result.insert("solution_status", scalar_to_mxArray(static_cast<double>(uno_result.solution_status)));
        result.insert("solution_objective", scalar_to_mxArray(uno_result.solution_objective));
        result.insert("solution_primal_feasibility", scalar_to_mxArray(uno_result.solution_primal_feasibility));
        result.insert("solution_stationarity", scalar_to_mxArray(uno_result.solution_stationarity));
        result.insert("solution_complementarity", scalar_to_mxArray(uno_result.solution_complementarity));
        result.insert("primal_solution", vector_to_mxArray(uno_result.primal_solution));
        result.insert("constraint_dual_solution", vector_to_mxArray(uno_result.constraint_dual_solution));
        result.insert("lower_bound_dual_solution", vector_to_mxArray(uno_result.lower_bound_dual_solution));
        result.insert("upper_bound_dual_solution", vector_to_mxArray(uno_result.upper_bound_dual_solution));
        result.insert("number_iterations", scalar_to_mxArray(static_cast<double>(uno_result.number_iterations)));
        result.insert("cpu_time", scalar_to_mxArray(uno_result.cpu_time));
        result.insert("number_objective_evaluations", scalar_to_mxArray(static_cast<double>(uno_result.number_objective_evaluations)));
        result.insert("number_constraint_evaluations", scalar_to_mxArray(static_cast<double>(uno_result.number_constraint_evaluations)));
        result.insert("number_objective_gradient_evaluations", scalar_to_mxArray(static_cast<double>(uno_result.number_objective_gradient_evaluations)));
        result.insert("number_jacobian_evaluations", scalar_to_mxArray(static_cast<double>(uno_result.number_jacobian_evaluations)));
        result.insert("number_hessian_evaluations", scalar_to_mxArray(static_cast<double>(uno_result.number_hessian_evaluations)));
        result.insert("number_subproblems_solved", scalar_to_mxArray(static_cast<double>(uno_result.number_subproblems_solved)));
        return result;
    }
    

};

#endif //UNO_MEX_UTILITIES_H