// Copyright (c) 2025 Stefano Lovato and Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_MEX_UTILITIES_H
#define UNO_MEX_UTILITIES_H

#include <stdint.h>
#include <map>
#include <string>
#include <algorithm>
#include <cstdio>
#include <cstdarg>
#include "mex.h"
#include "Uno.hpp"
#include "optimization/EvaluationErrors.hpp"
#include "linear_algebra/SparseVector.hpp"
#include "linear_algebra/Vector.hpp"
#include "options/Options.hpp"
#include "tools/Logger.hpp"
#include "tools/Logger.hpp"

namespace uno {
    
    // matlab error string formatter
    enum class ErrorType {
        NARGIN_NOTENOUGH,
        NARGIN_TOOMANY,
        NARGOUT_NOTENOUGH,
        NARGOUT_TOOMANY,
        INPUT_INVALID,
        INPUT_TYPE_STRING,
        INPUT_TYPE_STRUCT,
        INPUT_SCALAR,
        OUTPUT_TYPE,
        OUTPUT_SIZE,
        INVALID_HANDLE,
        EVALUATION,
        OPTION_TYPE,
        MISSING_FIELD,
        FIELD_CHAR,
        FIELD_POSITIVE_INTEGER,
        FIELD_UNITARY,
        FIELD_VECTOR,
        FIELD_VECTOR_INT,
        UNO
    }; 
    class ErrorString {
    public:
        inline static std::string format_error(ErrorType error_type, ...) {
            auto it = error_strings.find(error_type);
            if (it == error_strings.end()) {
                return "Unknown error.";
            }
            // format the string
            const std::string& error_fmt = it->second;
            va_list args1;
            va_start(args1, error_type);
            va_list args2;
            va_copy(args2, args1);
            int size = std::vsnprintf(nullptr, 0, error_fmt.c_str(), args2);
            va_end(args2);
            std::string error_msg(size, '\0');
            std::vsnprintf(error_msg.data(), error_msg.size() + 1, error_fmt.c_str(), args1);
            va_end(args1);
            return error_msg;
        }
    private:
        inline static const std::map<ErrorType, std::string> error_strings = {
            {ErrorType::NARGIN_NOTENOUGH, "Invalid argument list. Function requires %d more input(s)."},
            {ErrorType::NARGIN_TOOMANY, "Too many input arguments."},
            {ErrorType::NARGOUT_NOTENOUGH, "Invalid argument list. Function requires %d more output(s)."},
            {ErrorType::NARGOUT_TOOMANY, "Too many output arguments."},
            {ErrorType::INPUT_INVALID, "Invalid input argument at position %d."},
            {ErrorType::INPUT_TYPE_STRING, "Invalid argument at position %d. Value must be of type char or string."},
            {ErrorType::INPUT_TYPE_STRUCT, "Invalid argument at position %d. Value must be of type struct."},  
            {ErrorType::INPUT_SCALAR, "Invalid input argument at position %d. Value must be a scalar."},
            {ErrorType::OUTPUT_TYPE, "Invalid type of output argument at position %d."},
            {ErrorType::OUTPUT_SIZE, "Invalid size of output argument at position %d."},
            {ErrorType::INVALID_HANDLE, "Value must be a function handle."},
            {ErrorType::EVALUATION, "Error in function evaluation."},
            {ErrorType::OPTION_TYPE, "Invalid type of option '%s'."},
            {ErrorType::MISSING_FIELD, "Missing field '%s'."},
            {ErrorType::FIELD_CHAR, "Type of field '%s' must be char."},
            {ErrorType::FIELD_POSITIVE_INTEGER, "Field '%s' must be a positive integer value."},
            {ErrorType::FIELD_UNITARY, "Field '%s' must be either +1 or -1."},
            {ErrorType::FIELD_VECTOR, "Field '%s' must be a vector of %d element(s)."},
            {ErrorType::FIELD_VECTOR_INT, "Field '%s' must be a vector of %d integer element(s)."},
            {ErrorType::UNO, "Error in UNO. %s"}
        };
    };

    // matlab error
    struct MatlabFunctionError : EvaluationError {
        std::string msg;
        explicit MatlabFunctionError(const std::string& str) : msg(ErrorString::format_error(ErrorType::EVALUATION) + "\n" + str + "\n") {}
        [[nodiscard]] const char* what() const noexcept override {
            return this->msg.c_str();
        }
    };

    // std::map to store matlab struct
    class MxStruct {
    public:
        mxArray* operator[](const std::string& key) const {
            auto it = this->data.find(key);
            return (it != this->data.end()) ? it->second : nullptr;
        }
        void insert(const std::string& key, mxArray* value) {
            this->data[key] = value;
        }
        bool contains(const std::string& key) const {
            return this->data.find(key) != this->data.end();
        }
        auto begin() const { return this->data.begin(); }
        auto end() const { return this->data.end(); }
    private:
        std::map<std::string, mxArray*> data;
    };
    
    // matlab typedef for Matlab function
    typedef mxArray* handle_t;    

    // map C++ type to MATLAB class ID
    #define mxSTRING_CLASS 19 // cf. https://it.mathworks.com/matlabcentral/answers/2102376-how-do-i-process-a-string-class-in-a-mex-function
    template<typename T>
    mxClassID get_mxClassID();
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

    // validate mxArray pointer
    bool isvalid(const mxArray* arr) {
        return (arr != nullptr);
    }
    
    // validate mxArray empty
    bool isempty(const mxArray* arr) {
        return mxIsEmpty(arr);
    }

    // validate mxArray type
    template<typename T>
    bool isa(const mxArray* arr) {
        return mxGetClassID(arr) == get_mxClassID<T>();
    }

    // validate mxArray size (nrows, ncolumns)
    bool has_size(const mxArray* arr, const int32_t nrows, const int32_t ncolumns) {
        if (mxGetNumberOfDimensions(arr)>2) {
            return false;;
        }
        return (static_cast<int32_t>(mxGetM(arr))==nrows) && (static_cast<int32_t>(mxGetN(arr))==ncolumns);
    }

    // validate mxArray positive
    bool ispositive(const mxArray* arr) {
        const double value = mxGetScalar(arr);
        return value>=0;
    }

    // validate mxArray double integer
    bool isinteger(const mxArray* arr) {
        const double value = mxGetScalar(arr);
        if (!std::isfinite(value)) {
            return false;
        }
        double intpart;
        return std::modf(value, &intpart) == 0.0;
    }

    // validate mxArray unitary
    bool isunitary(const mxArray* arr) {
        const double value = mxGetScalar(arr);
        return (value==1.0) || (value==-1.0);
    }

    // validate mxArray scalar
    bool isscalar(const mxArray* arr) {
        return mxGetNumberOfElements(arr)==1;
    }

    // validate mxArray vector size (nrows, 1) or (1, ncolumns)
    bool isvector(const mxArray* arr, const int32_t len) {
        return has_size(arr, len, 1) || has_size(arr, 1, len);
    }

    // call matlab function with error trapping
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
            mexCallMATLAB(1, &errmsg, 1, &err, "getReport"); // get err message from MException 
            const std::string strerr = mxArrayToString(errmsg);
            mxDestroyArray(err);
            mxDestroyArray(errmsg);
            throw MatlabFunctionError(strerr);
            
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
    void mxArray_to_pointer(const mxArray* arr, T* pointer) {
        const size_t n = mxGetNumberOfElements(arr);
        const T* ptr = static_cast<T*>(mxGetData(arr));
        std::copy(ptr, ptr+n, pointer);
    }

    // copy raw pointer to mxArray
    template<typename T>
    mxArray* pointer_to_mxArray(const T* pointer, const size_t n) {
        const mxClassID classid = get_mxClassID<T>();
        mxArray* arr = mxCreateNumericMatrix(n, 1, classid, mxREAL);
        T* ptr = static_cast<T*>(mxGetData(arr));
        std::copy(pointer, pointer+n, ptr);
        return arr;
    }

    // convert mxArray to scalar
    template<typename T>
    T mxArray_to_scalar(const mxArray* arr) {
        return *(static_cast<T*>(mxGetData(arr)));
    }

    // convert scalar to mxArray
    template<typename T>
    mxArray* scalar_to_mxArray(const T value) {
        const mxClassID classid = get_mxClassID<T>();
        mxArray* arr = mxCreateNumericMatrix(1, 1, classid, mxREAL);
        T* ptr = static_cast<T*>(mxGetData(arr));
        *ptr = value;
        return arr;
    }

    // convert mxArray to string
    std::string mxArray_to_string(const mxArray* arr) {
        std::string str;
        if (!isvalid(arr)) {
            return str;
        }
        else if (mxIsChar(arr)) {
            char* cstr = mxArrayToString(arr);
            str = std::string(cstr);
            mxFree(cstr);
            return str;
        }
        else if (mxIsClass(arr,"string")) {
            // convert matlab string to char and call mxArray_to_string
            // TODO handle non-scalar matlab string
            mxArray* arr_str = mxDuplicateArray(arr);
            mxArray* arr_char;
            mexCallMATLAB(1, &arr_char, 1, &arr_str, "char");
            str = mxArray_to_string(arr_char);
            mxDestroyArray(arr_str);
            mxDestroyArray(arr_char);
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

    // convert Uno::Vector to mxArray
    template<typename T>
    mxArray* vector_to_mxArray(const Vector<T>& vec) {
        const mwSize n = static_cast<mwSize>(vec.size());
        const mxClassID classid = get_mxClassID<T>();
        mxArray* arr = mxCreateNumericMatrix(n, 1, classid, mxREAL);
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
                const OptionType option_type = uno_options.get_option_type(field_name);
                // OptionType::DOUBLE accepts double
                if (option_type == OptionType::DOUBLE) {
                    if (isa<double>(field_value)) {
                        uno_options.set_double(field_name, mxArray_to_scalar<double>(field_value));
                    }
                    else {
                        INFO << ErrorString::format_error(ErrorType::OPTION_TYPE, field_name.c_str()) << std::endl;
                    }
                }
                // OptionType::INTEGER accepts double
                else if (option_type == OptionType::INTEGER) {
                    if (isa<double>(field_value)) {
                        uno_options.set_integer(field_name, static_cast<int32_t>(mxArray_to_scalar<double>(field_value))); // cast to int
                    }
                    else {
                        INFO << ErrorString::format_error(ErrorType::OPTION_TYPE, field_name.c_str()) << std::endl;
                    }
                }
                // OptionType::BOOL accepts logical
                else if (option_type == OptionType::BOOL) {
                    if (isa<bool>(field_value)) { 
                        uno_options.set_bool(field_name, mxArray_to_scalar<bool>(field_value)); 
                    }
                    else {
                        INFO << ErrorString::format_error(ErrorType::OPTION_TYPE, field_name.c_str()) << std::endl;
                    }
                }
                // OptionType::STRING accepts char
                else if (option_type == OptionType::STRING) {
                    if (isa<char>(field_value) || isa<std::string>(field_value)) { 
                        uno_options.set_string(field_name, mxArray_to_string(field_value)); 
                    }
                    else {
                        INFO << ErrorString::format_error(ErrorType::OPTION_TYPE, field_name.c_str()) << std::endl;
                    }
                }
            } catch (const std::out_of_range&) {
                // set the option with type depending on mxArray class (only double, bool, char)
                if (isa<double>(field_value)) {
                    uno_options.set_double(field_name, mxArray_to_scalar<double>(field_value));
                }
                else if (isa<bool>(field_value)) {
                    uno_options.set_bool(field_name, mxArray_to_scalar<bool>(field_value));
                }
                else if (isa<char>(field_value) || isa<std::string>(field_value)) {
                    uno_options.set_string(field_name, mxArray_to_string(field_value));
                }
                else {
                    INFO << ErrorString::format_error(ErrorType::OPTION_TYPE, field_name.c_str()) << std::endl;
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
    
    // get the number of input and output arguments of Matlab function handle
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

    // validate Matlab function handle with no call
    template <typename InType, typename OutType>
    bool validate_matlab_handle_no_call(handle_t handle, int32_t expected_nin, int32_t expected_nout,
        const std::string field_name, std::string& errmsg) {
        if (!isvalid(handle)) {
            errmsg = ErrorString::format_error(ErrorType::MISSING_FIELD, field_name.c_str());
            return false;
        }
        if (!isa<handle_t>(handle)) {
            errmsg = ErrorString::format_error(ErrorType::INVALID_HANDLE);
            return false;
        }
        // nargin
        const int32_t nin = nargin(handle);
        if (nin > expected_nin) {
            errmsg = ErrorString::format_error(ErrorType::NARGIN_TOOMANY);
            return false;
        } 
        else if (nin < expected_nin) {
            errmsg = ErrorString::format_error(ErrorType::NARGIN_NOTENOUGH, expected_nin-nin);
            return false;
        }
        // nargout
        const int32_t nout = nargout(handle);
        if (nout > expected_nout) {
            errmsg = ErrorString::format_error(ErrorType::NARGOUT_TOOMANY);
            return false;
        } 
        else if (nout < expected_nout) {
            errmsg = ErrorString::format_error(ErrorType::NARGOUT_NOTENOUGH, expected_nout-nout);
            return false;
        }
        return true;
    }

    // validate Matlab function handle
    template <typename InType, typename OutType>
    bool validate_matlab_handle(handle_t handle,
        const std::vector<int32_t> input_dims, const std::vector<int32_t> output_dims, 
        const std::string field_name, std::string& errmsg) {
        const int32_t nin = static_cast<int32_t>(input_dims.size());
        const int32_t nout = static_cast<int32_t>(output_dims.size());
        if (!validate_matlab_handle_no_call<InType, OutType>(handle, nin, nout, field_name, errmsg)) {
            return false;
        }
        // test the function call
        std::vector<mxArray*> inputs(nin);
        std::vector<mxArray*> outputs(nout);
        for (size_t input_index: Range(nin)) {
            const int32_t nrows = input_dims[input_index];
            const Vector<InType> input(nrows, 0.);
            inputs[input_index] = vector_to_mxArray(input);
        }
        try {
            call_matlab_function(handle, inputs, outputs);
            destroy_mxArray_vector(inputs);
        } catch (const MatlabFunctionError& err) {
            errmsg = err.what(); 
            return false;
        }
        // validate the outputs
        for (size_t output_index: Range(nout)) {
            const mxArray* output = outputs[output_index];
            const int32_t nrows = output_dims[output_index];
            if (!isa<OutType>(output)) {
                errmsg = ErrorString::format_error(ErrorType::OUTPUT_TYPE, output_index+1);
                destroy_mxArray_vector(outputs);
                return false;
            }
            if (!isvector(output, nrows)) {
                errmsg = ErrorString::format_error(ErrorType::OUTPUT_SIZE, output_index+1);
                destroy_mxArray_vector(outputs);
                return false;
            }
        }
        destroy_mxArray_vector(outputs);
        return true;
    }

    bool validate_struct_input(const mxArray* arr, const int32_t position, std::string& errmsg) {
        if (!isvalid(arr) || !isa<MxStruct>(arr)) {
            errmsg = ErrorString::format_error(ErrorType::INPUT_TYPE_STRUCT, position);
            return false;
        }
        if (!isscalar(arr)) {
            errmsg = ErrorString::format_error(ErrorType::INPUT_SCALAR, position);
            return false;
        }
        return true;
    }

    bool validate_char(const mxArray* arr, const std::vector<char> chars, const std::string field_name, std::string& errmsg) {
        if (!isvalid(arr)) {
            errmsg = ErrorString::format_error(ErrorType::MISSING_FIELD, field_name.c_str());
            return false;
        }
        if (!isa<char>(arr) || !isscalar(arr)) {
            errmsg = ErrorString::format_error(ErrorType::FIELD_CHAR, field_name.c_str());
            return false;
        }
        if (chars.size()>0) {
            char c = mxArray_to_scalar<char>(arr);
            if (std::find(chars.begin(), chars.end(), c) == chars.end()) {
                errmsg = "Field '" + field_name + "' must be one of { ";
                for (size_t i: Range(chars.size())) {
                    if (i > 0) {
                        errmsg += ", ";
                    }
                    errmsg += "'" + std::string({chars[i]}) + "'";
                }
                errmsg += " }.";
                return false;
            }
        }
        return true;
    }

    bool validate_positive_integer(const mxArray* arr, const std::string field_name, std::string& errmsg) {
        if (!isvalid(arr)) {
            errmsg = ErrorString::format_error(ErrorType::MISSING_FIELD, field_name.c_str());
            return false;
        }
        if (!isa<double>(arr) || !isscalar(arr) || !ispositive(arr)  || !isinteger(arr) ) {
            errmsg = ErrorString::format_error(ErrorType::FIELD_POSITIVE_INTEGER, field_name.c_str());
            return false;
        }
        return true;
    }

    bool validate_unitary(const mxArray* arr, const std::string field_name, std::string& errmsg) {
        if (!isvalid(arr)) {
            errmsg = ErrorString::format_error(ErrorType::MISSING_FIELD, field_name.c_str());
            return false;
        }
        if (!isa<double>(arr) || !isscalar(arr) || !isinteger(arr) || !isunitary(arr) ) {
            errmsg = ErrorString::format_error(ErrorType::FIELD_UNITARY, field_name.c_str());
            return false;
        }
        return true;
    }

    bool validate_double_vector(const mxArray* arr, const int32_t len, const std::string field_name, std::string& errmsg) {
        if (!isvalid(arr)) {
            errmsg = ErrorString::format_error(ErrorType::MISSING_FIELD, field_name.c_str());
            return false;
        }
        if (!isa<double>(arr) || !isvector(arr, len) ) {
            errmsg = ErrorString::format_error(ErrorType::FIELD_VECTOR, field_name.c_str(), len);
            return false;
        }
        return true;
    }

    bool validate_integer_vector(const mxArray* arr, const int32_t len, const std::string field_name, std::string& errmsg) {
        if (!isvalid(arr)) {
            errmsg = ErrorString::format_error(ErrorType::MISSING_FIELD, field_name.c_str());
            return false;
        }
        if (!isa<double>(arr) || !isvector(arr, len) ) {
            errmsg = ErrorString::format_error(ErrorType::FIELD_VECTOR_INT, field_name.c_str(), len);
            return false;
        }
        const Vector vector = mxArray_to_vector<double>(arr);
        double intpart;
        for (auto value : vector) {
            if (std::modf(value, &intpart) != 0.0) {
                errmsg = ErrorString::format_error(ErrorType::FIELD_VECTOR_INT, field_name.c_str(), len);
                return false;
            }
        }
        return true;
    }

};

#endif //UNO_MEX_UTILITIES_H