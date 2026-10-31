// Copyright (c) 2025-2026 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_UNOPY_H
#define UNO_UNOPY_H

#include <functional>
#include <optional>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <vector>
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
#include "../UserModel.hpp"

namespace py = pybind11;

namespace uno {
   // typed, C-contiguous numpy array, used for all arrays exchanged with Python (a single type avoids
   // pybind11 re-wrapping an array_t of one flavor into another via PyArray_FromAny at each callback)
   template <typename T>
   using ndarray = py::array_t<T, py::array::c_style | py::array::forcecast>;

   // zero-copy numpy view of a buffer owned by C++. The array is created directly with numpy's C API: given a data
   // pointer, numpy does not own (and never frees) the memory, and the WRITEABLE flag is set at creation. This avoids
   // pybind11's generic constructor (heap-allocated shape/strides, copy unless a base object is given) and a
   // Python-level "flags.writeable = False" per callback
   template <typename T>
   ndarray<T> make_view(const T* data, size_t size, bool writeable) {
      using api = py::detail::npy_api;
      auto& numpy = api::get();
      Py_intptr_t shape = static_cast<Py_intptr_t>(size);
      // an empty buffer may have a null data pointer, which would make numpy allocate memory and ignore the flags
      static const T dummy{};
      PyObject* array = numpy.PyArray_NewFromDescr_(numpy.PyArray_Type_, py::dtype::of<T>().release().ptr(), 1, &shape,
         nullptr, const_cast<T*>(data != nullptr ? data : &dummy),
         api::NPY_ARRAY_C_CONTIGUOUS_ | (writeable ? api::NPY_ARRAY_WRITEABLE_ : 0), nullptr);
      if (array == nullptr) {
         throw py::error_already_set();
      }
      return py::reinterpret_steal<ndarray<T>>(array);
   }

   // input view (read-only) and output view (writeable)
   template <typename T>
   ndarray<T> input_view(const T* data, size_t size) { return make_view(data, size, false); }
   template <typename T>
   ndarray<T> output_view(T* data, size_t size) { return make_view(data, size, true); }

   // 1D array-like argument (numpy array, list, tuple, ...) converted to a std::vector by the caster below
   template <typename T>
   struct ArrayLike {
      std::vector<T> values;
   };

   template <typename T>
   std::vector<T> to_vector(ArrayLike<T>& array, size_t expected_size, const char* name) {
      if (array.values.size() != expected_size) {
         throw std::invalid_argument(std::string("Dimension mismatch in ") + name + ".");
      }
      return std::move(array.values);
   }
} // namespace

namespace pybind11::detail {
   template <typename T>
   struct type_caster<uno::ArrayLike<T>> {
      PYBIND11_TYPE_CASTER(uno::ArrayLike<T>, const_name("ArrayLike[") + make_caster<T>::name + const_name("]"));

      bool load(handle source, bool convert) {
         // numpy arrays: bulk conversion by numpy, much faster than an element-wise conversion
         if (isinstance<array>(source)) {
            if (!convert && !uno::ndarray<T>::check_(source)) {
               return false;
            }
            // like the conversion of Python sequences, reject non-integer values (e.g. floats) for integer arrays
            if constexpr (std::is_integral_v<T>) {
               const char kind = reinterpret_borrow<array>(source).dtype().kind();
               if (kind != 'i' && kind != 'u') {
                  return false;
               }
            }
            const auto array = uno::ndarray<T>::ensure(source);
            if (!array || array.ndim() != 1) {
               return false;
            }
            this->value.values.assign(array.data(), array.data() + array.size());
            return true;
         }
         // Python sequences: pybind11's element-wise conversion is faster than numpy's on lists
         make_caster<std::vector<T>> sequence_caster;
         if (!sequence_caster.load(source, convert)) {
            return false;
         }
         this->value.values = cast_op<std::vector<T>&&>(std::move(sequence_caster));
         return true;
      }
   };
} // namespace pybind11::detail

namespace uno {
   using Objective = std::function<double(ndarray<double>)>;

   using Constraints = std::function<void(ndarray<double>, ndarray<double>)>;

   using ObjectiveGradient = std::function<void(ndarray<double>, ndarray<double>)>;

   using Jacobian = std::function<void(ndarray<double>, ndarray<double>)>;

   using Hessian = std::function<void(ndarray<double>, double, ndarray<double>, ndarray<double>)>;

   using JacobianOperator = std::function<void(ndarray<double>, bool, ndarray<double>, ndarray<double>)>;

   using JacobianTransposedOperator = std::function<void(ndarray<double>, bool, ndarray<double>, ndarray<double>)>;

   using HessianOperator = std::function<void(ndarray<double>, bool, double, ndarray<double>,
      ndarray<double>, ndarray<double>)>;

   using NotifyAcceptableIterateCallback = std::function<void(ndarray<double> /*primals*/,
      ndarray<double> /*lower_bound_multipliers*/, ndarray<double> /*upper_bound_multipliers*/,
      ndarray<double> /*constraint_multipliers*/, double /*objective_multiplier*/, double /*primal_feasibility_residual*/,
      double /*stationarity_residual*/, double /*complementarity_residual*/)>;

   using TerminationCallback = std::function<bool(ndarray<double> /*primals*/, ndarray<double> /*lower_bound_multipliers*/,
      ndarray<double> /*upper_bound_multipliers*/, ndarray<double> /*constraint_multipliers*/, double /*objective_multiplier*/,
      double /*primal_feasibility_residual*/, double /*stationarity_residual*/, double /*complementarity_residual*/)>;

   using PythonUserModel = UserModel<std::optional<Objective>, std::optional<ObjectiveGradient>, std::optional<Constraints>,
      std::optional<Jacobian>, std::optional<JacobianOperator>, std::optional<JacobianTransposedOperator>,
      std::optional<Hessian>, std::optional<HessianOperator>, std::vector<double>, std::optional<py::object>>;
} // namespace

#endif // UNO_UNOPY_H