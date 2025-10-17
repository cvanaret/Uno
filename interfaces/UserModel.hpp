// Copyright (c) 2025 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_USERMODEL_H
#define UNO_USERMODEL_H

#include <iostream>
#include <streambuf>
#include <vector>
#include <cstring>
#include "C/Uno_C_API.h"

namespace uno {
   // UserModel contains the description of the model provided by the user
   template <typename Objective, typename ObjectiveGradient, typename Constraints, typename Jacobian,
      typename JacobianOperator, typename JacobianTransposedOperator, typename Hessian, typename HessianOperator,
      typename DoubleVector, typename UserDataType>
   class UserModel {
   public:
      UserModel(char problem_type, int32_t number_variables, int32_t base_indexing):
            problem_type(problem_type),
            base_indexing(base_indexing),
            number_variables(number_variables) {
      }

      ~UserModel() = default;

      const char problem_type; // 'L' for linear, 'Q' for quadratic, 'N' for nonlinear
      const int32_t base_indexing; // 0 for C-style indexing, 1 for Fortran-style indexing

      // variables
      const int32_t number_variables;
      DoubleVector variables_lower_bounds{};
      DoubleVector variables_upper_bounds{};

      // objective
      Objective objective_function{nullptr};
      ObjectiveGradient objective_gradient{nullptr};

      // constraints
      int32_t number_constraints{0};
      Constraints constraint_functions{nullptr};
      DoubleVector constraints_lower_bounds{};
      DoubleVector constraints_upper_bounds{};
      int32_t number_jacobian_nonzeros{0};
      std::vector<int32_t> jacobian_row_indices{};
      std::vector<int32_t> jacobian_column_indices{};
      Jacobian constraint_jacobian{nullptr};
      JacobianOperator jacobian_operator{nullptr};
      JacobianTransposedOperator jacobian_transposed_operator{nullptr};

      // Hessian
      int32_t number_hessian_nonzeros{0};
      // lower ('L') or upper ('U')
      char hessian_triangular_part{}; // default is empty
      std::vector<int32_t> hessian_row_indices{};
      std::vector<int32_t> hessian_column_indices{};
      Hessian lagrangian_hessian{nullptr};
      HessianOperator lagrangian_hessian_operator{nullptr};
      double lagrangian_sign_convention{UNO_MULTIPLIER_NEGATIVE};

      // User data
      UserDataType user_data{};

      // Optimization sense
      int32_t optimization_sense{UNO_MINIMIZE};

      // initial iterate
      DoubleVector initial_primal_iterate{};
      DoubleVector initial_dual_iterate{};
   };

   // base class for user stream callback
   class UserStreamCallback {
   public:
      UserStreamCallback() = default;
      virtual ~UserStreamCallback() = default;
      virtual int32_t operator()(const char* buf, int32_t len) const = 0;
   };

   // std::streambuf wrapper around UserStreamCallback
   class UserStreamBuffer : public std::streambuf {
   public:
      UserStreamBuffer(UserStreamCallback* user_stream_callback, std::size_t buffer_size) :
         user_stream_callback(user_stream_callback) {
         // allocate output buffer and set stream buffer pointer
         this->buffer = new char[buffer_size];
         this->setp(this->buffer, this->buffer + buffer_size - 1);
      }
      ~UserStreamBuffer() override {
         // flush remaining data and release buffer memory
         this->sync();
         // user_stream_callback has to be deleted explictly
         delete[] this->buffer;
      }

   protected:
      // called on buffer overflow
      int overflow(int character = EOF) override {
         if (character != EOF) {
               // insert the character into the buffer
               *this->pptr() = traits_type::to_char_type(character);
               this->pbump(1);
         }
         // return EOF for error
         return (this->flush_buffer() == 0) ? character : EOF;
      }

      // put char into buffer with \n triggering flush
      std::streamsize xsputn(const char* s, std::streamsize n) override {
         const char* p = s;
         const char* end = s + n;
         while (p < end) {
            const char* pos = static_cast<const char*>(std::memchr(p, '\n', end - p));
            if (pos) {
                  std::ptrdiff_t chunk_len = pos - p + 1;
                  if (this->epptr() - this->pptr() < chunk_len)
                     this->overflow();
                  std::memcpy(this->pptr(), p, chunk_len);
                  this->pbump(static_cast<int>(chunk_len));
                  this->overflow();
                  p = pos + 1;
            } else {
                  std::ptrdiff_t chunk_len = end - p;
                  if (this->epptr() - this->pptr() < chunk_len) {
                     overflow();
                     if (chunk_len > (this->epptr() - this->pptr()))
                        chunk_len = this->epptr() - this->pptr();
                  }
                  std::memcpy(pptr(), p, chunk_len);
                  this->pbump(static_cast<int>(chunk_len));
                  p += chunk_len;
            }
         }
         return n;
      }

      int sync() override {
         return this->flush_buffer();
      }

   private:
      UserStreamCallback* user_stream_callback;
      char* buffer;

      // flush buffer to the user_stream_callback
      int flush_buffer() {
         // check for invalid stream callback
         if (!this->user_stream_callback) {
            return -1;
         }
         std::ptrdiff_t current_used_buffer_size = this->pptr() - this->pbase();
         if (current_used_buffer_size > 0) {
            // call user stream callback
            const int32_t callback_result = (*this->user_stream_callback)(this->pbase(), static_cast<int32_t>(current_used_buffer_size));
            if (callback_result != static_cast<int32_t>(current_used_buffer_size)) {
               return -1;
            }
            // move buffer pointer
            this->pbump(static_cast<int>(-current_used_buffer_size));
         }
         return 0;
      }
   };

   // std::ostream wrapper around UserStreamCallback
   class UserOStream : public std::ostream {
   public:
      UserOStream(UserStreamCallback* user_stream_callback, std::size_t buffer_size = 1024) : // 1024 default buffer size, sufficient to store the entire line
         std::ostream(&this->buffer), buffer(user_stream_callback, buffer_size) { }
      
   private:
      // internal stream buffer that sends output to the LoggerStreamUserCallback
      UserStreamBuffer buffer;
   };

} // namespace

#endif // UNO_USERMODEL_H