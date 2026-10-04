#include "tape.hpp"
#include <algorithm>
#include <cmath>
#include <limits>

namespace cppasl::detail {

   double evaluate_piecewise_linear(const double* data, double x, double* slope) {
      // data = {k, s_0, b_0, s_1, b_1, ..., b_{k-2}, s_{k-1}}: slope s_j on (b_{j-1}, b_j], b_{-1} = -inf, b_{k-1} = +inf
      const int number_slopes = static_cast<int>(data[0]);
      const double infinity = std::numeric_limits<double>::infinity();
      double value = 0.;
      for (int j = 0; j < number_slopes; ++j) {
         const double segment_slope = data[1 + 2 * j];
         const double lower = (j == 0) ? -infinity : data[2 * j];
         const double upper = (j == number_slopes - 1) ? infinity : data[2 + 2 * j];
         if (slope != nullptr && lower < x && x <= upper) {
            *slope = segment_slope;
         }
         // integrate the slope between 0 and x
         if (x >= 0.) {
            const double from = std::max(lower, 0.), to = std::min(upper, x);
            if (to > from) value += segment_slope * (to - from);
         }
         else {
            const double from = std::max(lower, x), to = std::min(upper, 0.);
            if (to > from) value -= segment_slope * (to - from);
         }
      }
      return value;
   }

   namespace {
      /// operand selected by min, max or if-then-else
      inline std::int32_t selected_operand(const TapeInstruction& instruction, const TapeView& tape, const double* values) {
         const std::int32_t* operands = tape.arguments + instruction.first;
         if (instruction.operation == TapeOperation::IfThenElse) {
            return (values[operands[0]] != 0.) ? operands[1] : operands[2];
         }
         std::int32_t best = operands[0];
         for (std::int32_t k = 1; k < instruction.second; ++k) {
            const double candidate = values[operands[k]];
            if ((instruction.operation == TapeOperation::Minimum) ? (candidate < values[best]) : (candidate > values[best])) {
               best = operands[k];
            }
         }
         return best;
      }

      inline double round_to_decimals(double x, double decimals) {
         const double scale = std::pow(10., decimals);
         return std::round(x * scale) / scale;
      }

      inline double truncate_to_decimals(double x, double decimals) {
         const double scale = std::pow(10., decimals);
         return std::trunc(x * scale) / scale;
      }

      inline double round_to_significant_digits(double x, double digits) {
         if (x == 0. || !std::isfinite(x)) return x;
         const double exponent = std::floor(std::log10(std::fabs(x)));
         return round_to_decimals(x, digits - 1. - exponent);
      }
   } // namespace

   namespace {
      inline void sine_and_cosine(double u, double& sine, double& cosine) {
#if defined(__GLIBC__)
         ::sincos(u, &sine, &cosine);
#else
         sine = std::sin(u);
         cosine = std::cos(u);
#endif
      }

      /// value of one instruction (all operations)
      inline double evaluate_instruction(const TapeView& tape, const TapeInstruction& instruction, const double* x,
            const double* values) {
         using O = TapeOperation;
         const double c = instruction.constant;
         const int arity = operand_arity(instruction.operation);
         const double u = (arity >= 1) ? values[instruction.first] : 0.;
         const double w = (arity == 2) ? values[instruction.second] : 0.;
         double result;
         switch (instruction.operation) {
            case O::Constant: result = c; break;
            case O::Variable: result = x[instruction.first] + c; break; // shifted variable leaf (offset c)
            case O::DefinedValue: result = std::numeric_limits<double>::quiet_NaN(); break; // only in evaluation tapes
            case O::Negate: result = -u; break;
            case O::AddConstant: result = u + c; break;
            case O::MultiplyConstant: result = c * u; break;
            case O::Square: result = u * u; break;
            case O::PowerInteger: result = integer_power(u, static_cast<int>(c)); break;
            case O::PowerConstantExponent: result = std::pow(u, c); break;
            case O::PowerConstantBase: result = std::pow(c, u); break;
            case O::Abs: result = std::fabs(u); break;
            case O::Sqrt: result = std::sqrt(u); break;
            case O::Exp: result = std::exp(u); break;
            case O::Log: result = std::log(u); break;
            case O::Log10: result = std::log10(u); break;
            case O::Sin: result = std::sin(u); break;
            case O::Cos: result = std::cos(u); break;
            case O::Tan: result = std::tan(u); break;
            case O::Sinh: result = std::sinh(u); break;
            case O::Cosh: result = std::cosh(u); break;
            case O::Tanh: result = std::tanh(u); break;
            case O::Asin: result = std::asin(u); break;
            case O::Acos: result = std::acos(u); break;
            case O::Atan: result = std::atan(u); break;
            case O::Asinh: result = std::asinh(u); break;
            case O::Acosh: result = std::acosh(u); break;
            case O::Atanh: result = std::atanh(u); break;
            case O::Logistic: result = 1. / (1. + std::exp(-u)); break;
            case O::SignPowerConstant: result = std::copysign(std::pow(std::fabs(u), c), u); break;
            case O::PiecewiseLinear: result = evaluate_piecewise_linear(tape.constants + instruction.second, u, nullptr); break;
            case O::Add: result = u + w; break;
            case O::Subtract: result = u - w; break;
            case O::Multiply: result = u * w; break;
            case O::Divide: result = u / w; break;
            case O::Power: result = std::pow(u, w); break;
            case O::Atan2: result = std::atan2(u, w); break;
            case O::Less: result = (u > w) ? u - w : 0.; break;
            case O::Remainder: result = std::fmod(u, w); break;
            case O::Sum: {
               const std::int32_t* operands = tape.arguments + instruction.first;
               result = 0.;
               for (std::int32_t k = 0; k < instruction.second; ++k) result += values[operands[k]];
               break;
            }
            case O::Minimum: case O::Maximum: case O::IfThenElse:
               result = values[selected_operand(instruction, tape, values)];
               break;
            case O::Floor: result = std::floor(u); break;
            case O::Ceil: result = std::ceil(u); break;
            case O::Round: result = round_to_decimals(u, w); break;
            case O::Trunc: result = truncate_to_decimals(u, w); break;
            case O::Precision: result = round_to_significant_digits(u, w); break;
            case O::IntegerDivide: result = std::trunc(u / w); break;
            case O::LessThan: result = (u < w) ? 1. : 0.; break;
            case O::LessEqual: result = (u <= w) ? 1. : 0.; break;
            case O::Equal: result = (u == w) ? 1. : 0.; break;
            case O::GreaterEqual: result = (u >= w) ? 1. : 0.; break;
            case O::GreaterThan: result = (u > w) ? 1. : 0.; break;
            case O::NotEqual: result = (u != w) ? 1. : 0.; break;
            case O::And: result = (u != 0. && w != 0.) ? 1. : 0.; break;
            case O::Or: result = (u != 0. || w != 0.) ? 1. : 0.; break;
            case O::Not: result = (u == 0.) ? 1. : 0.; break;
            case O::AndList: case O::OrList: {
               const std::int32_t* operands = tape.arguments + instruction.first;
               const bool is_and = (instruction.operation == O::AndList);
               bool outcome = is_and;
               for (std::int32_t k = 0; k < instruction.second; ++k) {
                  const bool operand_true = (values[operands[k]] != 0.);
                  outcome = is_and ? (outcome && operand_true) : (outcome || operand_true);
               }
               result = outcome ? 1. : 0.;
               break;
            }
            default: result = std::numeric_limits<double>::quiet_NaN();
         }
         return result;
   }
   } // namespace

   void evaluate_tape_values(const TapeView& tape, const double* x, double* values) {
      for (std::size_t i = tape.begin; i < tape.end; ++i) values[i] = evaluate_instruction(tape, tape.instructions[i], x, values);
   }

   void evaluate_tape_values_and_partials(const TapeView& tape, const double* x, double* values, double* partials,
         double* operand_partials, const double* defined_values) {
      using O = TapeOperation;
      for (std::size_t i = tape.begin; i < tape.end; ++i) {
         const TapeInstruction& instruction = tape.instructions[i];
         const double c = instruction.constant;
         const int arity = operand_arity(instruction.operation);
         const double u = (arity >= 1) ? values[instruction.first] : 0.;
         const double w = (arity == 2) ? values[instruction.second] : 0.;
         double result = 0., pu = 0., pw = 0.;
         switch (instruction.operation) {
            case O::Constant: result = c; break;
            case O::Variable: result = x[instruction.first] + c; break; // shifted variable leaf (offset c)
            case O::Negate: result = -u; pu = -1.; break;
            case O::AddConstant: result = u + c; pu = 1.; break;
            case O::MultiplyConstant: result = c * u; pu = c; break;
            case O::Square: result = u * u; pu = 2. * u; break;
            case O::PowerInteger: {
               const int n = static_cast<int>(c);
               const double power = integer_power(u, n - 1);
               result = power * u;
               pu = c * power;
               break;
            }
            case O::PowerConstantExponent: {
               const double power = std::pow(u, c - 1.);
               result = power * u;
               pu = c * power;
               if (u == 0.) result = std::pow(u, c);
               break;
            }
            case O::PowerConstantBase: result = std::pow(c, u); pu = result * std::log(c); break;
            case O::Abs: result = std::fabs(u); pu = (u < 0.) ? -1. : 1.; break;
            case O::Sqrt: result = std::sqrt(u); pu = 0.5 / result; break;
            case O::Exp: result = std::exp(u); pu = result; break;
            case O::Log: result = std::log(u); pu = 1. / u; break;
            case O::Log10: result = std::log10(u); pu = 1. / (u * 2.302585092994045684); break;
            case O::Sin: sine_and_cosine(u, result, pu); break;
            case O::Cos: { double sine; sine_and_cosine(u, sine, result); pu = -sine; break; }
            case O::Tan: result = std::tan(u); pu = 1. + result * result; break;
            case O::Sinh: result = std::sinh(u); pu = std::cosh(u); break;
            case O::Cosh: result = std::cosh(u); pu = std::sinh(u); break;
            case O::Tanh: result = std::tanh(u); pu = 1. - result * result; break;
            case O::Asin: result = std::asin(u); pu = 1. / std::sqrt(1. - u * u); break;
            case O::Acos: result = std::acos(u); pu = -1. / std::sqrt(1. - u * u); break;
            case O::Atan: result = std::atan(u); pu = 1. / (1. + u * u); break;
            case O::Asinh: result = std::asinh(u); pu = 1. / std::sqrt(1. + u * u); break;
            case O::Acosh: result = std::acosh(u); pu = 1. / std::sqrt(u * u - 1.); break;
            case O::Atanh: result = std::atanh(u); pu = 1. / (1. - u * u); break;
            case O::Logistic: result = 1. / (1. + std::exp(-u)); pu = result * (1. - result); break;
            case O::SignPowerConstant: {
               const double magnitude = std::fabs(u);
               result = std::copysign(std::pow(magnitude, c), u);
               pu = c * std::pow(magnitude, c - 1.);
               break;
            }
            case O::PiecewiseLinear: result = evaluate_piecewise_linear(tape.constants + instruction.second, u, &pu); break;
            case O::Add: result = u + w; pu = 1.; pw = 1.; break;
            case O::Subtract: result = u - w; pu = 1.; pw = -1.; break;
            case O::Multiply: result = u * w; pu = w; pw = u; break;
            case O::Divide: pu = 1. / w; result = u * pu; pw = -result * pu; break;
            case O::Power: {
               const double power = std::pow(u, w - 1.);
               result = std::pow(u, w);
               pu = w * power;
               pw = result * ((u > 0.) ? std::log(u) : 0.);
               break;
            }
            case O::Atan2: {
               const double r = u * u + w * w;
               result = std::atan2(u, w);
               pu = w / r;
               pw = -u / r;
               break;
            }
            case O::Less: if (u > w) { result = u - w; pu = 1.; pw = -1.; } break;
            case O::Remainder: result = std::fmod(u, w); pu = 1.; pw = -std::trunc(u / w); break;
            case O::DefinedValue: result = defined_values[instruction.first]; break;
            case O::Minimum: case O::Maximum: case O::IfThenElse: {
               const std::int32_t selected = selected_operand(instruction, tape, values);
               const std::int32_t* operands = tape.arguments + instruction.first;
               bool found = false; // one unit partial even if the selected instruction is repeated among the operands
               for (std::int32_t k = 0; k < instruction.second; ++k) {
                  const bool is_selected = !found && operands[k] == selected;
                  operand_partials[instruction.first + k] = is_selected ? 1. : 0.;
                  found = found || is_selected;
               }
               result = values[selected];
               break;
            }
            default: result = evaluate_instruction(tape, instruction, x, values); // no stored partials
         }
         values[i] = result;
         partials[2 * i] = pu;
         partials[2 * i + 1] = pw;
      }
   }


   namespace {
      /// first (and optionally second) partial derivatives of a unary or binary instruction whose value is f
      inline void compute_instruction_partials(const TapeView& tape, const TapeInstruction& instruction,
            DerivativeCategory category, const double* values, double f, InstructionPartials& p, bool second_order) {
         using O = TapeOperation;
         const double u = values[instruction.first];
         const double c = instruction.constant;
         if (category == DerivativeCategory::Unary) {
            double first = 0., first_first = 0.;
            switch (instruction.operation) {
               case O::Negate: first = -1.; break;
               case O::AddConstant: first = 1.; break;
               case O::MultiplyConstant: first = c; break;
               case O::Square: first = 2. * u; first_first = 2.; break;
               case O::PowerInteger: {
                  const int n = static_cast<int>(c);
                  const double power_n_2 = integer_power(u, n - 2);
                  first = c * power_n_2 * u;
                  first_first = c * (c - 1.) * power_n_2;
                  if (u == 0. && n == 1) first = 1.;
                  break;
               }
               case O::PowerConstantExponent:
                  first = c * std::pow(u, c - 1.);
                  if (second_order) first_first = c * (c - 1.) * std::pow(u, c - 2.);
                  break;
               case O::PowerConstantBase: {
                  const double log_base = std::log(c);
                  first = f * log_base;
                  first_first = first * log_base;
                  break;
               }
               case O::Abs: first = (u < 0.) ? -1. : 1.; break;
               case O::Sqrt: first = 0.5 / f; first_first = -0.25 / (f * f * f); break;
               case O::Exp: first = f; first_first = f; break;
               case O::Log: first = 1. / u; first_first = -first * first; break;
               case O::Log10: first = 1. / (u * 2.302585092994045684); first_first = -first / u; break;
               case O::Sin: first = std::cos(u); first_first = -f; break;
               case O::Cos: first = -std::sin(u); first_first = -f; break;
               case O::Tan: first = 1. + f * f; first_first = 2. * f * first; break;
               case O::Sinh: first = std::cosh(u); first_first = f; break;
               case O::Cosh: first = std::sinh(u); first_first = f; break;
               case O::Tanh: first = 1. - f * f; first_first = -2. * f * first; break;
               case O::Asin: case O::Acos: {
                  const double t = 1. / std::sqrt(1. - u * u);
                  const double sign = (instruction.operation == O::Asin) ? 1. : -1.;
                  first = sign * t;
                  first_first = sign * u * t * t * t;
                  break;
               }
               case O::Atan: { const double t = 1. / (1. + u * u); first = t; first_first = -2. * u * t * t; break; }
               case O::Asinh: { const double t = 1. / std::sqrt(1. + u * u); first = t; first_first = -u * t * t * t; break; }
               case O::Acosh: { const double t = 1. / std::sqrt(u * u - 1.); first = t; first_first = -u * t * t * t; break; }
               case O::Atanh: { const double t = 1. / (1. - u * u); first = t; first_first = 2. * u * t * t; break; }
               case O::Logistic: first = f * (1. - f); first_first = first * (1. - 2. * f); break;
               case O::SignPowerConstant: {
                  const double magnitude = std::fabs(u);
                  first = c * std::pow(magnitude, c - 1.);
                  if (second_order) first_first = std::copysign(c * (c - 1.) * std::pow(magnitude, c - 2.), u);
                  break;
               }
               case O::PiecewiseLinear: evaluate_piecewise_linear(tape.constants + instruction.second, u, &first); break;
               default: break;
            }
            p.first = first;
            p.first_first = first_first;
         }
         else {
            const double w = values[instruction.second];
            double first = 0., second = 0., first_first = 0., first_second = 0., second_second = 0.;
            switch (instruction.operation) {
               case O::Add: first = 1.; second = 1.; break;
               case O::Subtract: first = 1.; second = -1.; break;
               case O::Multiply: first = w; second = u; first_second = 1.; break;
               case O::Divide:
                  first = 1. / w;
                  second = -f / w;
                  first_second = -first * first;
                  second_second = 2. * f / (w * w);
                  break;
               case O::Power: {
                  const double log_u = (u > 0.) ? std::log(u) : 0.;
                  const double power_minus_one = std::pow(u, w - 1.);
                  first = w * power_minus_one;
                  second = f * log_u;
                  if (second_order) {
                     first_first = w * (w - 1.) * std::pow(u, w - 2.);
                     first_second = power_minus_one * (1. + w * log_u);
                     second_second = second * log_u;
                  }
                  break;
               }
               case O::Atan2: {
                  const double r = u * u + w * w, r2 = r * r;
                  first = w / r;
                  second = -u / r;
                  first_first = -2. * u * w / r2;
                  second_second = 2. * u * w / r2;
                  first_second = (u * u - w * w) / r2;
                  break;
               }
               case O::Less: if (u > w) { first = 1.; second = -1.; } break;
               case O::Remainder: first = 1.; second = -std::trunc(u / w); break;
               default: break;
            }
            p.first = first;
            p.second = second;
            p.first_first = first_first;
            p.first_second = first_second;
            p.second_second = second_second;
         }
      }
   } // namespace

   // Conventions: values and first partials are indexed like the instructions (global indices: operands refer to
   // the same indices); adjoints and the scratch arrays of the second-order sweeps (partials, tangents,
   // second-order adjoints) are indexed from tape.begin.

   void evaluate_tape_partials(const TapeView& tape, const double* values, InstructionPartials* partials,
         bool second_order) {
      for (std::size_t i = tape.begin; i < tape.end; ++i) {
         const TapeInstruction& instruction = tape.instructions[i];
         const DerivativeCategory category = derivative_category(instruction.operation);
         if (category != DerivativeCategory::Unary && category != DerivativeCategory::Binary) continue;
         compute_instruction_partials(tape, instruction, category, values, values[i], partials[i - tape.begin], second_order);
      }
   }

   namespace {
      /// first-order reverse sweep: on entry, adjoints (indexed from tape.begin) hold the seeds
      void propagate_first_order(const TapeView& tape, const double* values, const InstructionPartials* partials, double* adjoints) {
         const std::size_t b = tape.begin, length = tape.end - tape.begin;
         for (std::size_t r = length; r-- > 0;) {
            const double adjoint = adjoints[r];
            if (adjoint == 0.) continue;
            const TapeInstruction& instruction = tape.instructions[b + r];
            switch (derivative_category(instruction.operation)) {
               case DerivativeCategory::Unary:
                  adjoints[instruction.first - b] += adjoint * partials[r].first;
                  break;
               case DerivativeCategory::Binary:
                  adjoints[instruction.first - b] += adjoint * partials[r].first;
                  adjoints[instruction.second - b] += adjoint * partials[r].second;
                  break;
               case DerivativeCategory::Sum: {
                  const std::int32_t* operands = tape.arguments + instruction.first;
                  for (std::int32_t k = 0; k < instruction.second; ++k) adjoints[operands[k] - b] += adjoint;
                  break;
               }
               case DerivativeCategory::Select:
                  adjoints[selected_operand(instruction, tape, values) - b] += adjoint;
                  break;
               default: break;
            }
         }
      }

      /// Forward tangents (leaf tangents given by leaf_tangent(instruction)) and second-order adjoints; on entry,
      /// second_order_adjoints hold their seeds. Scratch arrays are indexed from tape.begin.
      template <typename LeafTangent>
      void tangent_and_second_order_sweeps(const TapeView& tape, const double* values, const InstructionPartials* partials,
            const double* adjoints, double* tangents, double* second_order_adjoints, LeafTangent leaf_tangent) {
         const std::size_t b = tape.begin, length = tape.end - tape.begin;
         for (std::size_t r = 0; r < length; ++r) {
            const TapeInstruction& instruction = tape.instructions[b + r];
            double t = 0.;
            switch (derivative_category(instruction.operation)) {
               case DerivativeCategory::Leaf: t = leaf_tangent(instruction); break;
               case DerivativeCategory::Unary: t = partials[r].first * tangents[instruction.first - b]; break;
               case DerivativeCategory::Binary:
                  t = partials[r].first * tangents[instruction.first - b] + partials[r].second * tangents[instruction.second - b];
                  break;
               case DerivativeCategory::Sum: {
                  const std::int32_t* operands = tape.arguments + instruction.first;
                  for (std::int32_t k = 0; k < instruction.second; ++k) t += tangents[operands[k] - b];
                  break;
               }
               case DerivativeCategory::Select: t = tangents[selected_operand(instruction, tape, values) - b]; break;
               default: break;
            }
            tangents[r] = t;
         }
         for (std::size_t r = length; r-- > 0;) {
            const double adjoint = adjoints[r], s = second_order_adjoints[r];
            if (adjoint == 0. && s == 0.) continue;
            const TapeInstruction& instruction = tape.instructions[b + r];
            const InstructionPartials& p = partials[r];
            switch (derivative_category(instruction.operation)) {
               case DerivativeCategory::Unary:
                  second_order_adjoints[instruction.first - b] += s * p.first + adjoint * p.first_first * tangents[instruction.first - b];
                  break;
               case DerivativeCategory::Binary: {
                  const double tu = tangents[instruction.first - b], tw = tangents[instruction.second - b];
                  second_order_adjoints[instruction.first - b] += s * p.first + adjoint * (p.first_first * tu + p.first_second * tw);
                  second_order_adjoints[instruction.second - b] += s * p.second + adjoint * (p.first_second * tu + p.second_second * tw);
                  break;
               }
               case DerivativeCategory::Sum: {
                  const std::int32_t* operands = tape.arguments + instruction.first;
                  for (std::int32_t k = 0; k < instruction.second; ++k) second_order_adjoints[operands[k] - b] += s;
                  break;
               }
               case DerivativeCategory::Select: second_order_adjoints[selected_operand(instruction, tape, values) - b] += s; break;
               default: break;
            }
         }
      }
   } // namespace

   void evaluate_tape_adjoints(const TapeView& tape, const double* values, const InstructionPartials* partials,
         double* adjoints) {
      const std::size_t length = tape.end - tape.begin;
      std::fill(adjoints, adjoints + length, 0.);
      adjoints[length - 1] = 1.;
      propagate_first_order(tape, values, partials, adjoints);
   }

   namespace {
      /// Forward-over-reverse in vector mode: the K unit directions of the local variables
      /// first_direction, ..., first_direction + number_directions - 1 in one pair of sweeps
      /// (tangents and second_order_adjoints: K entries per instruction, indexed from tape.begin).
      template <int K>
      void hessian_directions(const TapeView& tape, const double* values, const InstructionPartials* partials,
            const double* adjoints, double* tangents, double* second_order_adjoints, const std::int32_t* direction_of_local,
            std::size_t first_direction, std::size_t number_directions, const DefinedHessianLinks* links) {
         const std::size_t b = tape.begin, length = tape.end - tape.begin;
         for (std::size_t r = 0; r < length; ++r) {
            const TapeInstruction& instruction = tape.instructions[b + r];
            double* t = tangents + K * r;
            switch (derivative_category(instruction.operation)) {
               case DerivativeCategory::Leaf: {
                  for (int d = 0; d < K; ++d) t[d] = 0.;
                  if (instruction.operation == TapeOperation::Variable) {
                     const std::int32_t direction = (direction_of_local != nullptr) ? direction_of_local[instruction.second] : instruction.second;
                     const auto d = static_cast<std::size_t>(direction);
                     if (direction >= 0 && d >= first_direction && d < first_direction + number_directions) t[d - first_direction] = 1.;
                  }
                  else if (instruction.operation == TapeOperation::DefinedValue) {
                     // tangent of a shared defined variable along e_d: its partial derivative w.r.t. the local variable d
                     const DefinedLeafLink& link = links->leaves[static_cast<std::size_t>(instruction.constant)];
                     for (std::size_t p = link.pair_begin; p < link.pair_begin + link.pair_count; ++p) {
                        const std::size_t d = links->pair_locals[p];
                        if (d >= first_direction && d < first_direction + number_directions) {
                           t[d - first_direction] = links->gradients[links->pair_gradients[p]];
                        }
                     }
                  }
                  break;
               }
               case DerivativeCategory::Unary: {
                  const double* tu = tangents + K * (instruction.first - b);
                  const double pu = partials[r].first;
                  for (int d = 0; d < K; ++d) t[d] = pu * tu[d];
                  break;
               }
               case DerivativeCategory::Binary: {
                  const double* tu = tangents + K * (instruction.first - b);
                  const double* tw = tangents + K * (instruction.second - b);
                  const double pu = partials[r].first, pw = partials[r].second;
                  for (int d = 0; d < K; ++d) t[d] = pu * tu[d] + pw * tw[d];
                  break;
               }
               case DerivativeCategory::Sum: {
                  const std::int32_t* operands = tape.arguments + instruction.first;
                  for (int d = 0; d < K; ++d) t[d] = 0.;
                  for (std::int32_t k = 0; k < instruction.second; ++k) {
                     const double* to = tangents + K * (operands[k] - b);
                     for (int d = 0; d < K; ++d) t[d] += to[d];
                  }
                  break;
               }
               case DerivativeCategory::Select: {
                  const double* ts = tangents + K * (selected_operand(instruction, tape, values) - b);
                  for (int d = 0; d < K; ++d) t[d] = ts[d];
                  break;
               }
               default:
                  for (int d = 0; d < K; ++d) t[d] = 0.;
            }
         }
         std::fill(second_order_adjoints, second_order_adjoints + K * length, 0.);
         for (std::size_t r = length; r-- > 0;) {
            const TapeInstruction& instruction = tape.instructions[b + r];
            const double* s = second_order_adjoints + K * r;
            const double adjoint = adjoints[r];
            const InstructionPartials& p = partials[r];
            switch (derivative_category(instruction.operation)) {
               case DerivativeCategory::Unary: {
                  double* su = second_order_adjoints + K * (instruction.first - b);
                  const double* tu = tangents + K * (instruction.first - b);
                  const double curvature = adjoint * p.first_first;
                  for (int d = 0; d < K; ++d) su[d] += s[d] * p.first + curvature * tu[d];
                  break;
               }
               case DerivativeCategory::Binary: {
                  double* su = second_order_adjoints + K * (instruction.first - b);
                  double* sw = second_order_adjoints + K * (instruction.second - b);
                  const double* tu = tangents + K * (instruction.first - b);
                  const double* tw = tangents + K * (instruction.second - b);
                  const double auu = adjoint * p.first_first, auw = adjoint * p.first_second, aww = adjoint * p.second_second;
                  for (int d = 0; d < K; ++d) {
                     const double tud = tu[d], twd = tw[d], sd = s[d];
                     su[d] += sd * p.first + auu * tud + auw * twd;
                     sw[d] += sd * p.second + auw * tud + aww * twd;
                  }
                  break;
               }
               case DerivativeCategory::Sum: {
                  const std::int32_t* operands = tape.arguments + instruction.first;
                  for (std::int32_t k = 0; k < instruction.second; ++k) {
                     double* so = second_order_adjoints + K * (operands[k] - b);
                     for (int d = 0; d < K; ++d) so[d] += s[d];
                  }
                  break;
               }
               case DerivativeCategory::Select: {
                  double* ss = second_order_adjoints + K * (selected_operand(instruction, tape, values) - b);
                  for (int d = 0; d < K; ++d) ss[d] += s[d];
                  break;
               }
               default: break;
            }
         }
      }
   } // namespace

   namespace {
      /// packed_lower += factor * Hessian of the range, for the directions given by the local variables
      /// local_of_direction[0 .. number_directions) (partials and adjoints of the range already computed)
      void accumulate_hessian(const TapeView& tape, const double* values, const InstructionPartials* partials,
            const double* adjoints, double* tangents, double* second_order_adjoints, const std::int32_t* direction_of_local,
            const std::int32_t* local_of_direction, std::size_t number_directions, std::size_t number_variables,
            double factor, double* packed_lower, const DefinedHessianLinks* links = nullptr) {
         for (std::size_t first = 0; first < number_directions; first += maximum_hessian_directions) {
            const std::size_t count = std::min(maximum_hessian_directions, number_directions - first);
            std::size_t width = 1;
            if (count > 4) {
               width = 8;
               hessian_directions<8>(tape, values, partials, adjoints, tangents, second_order_adjoints, direction_of_local, first, count, links);
            }
            else if (count > 2) {
               width = 4;
               hessian_directions<4>(tape, values, partials, adjoints, tangents, second_order_adjoints, direction_of_local, first, count, links);
            }
            else if (count == 2) {
               width = 2;
               hessian_directions<2>(tape, values, partials, adjoints, tangents, second_order_adjoints, direction_of_local, first, count, links);
            }
            else hessian_directions<1>(tape, values, partials, adjoints, tangents, second_order_adjoints, direction_of_local, first, count, links);
            // every occurrence of a variable (row a) contributes to the columns b = local_of_direction[first + d], a >= b
            for (std::size_t i = tape.begin; i < tape.end; ++i) {
               const TapeInstruction& instruction = tape.instructions[i];
               if (instruction.operation == TapeOperation::DefinedValue) {
                  // chain rule through a shared defined variable v: s_v(d) * dv/dx_a; its own curvature
                  // (adjoint * Hessian of v) is accumulated once for all functions in the weights
                  const DefinedLeafLink& link = links->leaves[static_cast<std::size_t>(instruction.constant)];
                  const double* s = second_order_adjoints + width * (i - tape.begin);
                  for (std::size_t p = link.pair_begin; p < link.pair_begin + link.pair_count; ++p) {
                     const std::size_t a = links->pair_locals[p];
                     const double gradient = factor * links->gradients[links->pair_gradients[p]];
                     for (std::size_t d = 0; d < count; ++d) {
                        const std::size_t b = first + d;
                        if (a >= b && s[d] != 0.) packed_lower[packed_lower_index(a, b, number_variables)] += s[d] * gradient;
                     }
                  }
                  if (first == 0) links->weights[instruction.first] += factor * adjoints[i - tape.begin];
                  continue;
               }
               if (instruction.operation != TapeOperation::Variable) continue;
               const auto a = static_cast<std::size_t>(instruction.second);
               const double* s = second_order_adjoints + width * (i - tape.begin);
               for (std::size_t d = 0; d < count; ++d) {
                  const auto b = static_cast<std::size_t>((local_of_direction != nullptr) ? local_of_direction[first + d]
                     : static_cast<std::int32_t>(first + d));
                  if (a >= b && s[d] != 0.) packed_lower[packed_lower_index(a, b, number_variables)] += factor * s[d];
               }
            }
         }
      }
   } // namespace

   void evaluate_tape_hessian(const TapeView& tape, const double* values, InstructionPartials* partials,
         double* adjoints, double* tangents, double* second_order_adjoints, std::size_t number_variables, double factor,
         double* packed_lower, const DefinedHessianLinks* links) {
      evaluate_tape_partials(tape, values, partials, true);
      evaluate_tape_adjoints(tape, values, partials, adjoints);
      accumulate_hessian(tape, values, partials, adjoints, tangents, second_order_adjoints, nullptr, nullptr,
         number_variables, number_variables, factor, packed_lower, links);
   }

   bool evaluate_group_hessian(const TapeView& tape, std::size_t group_sum,
         const std::pair<std::size_t, std::size_t>* operand_ranges, std::size_t number_ranges, const double* values,
         InstructionPartials* partials, double* adjoints, double* tangents, double* second_order_adjoints,
         std::int32_t* direction_of_local, std::int32_t* local_of_direction, double* sum_gradient,
         std::size_t number_variables, double factor, double* packed_lower) {
      // phi: the chain of unary instructions from the sum to the output
      std::size_t chain[64];
      std::size_t chain_length = 0;
      for (std::size_t i = tape.end - 1; i != group_sum;) {
         const TapeInstruction& instruction = tape.instructions[i];
         if (derivative_category(instruction.operation) != DerivativeCategory::Unary || chain_length == 64) return false;
         chain[chain_length++] = i;
         i = static_cast<std::size_t>(instruction.first);
         if (i < group_sum) return false;
      }
      double first_derivative = 1., second_derivative = 0.; // of phi at the sum
      for (std::size_t k = chain_length; k-- > 0;) {
         const TapeView single{tape.instructions, chain[k], chain[k] + 1, tape.arguments, tape.constants};
         evaluate_tape_partials(single, values, partials, true);
         second_derivative = partials[0].first_first * first_derivative * first_derivative + partials[0].first * second_derivative;
         first_derivative *= partials[0].first;
      }
      // H = phi'' g g^T + phi' sum_k H_k, with g the gradient of the sum and H_k the Hessians of its operands
      std::fill(sum_gradient, sum_gradient + number_variables, 0.);
      for (std::size_t r = 0; r < number_ranges; ++r) {
         const TapeView range{tape.instructions, operand_ranges[r].first, operand_ranges[r].second, tape.arguments, tape.constants};
         evaluate_tape_partials(range, values, partials, true);
         evaluate_tape_adjoints(range, values, partials, adjoints);
         std::size_t number_directions = 0;
         for (std::size_t i = range.begin; i < range.end; ++i) {
            const TapeInstruction& instruction = range.instructions[i];
            if (instruction.operation != TapeOperation::Variable) continue;
            sum_gradient[instruction.second] += adjoints[i - range.begin];
            if (direction_of_local[instruction.second] < 0) {
               direction_of_local[instruction.second] = static_cast<std::int32_t>(number_directions);
               local_of_direction[number_directions++] = instruction.second;
            }
         }
         if (first_derivative != 0. && number_directions > 0) {
            accumulate_hessian(range, values, partials, adjoints, tangents, second_order_adjoints, direction_of_local,
               local_of_direction, number_directions, number_variables, factor * first_derivative, packed_lower);
         }
         for (std::size_t d = 0; d < number_directions; ++d) direction_of_local[local_of_direction[d]] = -1;
      }
      const double curvature = factor * second_derivative;
      if (curvature != 0.) {
         double* column = packed_lower;
         for (std::size_t b = 0; b < number_variables; ++b) { // rank-1 update, column by column
            const double gb = curvature * sum_gradient[b];
            const std::size_t height = number_variables - b;
            const double* g = sum_gradient + b;
            if (gb != 0.) for (std::size_t r = 0; r < height; ++r) column[r] += gb * g[r];
            column += height;
         }
      }
      return true;
   }

   void evaluate_tape_hessian_vector_product(const TapeView& tape, const double* values, InstructionPartials* partials,
         double* adjoints, double* tangents, double* second_order_adjoints, const double* local_direction, double factor,
         double* result, const DefinedVariableLinks* links) {
      evaluate_tape_partials(tape, values, partials, true);
      evaluate_tape_adjoints(tape, values, partials, adjoints);
      const std::size_t b = tape.begin, length = tape.end - tape.begin;
      std::fill(second_order_adjoints, second_order_adjoints + length, 0.);
      tangent_and_second_order_sweeps(tape, values, partials, adjoints, tangents, second_order_adjoints,
         [&](const TapeInstruction& instruction) {
            if (instruction.operation == TapeOperation::Variable) return local_direction[instruction.second];
            if (instruction.operation == TapeOperation::DefinedValue) return links->tangents[instruction.first];
            return 0.;
         });
      for (std::size_t r = 0; r < length; ++r) {
         const TapeInstruction& instruction = tape.instructions[b + r];
         if (instruction.operation == TapeOperation::Variable) result[instruction.first] += factor * second_order_adjoints[r];
         else if (instruction.operation == TapeOperation::DefinedValue) {
            links->adjoint_seeds[instruction.first] += factor * adjoints[r];
            links->second_order_seeds[instruction.first] += factor * second_order_adjoints[r];
         }
      }
   }

   void evaluate_seeded_hessian_vector_product(const TapeView& tape, const double* values, InstructionPartials* partials,
         double* adjoints, double* tangents, double* second_order_adjoints, const double* direction, double* result) {
      evaluate_tape_partials(tape, values, partials, true);
      propagate_first_order(tape, values, partials, adjoints);
      tangent_and_second_order_sweeps(tape, values, partials, adjoints, tangents, second_order_adjoints,
         [&](const TapeInstruction& instruction) {
            return (instruction.operation == TapeOperation::Variable) ? direction[instruction.first] : 0.;
         });
      for (std::size_t r = 0; r < tape.end - tape.begin; ++r) {
         const TapeInstruction& instruction = tape.instructions[tape.begin + r];
         if (instruction.operation == TapeOperation::Variable) result[instruction.first] += second_order_adjoints[r];
      }
   }

} // namespace cppasl::detail
