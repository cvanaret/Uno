// Reader of AMPL .nl files (text 'g' and binary 'b'/'z'/'h' formats).
// Reference: D. M. Gay, "Writing .nl Files" (2005), and ASL's solvers2/pfghread.c.
#include <charconv>
#include <cstdio>
#include <cstring>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>
#include "cppasl/nl_model.hpp"
#include "model_builder.hpp"

#if defined(__unix__) || defined(__APPLE__)
#include <fcntl.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <unistd.h>
#define CPPASL_HAS_MMAP 1
#endif

namespace cppasl {

   namespace detail {
      namespace {
         constexpr double infinity = std::numeric_limits<double>::infinity();

         [[noreturn]] void parse_error(std::size_t offset, const std::string& message) {
            throw std::runtime_error("cppasl: .nl parse error at byte " + std::to_string(offset) + ": " + message);
         }

         // -------------------------------------------------------------------------------------------------------------
         // header
         // -------------------------------------------------------------------------------------------------------------

         /// Reads the integers at the beginning of a header line (stops at the first non-integer token).
         std::size_t read_header_integers(const char*& cursor, const char* end, long long* values, std::size_t capacity) {
            const char* line_end = static_cast<const char*>(std::memchr(cursor, '\n', static_cast<std::size_t>(end - cursor)));
            if (line_end == nullptr) line_end = end;
            std::size_t count = 0;
            const char* p = cursor;
            while (count < capacity) {
               while (p < line_end && (*p == ' ' || *p == '\t')) ++p;
               long long value = 0;
               const auto [next, error] = std::from_chars(p, line_end, value);
               if (error != std::errc()) break;
               values[count++] = value;
               p = next;
            }
            cursor = (line_end < end) ? line_end + 1 : end;
            return count;
         }

         NlHeader parse_header(const char*& cursor, const char* end) {
            NlHeader header;
            if (cursor >= end) parse_error(0, "empty file");
            switch (*cursor) {
               case 'g': case 'G': header.format = NlFormat::Text; break;
               case 'b': case 'B': header.format = NlFormat::Binary; break;
               case 'z': case 'Z': header.format = NlFormat::BinaryShortOpcodes; break;
               case 'h': case 'H': header.format = NlFormat::BinaryLongIntegers; break;
               default: parse_error(0, "unknown format letter");
            }
            const char* line_begin = cursor;
            ++cursor;
            long long values[16]{};
            // line 1: number of options, the options, and possibly the variable-bound tolerance
            {
               const char* line_end = static_cast<const char*>(std::memchr(cursor, '\n', static_cast<std::size_t>(end - cursor)));
               if (line_end == nullptr) parse_error(0, "truncated header");
               const char* p = cursor;
               auto next_integer = [&](long long& value) {
                  while (p < line_end && (*p == ' ' || *p == '\t')) ++p;
                  const auto [next, error] = std::from_chars(p, line_end, value);
                  if (error != std::errc()) return false;
                  p = next;
                  return true;
               };
               long long number_options = 0;
               if (next_integer(number_options)) {
                  header.ampl_options[0] = static_cast<int>(std::min<long long>(number_options, 9));
                  for (int k = 1; k <= header.ampl_options[0]; ++k) {
                     long long option = 0;
                     if (!next_integer(option)) break;
                     header.ampl_options[k] = static_cast<int>(option);
                  }
                  if (header.ampl_options[2] == 3) {
                     while (p < line_end && (*p == ' ' || *p == '\t')) ++p;
                     std::from_chars(p, line_end, header.variable_bound_tolerance);
                  }
               }
               cursor = line_end + 1;
            }
            auto line = [&](std::size_t minimum, std::size_t capacity) {
               const char* start = cursor;
               std::fill(values, values + 16, 0);
               if (read_header_integers(cursor, end, values, capacity) < minimum) {
                  parse_error(static_cast<std::size_t>(start - line_begin), "malformed header line");
               }
            };
            line(3, 6);
            header.number_variables = static_cast<int>(values[0]);
            header.number_constraints = static_cast<int>(values[1]);
            header.number_objectives = static_cast<int>(values[2]);
            header.number_ranges = static_cast<int>(values[3]);
            header.number_equations = static_cast<int>(values[4]);
            header.number_logical_constraints = static_cast<int>(values[5]);
            line(2, 6);
            header.number_nonlinear_constraints = static_cast<int>(values[0]);
            header.number_nonlinear_objectives = static_cast<int>(values[1]);
            header.number_linear_complementarities = static_cast<int>(values[2]);
            header.number_nonlinear_complementarities = static_cast<int>(values[3]);
            header.number_double_inequality_complementarities = static_cast<int>(values[4]);
            header.number_nonzero_lower_bound_complementarities = static_cast<int>(values[5]);
            line(2, 2);
            header.number_nonlinear_network_constraints = static_cast<int>(values[0]);
            header.number_linear_network_constraints = static_cast<int>(values[1]);
            line(2, 3);
            header.number_nonlinear_variables_in_constraints = static_cast<int>(values[0]);
            header.number_nonlinear_variables_in_objectives = static_cast<int>(values[1]);
            header.number_nonlinear_variables_in_both = static_cast<int>(values[2]);
            line(2, 4);
            header.number_linear_network_variables = static_cast<int>(values[0]);
            header.number_imported_functions = static_cast<int>(values[1]);
            header.arithmetic_kind = static_cast<int>(values[2]);
            header.flags = static_cast<int>(values[3]);
            line(5, 5);
            header.number_binary_variables = static_cast<int>(values[0]);
            header.number_integer_variables = static_cast<int>(values[1]);
            header.number_nonlinear_integer_variables_in_both = static_cast<int>(values[2]);
            header.number_nonlinear_integer_variables_in_constraints = static_cast<int>(values[3]);
            header.number_nonlinear_integer_variables_in_objectives = static_cast<int>(values[4]);
            line(2, 2);
            header.number_jacobian_nonzeros = static_cast<std::size_t>(values[0]);
            header.number_objective_gradient_nonzeros = static_cast<std::size_t>(values[1]);
            line(2, 2);
            header.maximum_constraint_name_length = static_cast<int>(values[0]);
            header.maximum_variable_name_length = static_cast<int>(values[1]);
            line(5, 5);
            header.number_common_expressions_in_both = static_cast<int>(values[0]);
            header.number_common_expressions_in_constraints = static_cast<int>(values[1]);
            header.number_common_expressions_in_objectives = static_cast<int>(values[2]);
            header.number_common_expressions_in_single_constraint = static_cast<int>(values[3]);
            header.number_common_expressions_in_single_objective = static_cast<int>(values[4]);
            if (header.number_variables <= 0 || header.number_constraints < 0 || header.number_objectives < 0) {
               parse_error(0, "invalid problem dimensions");
            }
            return header;
         }

         // -------------------------------------------------------------------------------------------------------------
         // tokenizers: one "item" = a letter followed by values (in text, one line; the rest of the line is ignored)
         // -------------------------------------------------------------------------------------------------------------

         class TextTokenizer {
         public:
            TextTokenizer(const char* begin, const char* cursor, const char* end): begin(begin), cursor(cursor), end(end) {}

            /// next item letter, or -1 at the end of the file
            int read_letter() {
               while (this->cursor < this->end && is_space(*this->cursor)) ++this->cursor;
               if (this->cursor >= this->end) return -1;
               return static_cast<unsigned char>(*this->cursor++);
            }
            int read_int() { return static_cast<int>(this->read_integer<long long>()); }
            long long read_long() { return this->read_integer<long long>(); }
            int read_opcode() { return this->read_int(); }
            double read_short_constant() { return static_cast<double>(this->read_long()); }
            double read_long_constant() { return static_cast<double>(this->read_long()); }
            double read_double() {
               this->skip_blanks();
               double value = 0.;
               const auto [next, error] = std::from_chars(this->cursor, this->end, value);
               if (error != std::errc()) this->fail("expected a real number");
               this->cursor = next;
               return value;
            }
            std::string read_name() {
               this->skip_blanks();
               const char* start = this->cursor;
               while (this->cursor < this->end && !is_space(*this->cursor)) ++this->cursor;
               return std::string(start, this->cursor);
            }
            /// ignores the rest of the current line
            void end_item() {
               const void* newline = std::memchr(this->cursor, '\n', static_cast<std::size_t>(this->end - this->cursor));
               this->cursor = (newline != nullptr) ? static_cast<const char*>(newline) + 1 : this->end;
            }
            [[noreturn]] void fail(const std::string& message) const {
               parse_error(static_cast<std::size_t>(this->cursor - this->begin), message);
            }

         private:
            const char* begin;
            const char* cursor;
            const char* end;

            static bool is_space(char c) { return c == ' ' || c == '\t' || c == '\n' || c == '\r'; }
            void skip_blanks() {
               while (this->cursor < this->end && (*this->cursor == ' ' || *this->cursor == '\t')) ++this->cursor;
            }
            /// hand-written decimal parser (the hot path of text files: opcodes, variable indices and counts)
            template <typename Integer>
            Integer read_integer() {
               this->skip_blanks();
               const char* p = this->cursor;
               const bool negative = (p < this->end && *p == '-');
               if (negative) ++p;
               if (p >= this->end || static_cast<unsigned>(*p - '0') > 9u) this->fail("expected an integer");
               Integer value = 0;
               while (p < this->end && static_cast<unsigned>(*p - '0') <= 9u) value = 10 * value + (*p++ - '0');
               this->cursor = p;
               return negative ? -value : value;
            }
         };

         /// Native-endian binary encoding: letters are single bytes, integers IntegerType, opcodes OpcodeType.
         template <typename IntegerType, typename OpcodeType>
         class BinaryTokenizer {
         public:
            BinaryTokenizer(const char* begin, const char* cursor, const char* end): begin(begin), cursor(cursor), end(end) {}

            int read_letter() {
               if (this->cursor >= this->end) return -1;
               return static_cast<unsigned char>(*this->cursor++);
            }
            int read_int() { return static_cast<int>(this->read_raw<IntegerType>()); }
            long long read_long() { return static_cast<long long>(this->read_raw<IntegerType>()); }
            int read_opcode() { return static_cast<int>(this->read_raw<OpcodeType>()); }
            double read_short_constant() { return static_cast<double>(this->read_raw<std::int16_t>()); }
            double read_long_constant() { return static_cast<double>(this->read_raw<IntegerType>()); }
            double read_double() { return this->read_raw<double>(); }
            std::string read_name() {
               const auto length = static_cast<std::size_t>(this->read_raw<std::int32_t>());
               if (length > static_cast<std::size_t>(this->end - this->cursor)) this->fail("truncated string");
               std::string name(this->cursor, length);
               this->cursor += length;
               return name;
            }
            void end_item() {}
            [[noreturn]] void fail(const std::string& message) const {
               parse_error(static_cast<std::size_t>(this->cursor - this->begin), message);
            }

         private:
            const char* begin;
            const char* cursor;
            const char* end;

            template <typename T>
            T read_raw() {
               if (static_cast<std::size_t>(this->end - this->cursor) < sizeof(T)) this->fail("unexpected end of file");
               T value;
               std::memcpy(&value, this->cursor, sizeof(T));
               this->cursor += sizeof(T);
               return value;
            }
         };

         // -------------------------------------------------------------------------------------------------------------
         // segments
         // -------------------------------------------------------------------------------------------------------------

         template <typename Tokenizer>
         class SegmentReader {
         public:
            SegmentReader(Tokenizer tokenizer, const NlHeader& header, ModelBuilder& builder):
                  tokenizer(tokenizer), header(header), builder(builder), model(builder.model()),
                  number_variables(header.number_variables),
                  number_expression_variables(header.number_variables + header.number_defined_variables()) {
            }

            void read_all_segments() {
               for (int letter = this->tokenizer.read_letter(); letter >= 0; letter = this->tokenizer.read_letter()) {
                  switch (letter) {
                     case 'C': this->read_constraint(); break;
                     case 'O': this->read_objective(); break;
                     case 'V': this->read_defined_variable(); break;
                     case 'J': this->read_jacobian_row(); break;
                     case 'G': this->read_objective_gradient(); break;
                     case 'k': case 'K': this->read_column_counts(); break;
                     case 'r': this->read_constraint_bounds(); break;
                     case 'b': this->read_variable_bounds(); break;
                     case 'x': this->read_initial_point(this->model.initial_primal_point, this->number_variables); break;
                     case 'd': this->read_initial_point(this->model.initial_dual_point, this->model.number_constraints); break;
                     case 'S': this->read_suffix(); break;
                     case 'F': this->tokenizer.fail("imported functions (F segments) are not supported");
                     case 'L': this->tokenizer.fail("logical constraints (L segments) are not supported");
                     default: this->tokenizer.fail(std::string("unknown segment '") + static_cast<char>(letter) + "'");
                  }
               }
            }

         private:
            Tokenizer tokenizer;
            const NlHeader& header;
            ModelBuilder& builder;
            ModelData& model;
            const int number_variables;
            const int number_expression_variables; ///< variables + defined variables
            bool column_counts_seen{false};
            struct Frame {
               int opcode;
               std::uint32_t operand_count;
               std::size_t first_operand; ///< in operand_stack
            };
            std::vector<Frame> frames;
            std::vector<NodeIndex> operand_stack;
            std::vector<LinearTerm> linear_scratch;

            int read_index(int upper_bound, const char* what) {
               const int index = this->tokenizer.read_int();
               if (index < 0 || index >= upper_bound) this->tokenizer.fail(std::string("invalid ") + what + " index");
               return index;
            }

            double read_number_item() {
               switch (this->tokenizer.read_letter()) {
                  case 'n': return this->tokenizer.read_double();
                  case 's': return this->tokenizer.read_short_constant();
                  case 'l': return this->tokenizer.read_long_constant();
                  default: this->tokenizer.fail("expected a constant");
               }
            }

            /// Iterative (explicit stack) prefix-notation expression parser: deep expressions cannot overflow the stack.
            NodeIndex read_expression() {
               ExpressionArena& arena = this->builder.expression_arena();
               this->frames.clear();
               this->operand_stack.clear();
               for (;;) {
                  NodeIndex node = 0;
                  const int letter = this->tokenizer.read_letter();
                  switch (letter) {
                     case 'n': node = arena.add_constant(this->tokenizer.read_double()); break;
                     case 's': node = arena.add_constant(this->tokenizer.read_short_constant()); break;
                     case 'l': node = arena.add_constant(this->tokenizer.read_long_constant()); break;
                     case 'v': node = arena.add_variable(this->read_index(this->number_expression_variables, "variable")); break;
                     case 'o': {
                        const int opcode = this->tokenizer.read_opcode();
                        this->tokenizer.end_item();
                        std::uint32_t count = 0;
                        switch (operand_layout(opcode)) {
                           case OperandLayout::Unary: count = 1; break;
                           case OperandLayout::Binary: count = 2; break;
                           case OperandLayout::Ternary: count = 3; break;
                           case OperandLayout::CountedList: case OperandLayout::SumList: {
                              const int list_size = this->tokenizer.read_int();
                              this->tokenizer.end_item();
                              if (list_size < 1) this->tokenizer.fail("invalid operand count");
                              count = static_cast<std::uint32_t>(list_size);
                              break;
                           }
                           case OperandLayout::PiecewiseLinear: {
                              const int number_slopes = this->tokenizer.read_int();
                              this->tokenizer.end_item();
                              if (number_slopes < 2) this->tokenizer.fail("invalid piecewise-linear term");
                              this->frames.push_back({opcode, static_cast<std::uint32_t>(2 * number_slopes), this->operand_stack.size()});
                              for (int k = 0; k < 2 * number_slopes - 1; ++k) {
                                 this->operand_stack.push_back(arena.add_constant(this->read_number_item()));
                                 this->tokenizer.end_item();
                              }
                              continue;
                           }
                           default:
                              this->tokenizer.fail("unknown operation o" + std::to_string(opcode));
                        }
                        this->frames.push_back({opcode, count, this->operand_stack.size()});
                        continue;
                     }
                     case 'f': this->tokenizer.fail("imported function calls are not supported");
                     case 'h': this->tokenizer.fail("string arguments are not supported");
                     default: this->tokenizer.fail("invalid expression token");
                  }
                  this->tokenizer.end_item();
                  // a leaf was read: complete every operation whose operands are all available
                  this->operand_stack.push_back(node);
                  while (!this->frames.empty() &&
                        this->operand_stack.size() - this->frames.back().first_operand == this->frames.back().operand_count) {
                     const Frame frame = this->frames.back();
                     this->frames.pop_back();
                     node = arena.add_operation(frame.opcode, this->operand_stack.data() + frame.first_operand, frame.operand_count);
                     this->operand_stack.resize(frame.first_operand);
                     this->operand_stack.push_back(node);
                  }
                  if (this->frames.empty()) return this->operand_stack.back();
               }
            }

            void read_constraint() {
               const int constraint = this->read_index(this->model.number_constraints, "constraint");
               this->tokenizer.end_item();
               const ExpressionArena::Mark mark = this->builder.expression_arena().mark();
               const NodeIndex root = this->read_expression();
               this->builder.set_constraint_body(constraint, root, mark);
            }

            void read_objective() {
               const int objective = this->read_index(this->model.number_objectives, "objective");
               const int sense = this->tokenizer.read_int();
               this->tokenizer.end_item();
               const ExpressionArena::Mark mark = this->builder.expression_arena().mark();
               const NodeIndex root = this->read_expression();
               this->builder.set_objective(objective, sense != 0, root, mark);
            }

            void read_defined_variable() {
               const int index = this->tokenizer.read_int();
               const int number_linear_terms = this->tokenizer.read_int();
               this->tokenizer.read_int(); // kind (where the variable is used): not needed
               this->tokenizer.end_item();
               if (index < this->number_variables || index >= this->number_expression_variables) {
                  this->tokenizer.fail("invalid defined variable index");
               }
               this->linear_scratch.clear();
               for (int k = 0; k < number_linear_terms; ++k) {
                  const int variable = this->read_index(this->number_variables, "variable");
                  const double coefficient = this->tokenizer.read_double();
                  this->tokenizer.end_item();
                  this->linear_scratch.push_back({variable, coefficient});
               }
               const NodeIndex root = this->read_expression();
               this->builder.define_variable(index - this->number_variables, this->linear_scratch.data(),
                  this->linear_scratch.size(), root);
            }

            void read_jacobian_row() {
               const int constraint = this->read_index(this->model.number_constraints, "constraint");
               const int count = this->tokenizer.read_int();
               this->tokenizer.end_item();
               for (int k = 0; k < count; ++k) {
                  const int variable = this->read_index(this->number_variables, "variable");
                  if (!this->column_counts_seen) this->tokenizer.read_int(); // old format: explicit goff
                  const double coefficient = this->tokenizer.read_double();
                  this->tokenizer.end_item();
                  this->builder.add_jacobian_entry(constraint, variable, coefficient);
               }
            }

            void read_objective_gradient() {
               const int objective = this->read_index(this->model.number_objectives, "objective");
               const int count = this->tokenizer.read_int();
               this->tokenizer.end_item();
               for (int k = 0; k < count; ++k) {
                  const int variable = this->read_index(this->number_variables, "variable");
                  const double coefficient = this->tokenizer.read_double();
                  this->tokenizer.end_item();
                  this->builder.add_objective_gradient_entry(objective, variable, coefficient);
               }
            }

            /// cumulative Jacobian column counts: the Jacobian pattern is rebuilt row-wise, so they are only checked
            void read_column_counts() {
               const int count = this->tokenizer.read_int();
               this->tokenizer.end_item();
               if (count != this->number_variables - 1) this->tokenizer.fail("invalid column count segment");
               for (int k = 0; k < count; ++k) {
                  this->tokenizer.read_long();
                  this->tokenizer.end_item();
               }
               this->column_counts_seen = true;
            }

            void read_bounds(std::vector<double>& lower, std::vector<double>& upper, bool is_constraint) {
               this->tokenizer.end_item();
               const std::size_t count = lower.size();
               for (std::size_t i = 0; i < count; ++i) {
                  const int kind = this->tokenizer.read_letter() - '0';
                  switch (kind) {
                     case 0: lower[i] = this->tokenizer.read_double(); upper[i] = this->tokenizer.read_double(); break;
                     case 1: lower[i] = -infinity; upper[i] = this->tokenizer.read_double(); break;
                     case 2: lower[i] = this->tokenizer.read_double(); upper[i] = infinity; break;
                     case 3: lower[i] = -infinity; upper[i] = infinity; break;
                     case 4: lower[i] = upper[i] = this->tokenizer.read_double(); break;
                     case 5: {
                        if (!is_constraint) this->tokenizer.fail("complementarity in variable bounds");
                        const int flags = this->tokenizer.read_int();
                        const int variable = this->tokenizer.read_int();
                        if (variable < 1 || variable > this->number_variables) this->tokenizer.fail("invalid complementary variable");
                        this->model.complementary_variables[i] = variable - 1;
                        lower[i] = (flags & 2) ? -infinity : 0.;
                        upper[i] = (flags & 1) ? infinity : 0.;
                        break;
                     }
                     default: this->tokenizer.fail("invalid bound type");
                  }
                  this->tokenizer.end_item();
               }
            }
            void read_constraint_bounds() {
               this->read_bounds(this->model.constraint_lower_bounds, this->model.constraint_upper_bounds, true);
            }
            void read_variable_bounds() {
               this->read_bounds(this->model.variable_lower_bounds, this->model.variable_upper_bounds, false);
            }

            void read_initial_point(std::vector<double>& point, int dimension) {
               const int count = this->tokenizer.read_int();
               this->tokenizer.end_item();
               for (int k = 0; k < count; ++k) {
                  const int index = this->read_index(dimension, "initial point");
                  point[static_cast<std::size_t>(index)] = this->tokenizer.read_double();
                  this->tokenizer.end_item();
               }
            }

            void read_suffix() {
               const int kind = this->tokenizer.read_int();
               const int count = this->tokenizer.read_int();
               Suffix suffix;
               suffix.name = this->tokenizer.read_name();
               this->tokenizer.end_item();
               if (kind < 0 || kind > 7 || count <= 0) this->tokenizer.fail("invalid suffix header");
               suffix.target = static_cast<SuffixTarget>(kind & 3);
               suffix.is_real = (kind & 4) != 0;
               const int dimensions[4] = {this->number_variables, this->model.number_constraints + this->header.number_logical_constraints,
                  this->model.number_objectives, 1};
               const int dimension = dimensions[kind & 3];
               suffix.values.assign(static_cast<std::size_t>(dimension), 0.);
               for (int k = 0; k < count; ++k) {
                  const int index = this->read_index(dimension, "suffix");
                  suffix.values[static_cast<std::size_t>(index)] =
                     suffix.is_real ? this->tokenizer.read_double() : static_cast<double>(this->tokenizer.read_long());
                  this->tokenizer.end_item();
               }
               this->model.suffixes.push_back(std::move(suffix));
            }
         };

         template <typename Tokenizer>
         std::shared_ptr<const ModelData> read_body(const char* begin, const char* body, const char* end,
               const NlHeader& header, const ReaderOptions& options) {
            ModelBuilder builder(header, options);
            SegmentReader<Tokenizer> reader(Tokenizer(begin, body, end), header, builder);
            reader.read_all_segments();
            return builder.finalize();
         }

         std::shared_ptr<const ModelData> read_model(const char* contents, std::size_t size, const ReaderOptions& options) {
            const char* cursor = contents;
            const char* end = contents + size;
            const NlHeader header = parse_header(cursor, end);
            if (header.number_logical_constraints > 0) parse_error(0, "logical constraints are not supported");
            switch (header.format) {
               case NlFormat::Text:
                  return read_body<TextTokenizer>(contents, cursor, end, header, options);
               case NlFormat::Binary:
                  return read_body<BinaryTokenizer<std::int32_t, std::int32_t>>(contents, cursor, end, header, options);
               case NlFormat::BinaryShortOpcodes:
                  return read_body<BinaryTokenizer<std::int32_t, std::int16_t>>(contents, cursor, end, header, options);
               case NlFormat::BinaryLongIntegers:
                  return read_body<BinaryTokenizer<std::int64_t, std::int16_t>>(contents, cursor, end, header, options);
            }
            parse_error(0, "unknown format");
         }
      } // namespace
   } // namespace detail

   NlModel NlModel::read_from_memory(const char* contents, std::size_t size, const ReaderOptions& options) {
      return NlModel(detail::read_model(contents, size, options));
   }

   NlModel NlModel::read_file(const std::string& file_name, const ReaderOptions& options) {
#ifdef CPPASL_HAS_MMAP
      const int descriptor = ::open(file_name.c_str(), O_RDONLY);
      if (descriptor < 0) throw std::runtime_error("cppasl: cannot open " + file_name);
      struct stat status {};
      if (::fstat(descriptor, &status) != 0 || status.st_size <= 0) {
         ::close(descriptor);
         throw std::runtime_error("cppasl: cannot read " + file_name);
      }
      const auto size = static_cast<std::size_t>(status.st_size);
      void* mapping = ::mmap(nullptr, size, PROT_READ, MAP_PRIVATE, descriptor, 0);
      ::close(descriptor);
      if (mapping == MAP_FAILED) throw std::runtime_error("cppasl: cannot map " + file_name);
      ::madvise(mapping, size, MADV_SEQUENTIAL);
      struct Unmap {
         void* address;
         std::size_t size;
         ~Unmap() { ::munmap(this->address, this->size); }
      } unmap{mapping, size};
      return NlModel::read_from_memory(static_cast<const char*>(mapping), size, options);
#else
      std::FILE* file = std::fopen(file_name.c_str(), "rb");
      if (file == nullptr) throw std::runtime_error("cppasl: cannot open " + file_name);
      std::vector<char> contents;
      char chunk[1 << 16];
      std::size_t read = 0;
      while ((read = std::fread(chunk, 1, sizeof(chunk), file)) > 0) contents.insert(contents.end(), chunk, chunk + read);
      std::fclose(file);
      return NlModel::read_from_memory(contents.data(), contents.size(), options);
#endif
   }

} // namespace cppasl
