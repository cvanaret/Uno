# Script to parse Uno headers and generate Fortran interfaces.
using Clang
using Clang.Generators
using Clang.LibClang

# Set of callback typedef names (used to map to type(c_funptr) in function arguments)
const CALLBACKS = String[
  "uno_objective_callback",
  "uno_objective_gradient_callback",
  "uno_constraints_callback",
  "uno_constraints_jacobian_callback",
  "uno_lagrangian_hessian_callback",
  "uno_constraints_jacobian_operator_callback",
  "uno_constraints_jacobian_transposed_operator_callback",
  "uno_lagrangian_hessian_operator_callback",
  "uno_notify_acceptable_iterate_callback",
  "uno_termination_callback",
  "uno_logger_stream_callback",
]

# ----------------------------------------------------------------
# Hardcoded prologue for uno_c.f90:
#   - uno_int definition
#   - all constants
#   - all abstract interfaces for callbacks
#
# These require semantic knowledge that Clang cannot fully provide.
# ----------------------------------------------------------------

# ----------------------------------------------------------------
# Pre-computed type info structs (must be defined before any use)
# ----------------------------------------------------------------
struct FortranArgInfo
  name      :: String
  ftype     :: String   # e.g. "real(c_double)"
  suffix    :: String   # e.g. ", value" or "(*)" or ""
  is_string :: Bool     # true if original C type is char*
end

struct FortranFuncInfo
  name      :: String
  ret_ftype :: String   # e.g. "logical(c_bool)" or "void"
  is_void   :: Bool
  args      :: Vector{FortranArgInfo}
end

const C_TO_FORTRAN_KIND = Dict(
  "int8_t"   => "c_int8_t",
  "int16_t"  => "c_int16_t",
  "int32_t"  => "c_int32_t",
  "int64_t"  => "c_int64_t",
  "uint8_t"  => "c_int8_t",
  "uint16_t" => "c_int16_t",
  "uint32_t" => "c_int32_t",
  "uint64_t" => "c_int64_t",
  "int"      => "c_int",
  "long"     => "c_long",
)

function extract_uno_int_kind(include_dir)
  uno_int_h = joinpath(include_dir, "uno_int.h")
  text = read(uno_int_h, String)
  m = match(r"typedef\s+(\w+)\s+uno_int\s*;", text)
  m === nothing && error("Cannot find 'typedef ... uno_int' in $uno_int_h")
  c_type = m[1]
  haskey(C_TO_FORTRAN_KIND, c_type) || error("Unknown C type '$c_type' for uno_int")
  return C_TO_FORTRAN_KIND[c_type]
end

function is_string_arg(t)
  k = Clang.kind(t)
  if k == Clang.CXType_Pointer
    pt = Clang.getPointeeType(t)
    pk = Clang.kind(pt)
    return pk in (Clang.CXType_Char_S, Clang.CXType_Char_U)
  end
  return false
end

# Write a Fortran statement, wrapping long lines at max_col characters using '&'
const FORTRAN_MAX_COL = 90

function println_wrapped(io, line::String, continuation_indent::String = "   ")
  if length(line) <= FORTRAN_MAX_COL
    println(io, line)
    return
  end
  pos = FORTRAN_MAX_COL - 2
  while pos > 1 && line[pos] != ','
    pos -= 1
  end
  if pos <= 1
    println(io, line)
    return
  end
  println(io, line[1:pos] * " &")
  println_wrapped(io, continuation_indent * lstrip(line[pos+1:end]), continuation_indent)
end

# Emit a function/subroutine statement, handling line breaks correctly.
# Breaks at result(...) if the full line is too long but args alone fit,
# then falls back to comma-breaking within the arg list.
# cont_indent : indent for bind(C,...) and result(...) continuation lines
# with_amp    : whether to append ' &' (set false for wrapper subprograms)
function emit_subprogram_header(io, base_indent, name, joined, rname, cont_indent; with_amp=true)
  kind        = rname !== nothing ? "function" : "subroutine"
  args_indent = " " ^ length("$base_indent$kind $name(")
  amp         = with_amp ? " &" : ""

  if rname !== nothing && with_amp
    # bind(C,...) follows: always put result(...) on its own line for consistency
    line1 = "$base_indent$kind $name($joined) &"
    println_wrapped(io, line1, args_indent)
    println(io, "$(cont_indent)result($rname)$amp")
  else
    # Wrapper (no bind): break only when needed
    result_clause = rname !== nothing ? " result($rname)" : ""
    full = "$base_indent$kind $name($joined)$result_clause$amp"
    if length(full) <= FORTRAN_MAX_COL
      println(io, full)
      return
    end
    line1 = "$base_indent$kind $name($joined) &"
    println_wrapped(io, line1, args_indent)
    rname !== nothing && println(io, "$(cont_indent)result($rname)$amp")
  end
end

function name_to_result(fname, ret_ftype="")
  r = replace(fname, r"^uno_" => "")
  # uno_get_xxx → strip get_ prefix (takes priority over success naming)
  r2 = replace(r, r"^get_" => "")
  r2 != r && return r2
  ret_ftype == "logical(c_bool)" && return "success"
  # functions returning a handle: keep only the noun
  if ret_ftype == "type(c_ptr)"
    for noun in ("solver", "model")
      endswith(r, noun) && return noun
    end
  end
  return r
end

# Extract the import name(s) from a Fortran type string
# e.g. "real(c_double)" → "c_double", "type(c_ptr)" → "c_ptr"
function ftype_to_import(ftype::String)
  m = match(r"\((\w+)\)", ftype)
  m === nothing && return nothing
  return m[1]
end

# Build a minimal import list in order of first appearance in arguments.
# Works on pre-computed FortranArgInfo strings (no Clang calls needed).
function build_import(args::Vector{FortranArgInfo}, ret_ftype, is_void)
  order = String[]
  seen  = Set{String}()
  add = function(k)
    k !== nothing && k ∉ seen && (push!(order, k); push!(seen, k))
  end
  for a in args
    ftype = a.is_string ? "character(c_char)" : a.ftype
    add(ftype_to_import(ftype))
  end
  is_void || add(ftype_to_import(ret_ftype))
  return join(order, ", ")
end

# ----------------------------------------------------------------
# Map CLType → Fortran type string (for return values and wrappers)
# ----------------------------------------------------------------
function cltype_to_fortran(t)
  k = Clang.kind(t)
  spell = Clang.spelling(t)

  if k == Clang.CXType_Elaborated
    if spell in CALLBACKS
      return ("type(c_funptr)", true)
    end
    ct = Clang.getCanonicalType(t)
    ck = Clang.kind(ct)
    if ck in (Clang.CXType_Int, Clang.CXType_UInt,
               Clang.CXType_Long, Clang.CXType_ULong,
               Clang.CXType_LongLong, Clang.CXType_ULongLong)
      return ("integer(uno_int)", true)
    end
    if ck == Clang.CXType_Double;  return ("real(c_double)", true);   end
    if ck == Clang.CXType_Bool;    return ("logical(c_bool)", true);   end
    if ck == Clang.CXType_Pointer
      pt = Clang.getPointeeType(ct)
      if Clang.kind(pt) == Clang.CXType_FunctionProto
        return ("type(c_funptr)", true)
      end
      return ("type(c_ptr)", true)
    end
    return ("type(c_ptr)", true)
  end

  if k in (Clang.CXType_Int, Clang.CXType_UInt); return ("integer(c_int)", true); end
  if k == Clang.CXType_Double;  return ("real(c_double)", true);  end
  if k == Clang.CXType_Bool;    return ("logical(c_bool)", true); end
  if k in (Clang.CXType_Char_S, Clang.CXType_Char_U)
    return ("character(c_char)", true)
  end
  if k == Clang.CXType_Void; return ("void", true); end

  if k == Clang.CXType_Pointer
    pt = Clang.getPointeeType(t)
    pk = Clang.kind(pt)
    if pk == Clang.CXType_Void;   return ("type(c_ptr)", true);    end
    if pk == Clang.CXType_Double; return ("real(c_double)", false); end
    if pk in (Clang.CXType_Char_S, Clang.CXType_Char_U)
      return ("character(c_char)", false)
    end
    if pk == Clang.CXType_FunctionProto; return ("type(c_funptr)", true); end
    if pk == Clang.CXType_Elaborated
      ptspell = Clang.spelling(pt)
      if ptspell in CALLBACKS; return ("type(c_funptr)", true); end
      ct = Clang.getCanonicalType(pt)
      if Clang.kind(ct) in (Clang.CXType_Int, Clang.CXType_UInt)
        return ("integer(uno_int)", false)
      end
      if Clang.kind(ct) == Clang.CXType_FunctionProto
        return ("type(c_funptr)", true)
      end
    end
    return ("type(c_ptr)", true)
  end

  return ("type(c_ptr)", true)
end

# ----------------------------------------------------------------
# Emit arg declarations from pre-computed FortranArgInfo (no Clang calls)
# ----------------------------------------------------------------
function emit_bind_c_arg(io, indent, a::FortranArgInfo)
  s = a.suffix
  if s == "(array)"
    println(io, "$(indent)$(a.ftype) :: $(a.name)(*)")
  elseif s == "(scalar)"
    println(io, "$(indent)$(a.ftype) :: $(a.name)")
  elseif s == "(*)"
    println(io, "$(indent)$(a.ftype) :: $(a.name)(*)")
  else
    println(io, "$(indent)$(a.ftype)$s :: $(a.name)")
  end
end

function emit_wrapper_arg(io, a::FortranArgInfo)
  s = a.suffix
  if s in ("(*)", "(array)", "(scalar)")
    println(io, "   $(a.ftype) :: $(a.name)$(s == "(scalar)" ? "" : "(*)")")
  else
    println(io, "   $(a.ftype)$s :: $(a.name)")
  end
end

# ----------------------------------------------------------------
# Collect function declarations in ONE pass over the cursor tree
# ----------------------------------------------------------------
# Compute FortranArgInfo from a cursor arg — must be called while Clang ctx is alive
function make_arg_info(arg_cursor, arg_type)
  aname = Clang.spelling(arg_cursor)
  t = arg_type
  k = Clang.kind(t)
  spell = Clang.spelling(t)

  is_str = is_string_arg(t)

  # Determine ftype and suffix
  ftype, suffix = if k in (Clang.CXType_Char_S, Clang.CXType_Char_U)
    "character(c_char)", ", value"
  elseif k == Clang.CXType_Bool
    "logical(c_bool)", ", value"
  elseif k == Clang.CXType_Double
    "real(c_double)", ", value"
  elseif k in (Clang.CXType_Int, Clang.CXType_UInt)
    "integer(c_int)", ", value"
  elseif k == Clang.CXType_Elaborated
    if spell in CALLBACKS
      "type(c_funptr)", ", value"
    else
      ct = Clang.getCanonicalType(t)
      ck = Clang.kind(ct)
      if ck in (Clang.CXType_Int, Clang.CXType_UInt)
        "integer(uno_int)", ", value"
      elseif ck == Clang.CXType_Bool
        "logical(c_bool)", ", value"
      elseif ck == Clang.CXType_Double
        "real(c_double)", ", value"
      elseif ck == Clang.CXType_Pointer
        pt = Clang.getPointeeType(ct)
        if Clang.kind(pt) == Clang.CXType_FunctionProto
          "type(c_funptr)", ", value"
        else
          "type(c_ptr)", ", value"
        end
      else
        "type(c_ptr)", ", value"
      end
    end
  elseif k == Clang.CXType_Pointer
    pt = Clang.getPointeeType(t)
    pk = Clang.kind(pt)
    if pk == Clang.CXType_Void
      "type(c_ptr)", ", value"
    elseif pk == Clang.CXType_Double
      "real(c_double)", ""   # array: suffix applied as (*) below
    elseif pk in (Clang.CXType_Char_S, Clang.CXType_Char_U)
      "character(c_char)", ""
    elseif pk == Clang.CXType_FunctionProto
      "type(c_funptr)", ", value"
    elseif pk == Clang.CXType_Elaborated
      ptspell = Clang.spelling(pt)
      if ptspell in CALLBACKS
        "type(c_funptr)", ", value"
      else
        ct = Clang.getCanonicalType(pt)
        if Clang.kind(ct) in (Clang.CXType_Int, Clang.CXType_UInt)
          # const uno_int* → array input; mutable → scalar output
          if occursin("const", spell)
            "integer(uno_int)", ""
          else
            "integer(uno_int)", ""   # scalar: no (*) added for scalars by convention below
          end
        elseif Clang.kind(ct) == Clang.CXType_FunctionProto
          "type(c_funptr)", ", value"
        else
          "type(c_ptr)", ", value"
        end
      end
    else
      "type(c_ptr)", ", value"
    end
  else
    "type(c_ptr)", ", value"
  end

  # For pointer-to-int: const → array, non-const → scalar
  final_suffix = if ftype == "integer(uno_int)" && suffix == "" && k == Clang.CXType_Pointer
    occursin("const", spell) ? "(array)" : "(scalar)"
  elseif suffix == ""
    "(*)"
  else
    suffix
  end

  FortranArgInfo(aname, ftype, final_suffix, is_str)
end

function collect_funcs(root)
  funcs = FortranFuncInfo[]
  for child in Clang.children(root)
    k = Clang.kind(child)
    name = Clang.spelling(child)
    if k == Clang.CXCursor_FunctionDecl && startswith(name, "uno_")
      ret_cltype = Clang.getCursorResultType(child)
      is_void = Clang.kind(ret_cltype) == Clang.CXType_Void
      ret_ftype = is_void ? "void" : first(cltype_to_fortran(ret_cltype))
      clang_args = Clang.get_function_args(child)
      clang_types = [Clang.getCursorType(a) for a in clang_args]
      args = [make_arg_info(clang_args[i], clang_types[i]) for i in eachindex(clang_args)]
      push!(funcs, FortranFuncInfo(name, ret_ftype, is_void, args))
    end
  end
  return funcs
end

# ----------------------------------------------------------------
# Generate uno_c.f90
# ----------------------------------------------------------------
function gen_uno_c(io, include_dir, funcs)
  println(io, "! Copyright (c) 2026 Alexis Montoison and Charlie Vanaret")
  println(io, "! Licensed under the MIT license. See LICENSE file in the project directory for details.")
  println(io, "")
  println(io, "!==============================================================")
  println(io, "! Fortran interfaces -- uno_c.f90")
  println(io, "!==============================================================")
  println(io, "")

  # Prologue: read prologue_fortran.f90 and inject uno_int at the marker
  uno_int_kind = extract_uno_int_kind(include_dir)
  prologue = read(joinpath(@__DIR__, "prologue_fortran.f90"), String)
  marker = "!--- UNO_INT_KIND ---\n"
  idx = findfirst(marker, prologue)
  idx === nothing && error("Marker '$marker' not found in prologue_fortran.f90")
  print(io, prologue[1:first(idx)-1])
  println(io, "integer, parameter :: uno_int = $uno_int_kind")
  println(io, "")
  print(io, prologue[last(idx)+1:end])
  println(io, "")

  # C bind(C) interfaces — skip functions that have string args
  non_string = [f for f in funcs if !any(a.is_string for a in f.args)]
  for (i, f) in enumerate(non_string)
    gen_one_c_interface(io, f, i < length(non_string))
  end
end

function gen_one_c_interface(io, f::FortranFuncInfo, trailing_blank=true)
  name   = f.name
  joined = join((a.name for a in f.args), ", ")
  rname  = name_to_result(name, f.ret_ftype)
  imports = build_import(f.args, f.ret_ftype, f.is_void)

  println(io, "!---------------------------------------------")
  println(io, "! $name")
  println(io, "!---------------------------------------------")
  println(io, "interface")
  emit_subprogram_header(io, "   ", name, joined, f.is_void ? nothing : rname, "      ")
  println(io, "      bind(C, name=\"$name\")")
  !isempty(imports) && println(io, "      import :: $imports")

  for a in f.args
    emit_bind_c_arg(io, "      ", a)
  end

  if !f.is_void
    println(io, "      $(f.ret_ftype) :: $rname")
    println(io, "   end function $name")
  else
    println(io, "   end subroutine $name")
  end
  println(io, "end interface")
  trailing_blank && println(io, "")
end

# ----------------------------------------------------------------
# Generate uno_fortran.f90: Fortran-friendly wrappers for string arguments
# ----------------------------------------------------------------
function gen_uno_fortran(io, funcs)
  println(io, "! Copyright (c) 2026 Alexis Montoison and Charlie Vanaret")
  println(io, "! Licensed under the MIT license. See LICENSE file in the project directory for details.")
  println(io, "")
  println(io, "!==============================================================")
  println(io, "! Fortran interfaces -- uno_fortran.f90")
  println(io, "!==============================================================")

  for f in funcs
    string_arg_names = [a.name for a in f.args if a.is_string]
    if isempty(string_arg_names)
      continue
    end
    if f.name == "uno_get_solver_string_option"
      continue
    end
    gen_one_string_wrapper(io, f, string_arg_names)
  end

  # Special wrapper for uno_get_solver_string_option (returns const char*)
  gen_get_string_option_fortran_wrapper(io)
end

function gen_one_string_wrapper(io, f::FortranFuncInfo, string_arg_names)
  name        = f.name
  joined      = join((a.name for a in f.args), ", ")
  result_name = name_to_result(name, f.ret_ftype)
  c_name      = name * "_c"
  imports     = build_import(f.args, f.ret_ftype, f.is_void)

  println(io, "")
  println(io, "!---------------------------------------------")
  println(io, "! $name")
  println(io, "!---------------------------------------------")

  emit_subprogram_header(io, "", name, joined, f.is_void ? nothing : result_name, "   "; with_amp=false)

  for a in f.args
    if a.is_string
      println(io, "   character(len=*) :: $(a.name)")
    else
      emit_wrapper_arg(io, a)
    end
  end

  !f.is_void && println(io, "   $(f.ret_ftype) :: $result_name")

  for sarg in string_arg_names
    println(io, "   character(c_char), allocatable :: $(sarg)_c(:)")
  end
  println(io, "   integer :: i, n")
  println(io, "")

  # Inner C interface
  println(io, "   interface")
  if f.is_void
    emit_subprogram_header(io, "      ", c_name, joined, nothing, "         ")
    println(io, "         bind(C, name=\"$name\")")
    !isempty(imports) && println(io, "         import :: $imports")
    for a in f.args
      if a.is_string
        println(io, "         character(c_char) :: $(a.name)(*)")
      else
        emit_bind_c_arg(io, "         ", a)
      end
    end
    println(io, "      end subroutine $c_name")
  else
    emit_subprogram_header(io, "      ", c_name, joined, result_name, "         ")
    println(io, "         bind(C, name=\"$name\")")
    !isempty(imports) && println(io, "         import :: $imports")
    for a in f.args
      if a.is_string
        println(io, "         character(c_char) :: $(a.name)(*)")
      else
        emit_bind_c_arg(io, "         ", a)
      end
    end
    println(io, "         $(f.ret_ftype) :: $result_name")
    println(io, "      end function $c_name")
  end
  println(io, "   end interface")
  println(io, "")

  # Convert each string arg to C char array
  for sarg in string_arg_names
    println(io, "   n = len_trim($sarg)")
    println(io, "   allocate($(sarg)_c(n+1))")
    println(io, "   do i = 1, n")
    println(io, "      $(sarg)_c(i) = $sarg(i:i)")
    println(io, "   end do")
    println(io, "   $(sarg)_c(n+1) = c_null_char")
  end

  c_args   = [a.is_string ? "$(a.name)_c" : a.name for a in f.args]
  c_joined = join(c_args, ", ")

  if f.is_void
    println_wrapped(io, "   call $c_name($c_joined)", " " ^ (length("   call $c_name(")))
  else
    println_wrapped(io, "   $result_name = $c_name($c_joined)", " " ^ (length("   $result_name = $c_name(")))
  end

  for sarg in string_arg_names
    println(io, "   deallocate($(sarg)_c)")
  end

  f.is_void ? println(io, "end subroutine $name") : println(io, "end function $name")
end

function gen_get_string_option_fortran_wrapper(io)
  println(io, "")
  println(io, "!---------------------------------------------")
  println(io, "! uno_get_solver_string_option")
  println(io, "!---------------------------------------------")
  println(io, "function uno_get_solver_string_option(solver, option_name) &")
  println(io, "   result(solver_string_option)")
  println(io, "   type(c_ptr), value :: solver")
  println(io, "   character(len=*) :: option_name")
  println(io, "   character(:), allocatable :: solver_string_option")
  println(io, "   character(c_char), allocatable :: option_name_c(:)")
  println(io, "   type(c_ptr) :: ptr_solver_string_option_c")
  println(io, "   integer :: i, n")
  println(io, "   character(c_char), pointer :: solver_string_option_c(:)")
  println(io, "")
  println(io, "   interface")
  println(io, "      function uno_get_solver_string_option_c(solver, option_name) &")
  println(io, "         result(solver_string_option) &")
  println(io, "         bind(C, name=\"uno_get_solver_string_option\")")
  println(io, "         import :: c_ptr, c_char")
  println(io, "         type(c_ptr), value :: solver")
  println(io, "         character(c_char) :: option_name(*)")
  println(io, "         type(c_ptr) :: solver_string_option")
  println(io, "      end function uno_get_solver_string_option_c")
  println(io, "   end interface")
  println(io, "")
  println(io, "   n = len_trim(option_name)")
  println(io, "   allocate(option_name_c(n+1))")
  println(io, "   do i = 1, n")
  println(io, "      option_name_c(i) = option_name(i:i)")
  println(io, "   end do")
  println(io, "   option_name_c(n+1) = c_null_char")
  println(io, "   ptr_solver_string_option_c = uno_get_solver_string_option_c(solver, option_name_c)")
  println(io, "   call c_f_pointer(ptr_solver_string_option_c, solver_string_option_c, [0])")
  println(io, "   n = 0")
  println(io, "   do while (solver_string_option_c(n+1) /= c_null_char)")
  println(io, "      n = n + 1")
  println(io, "   end do")
  println(io, "   allocate(character(len=n)::solver_string_option)")
  println(io, "   do i = 1, n")
  println(io, "      solver_string_option(i:i) = solver_string_option_c(i)")
  println(io, "   end do")
  println(io, "   deallocate(option_name_c)")
  println(io, "end function uno_get_solver_string_option")
end

function main_fortran()
  include_dir = joinpath(@__DIR__, "..", "..", "C")
  header_path = joinpath(include_dir, "Uno_C_API.h")

  args = get_default_args()
  push!(args, "-I$include_dir")

  ctx = create_context([header_path], args)
  tu = ctx.trans_units[1]
  root = Clang.getTranslationUnitCursor(tu)

  # Single pass — keep ctx alive via GC.@preserve
  funcs = GC.@preserve ctx tu collect_funcs(root)

  fortran_dir = joinpath(@__DIR__, "..", "..", "Fortran")

  out_c = joinpath(fortran_dir, "uno_c.f90")
  open(out_c, "w") do io
    gen_uno_c(io, include_dir, funcs)
  end
  println("Generated: $out_c")

  out_fortran = joinpath(fortran_dir, "uno_fortran.f90")
  open(out_fortran, "w") do io
    gen_uno_fortran(io, funcs)
  end
  println("Generated: $out_fortran")
end

# If we want to use the file as a script with `julia wrapper_fortran.jl`
if abspath(PROGRAM_FILE) == @__FILE__
  main_fortran()
end
