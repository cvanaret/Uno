# Script to parse Uno headers and generate Rust FFI bindings (ffi.rs).
using Clang
using Clang.Generators
using Clang.LibClang

# Callback typedef names — excluded from extern "C" (defined in prologue_rust.rs)
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

# Callbacks that may be passed as NULL (wrapped in Option<> in Rust)
const OPTIONAL_CALLBACKS = Set{String}([
  "uno_notify_acceptable_iterate_callback",
  "uno_termination_callback",
])

const C_TO_RUST_INT = Dict(
  "int8_t"   => "i8",
  "int16_t"  => "i16",
  "int32_t"  => "i32",
  "int64_t"  => "i64",
  "uint8_t"  => "u8",
  "uint16_t" => "u16",
  "uint32_t" => "u32",
  "uint64_t" => "u64",
  "int"      => "c_int",
  "long"     => "c_long",
)

# ----------------------------------------------------------------
# Read uno_int.h and return the Rust primitive type for uno_int
# ----------------------------------------------------------------
function extract_uno_int_rust(include_dir)
  uno_int_h = joinpath(include_dir, "uno_int.h")
  text = read(uno_int_h, String)
  m = match(r"typedef\s+(\w+)\s+uno_int\s*;", text)
  m === nothing && error("Cannot find 'typedef ... uno_int' in $uno_int_h")
  c_type = m[1]
  haskey(C_TO_RUST_INT, c_type) || error("Unknown C type '$c_type' for uno_int")
  return C_TO_RUST_INT[c_type]
end

# ----------------------------------------------------------------
# Map a Clang type → Rust type string for extern "C" declarations
# ----------------------------------------------------------------
function cltype_to_rust(t, spell::String)
  k = Clang.kind(t)

  # Elaborated types: uno_int or callback typedefs
  if k == Clang.CXType_Elaborated
    if spell in CALLBACKS
      return spell in OPTIONAL_CALLBACKS ? "Option<$spell>" : spell
    end
    ct = Clang.getCanonicalType(t)
    ck = Clang.kind(ct)
    if ck in (Clang.CXType_Int, Clang.CXType_UInt,
               Clang.CXType_Long, Clang.CXType_ULong,
               Clang.CXType_LongLong, Clang.CXType_ULongLong)
      return "uno_int"
    end
    if ck == Clang.CXType_Double;  return "f64";  end
    if ck == Clang.CXType_Bool;    return "bool"; end
    if ck == Clang.CXType_Pointer
      pt = Clang.getPointeeType(ct)
      Clang.kind(pt) == Clang.CXType_FunctionProto && return spell
      return "*mut c_void"
    end
    return "*mut c_void"
  end

  # Primitive by-value types
  if k in (Clang.CXType_Int, Clang.CXType_UInt); return "c_int"; end
  if k == Clang.CXType_Double;                    return "f64";   end
  if k == Clang.CXType_Bool;                      return "bool";  end
  if k in (Clang.CXType_Char_S, Clang.CXType_Char_U); return "c_char"; end

  # Pointer types
  if k == Clang.CXType_Pointer
    pt       = Clang.getPointeeType(t)
    pk       = Clang.kind(pt)
    is_const = occursin("const", spell)

    if pk == Clang.CXType_Void
      return is_const ? "*const c_void" : "*mut c_void"
    end
    if pk == Clang.CXType_Double
      return is_const ? "*const f64" : "*mut f64"
    end
    if pk in (Clang.CXType_Char_S, Clang.CXType_Char_U)
      return is_const ? "*const c_char" : "*mut c_char"
    end
    if pk == Clang.CXType_Bool
      return is_const ? "*const bool" : "*mut bool"
    end
    if pk in (Clang.CXType_Int, Clang.CXType_UInt)
      return is_const ? "*const c_int" : "*mut c_int"
    end
    if pk == Clang.CXType_Elaborated
      ptspell = Clang.spelling(pt)
      if ptspell in CALLBACKS
        return ptspell in OPTIONAL_CALLBACKS ? "Option<$ptspell>" : ptspell
      end
      ct  = Clang.getCanonicalType(pt)
      ck  = Clang.kind(ct)
      if ck in (Clang.CXType_Int, Clang.CXType_UInt,
                 Clang.CXType_Long, Clang.CXType_ULong,
                 Clang.CXType_LongLong, Clang.CXType_ULongLong)
        return is_const ? "*const uno_int" : "*mut uno_int"
      end
      if ck == Clang.CXType_FunctionProto
        return ptspell in OPTIONAL_CALLBACKS ? "Option<$ptspell>" : ptspell
      end
    end
    if pk == Clang.CXType_FunctionProto
      return "*mut c_void"  # shouldn't happen for Uno
    end
    return is_const ? "*const c_void" : "*mut c_void"
  end

  return "*mut c_void"
end

# ----------------------------------------------------------------
# Map a Clang return type → Rust return type string (nothing = void)
# ----------------------------------------------------------------
function rettype_to_rust(t)
  k = Clang.kind(t)
  if k == Clang.CXType_Void;   return nothing; end
  if k == Clang.CXType_Bool;   return "bool";  end
  if k == Clang.CXType_Double; return "f64";   end
  if k in (Clang.CXType_Int, Clang.CXType_UInt); return "c_int"; end

  if k == Clang.CXType_Elaborated
    ct = Clang.getCanonicalType(t)
    ck = Clang.kind(ct)
    if ck in (Clang.CXType_Int, Clang.CXType_UInt,
               Clang.CXType_Long, Clang.CXType_ULong,
               Clang.CXType_LongLong, Clang.CXType_ULongLong)
      return "uno_int"
    end
    if ck == Clang.CXType_Bool;   return "bool"; end
    if ck == Clang.CXType_Double; return "f64";  end
    if ck == Clang.CXType_Pointer; return "*mut c_void"; end
  end

  if k == Clang.CXType_Pointer
    pt       = Clang.getPointeeType(t)
    pk       = Clang.kind(pt)
    spell    = Clang.spelling(t)
    is_const = occursin("const", spell)
    if pk == Clang.CXType_Void
      return is_const ? "*const c_void" : "*mut c_void"
    end
    if pk in (Clang.CXType_Char_S, Clang.CXType_Char_U)
      return is_const ? "*const c_char" : "*mut c_char"
    end
    if pk == Clang.CXType_Double
      return is_const ? "*const f64" : "*mut f64"
    end
  end

  return "*mut c_void"
end

# ----------------------------------------------------------------
# Pre-computed structs
# ----------------------------------------------------------------
struct ArgInfo
  name  :: String
  rtype :: String
end

struct FuncInfo
  name :: String
  ret  :: Union{String, Nothing}  # nothing = void
  args :: Vector{ArgInfo}
end

# ----------------------------------------------------------------
# Collect all uno_* function declarations in one pass
# ----------------------------------------------------------------
function collect_funcs(root)
  funcs = FuncInfo[]
  for child in Clang.children(root)
    k    = Clang.kind(child)
    name = Clang.spelling(child)
    if k == Clang.CXCursor_FunctionDecl && startswith(name, "uno_")
      ret_t  = Clang.getCursorResultType(child)
      ret    = rettype_to_rust(ret_t)
      clargs = Clang.get_function_args(child)
      args   = map(clargs) do a
        aname = Clang.spelling(a)
        at    = Clang.getCursorType(a)
        spell = Clang.spelling(at)
        ArgInfo(aname, cltype_to_rust(at, spell))
      end
      push!(funcs, FuncInfo(name, ret, collect(args)))
    end
  end
  return funcs
end

# ----------------------------------------------------------------
# Emit one extern "C" function declaration
# ----------------------------------------------------------------
const RUST_MAX_COL = 100

function gen_one_fn(io, f::FuncInfo, trailing_blank::Bool)
  println(io, "    // $(f.name)")
  ret_str = f.ret === nothing ? "" : " -> $(f.ret)"

  if isempty(f.args)
    println(io, "    pub fn $(f.name)()$ret_str;")
  else
    args_inline = join(["$(a.name): $(a.rtype)" for a in f.args], ", ")
    oneliner    = "    pub fn $(f.name)($args_inline)$ret_str;"
    if length(oneliner) <= RUST_MAX_COL
      println(io, oneliner)
    else
      println(io, "    pub fn $(f.name)(")
      for a in f.args
        println(io, "        $(a.name): $(a.rtype),")
      end
      println(io, "    )$ret_str;")
    end
  end

  trailing_blank && println(io, "")
end

# ----------------------------------------------------------------
# Generate the extern "C" block
# ----------------------------------------------------------------
function gen_extern_c(io, funcs)
  println(io, "// ─────────────────────────────────────────────")
  println(io, "// extern \"C\" declarations")
  println(io, "// ─────────────────────────────────────────────")
  println(io, "")
  println(io, "#[link(name = \"uno\")]")
  println(io, "extern \"C\" {")
  for (i, f) in enumerate(funcs)
    gen_one_fn(io, f, i < length(funcs))
  end
  println(io, "}")
end

# ----------------------------------------------------------------
# Entry point
# ----------------------------------------------------------------
function main_rust()
  include_dir = joinpath(@__DIR__, "..", "..", "C")
  header_path = joinpath(include_dir, "Uno_C_API.h")

  uno_int_rust = extract_uno_int_rust(include_dir)

  args = get_default_args()
  push!(args, "-I$include_dir")

  ctx  = create_context([header_path], args)
  tu   = ctx.trans_units[1]
  root = Clang.getTranslationUnitCursor(tu)

  funcs = GC.@preserve ctx tu collect_funcs(root)

  rust_dir = joinpath(@__DIR__, "..", "..", "Rust", "src")
  out_path = joinpath(rust_dir, "ffi.rs")

  open(out_path, "w") do io
    # Inject prologue with UnoInt type at the marker
    prologue = read(joinpath(@__DIR__, "prologue_rust.rs"), String)
    marker   = "//--- UNO_INT_TYPE ---\n"
    idx      = findfirst(marker, prologue)
    idx === nothing && error("Marker '$marker' not found in prologue_rust.rs")
    print(io, prologue[1:first(idx)-1])
    println(io, "pub type uno_int = $uno_int_rust;")
    print(io, prologue[last(idx)+1:end])
    println(io, "")
    gen_extern_c(io, funcs)
  end

  println("Generated: $out_path")
end

# If run as a script with `julia wrapper_rust.jl`
if abspath(PROGRAM_FILE) == @__FILE__
  main_rust()
end
