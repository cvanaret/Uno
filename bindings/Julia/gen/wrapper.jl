# Script to parse Uno headers and generate Julia wrappers.
using Uno_jll
using Clang
using Clang.Generators
using JuliaFormatter

function main()
  # Use the header Uno_C_API.h provided by Uno_jll.jl
  # include_dir = joinpath(Uno_jll.artifact_dir, "include")
  include_dir = joinpath(@__DIR__, "..", "..", "..", "C")
  headers = [joinpath(include_dir, "Uno_C_API.h")]

  options_path = joinpath(@__DIR__, "uno.toml")
  options = load_options(options_path)

  callbacks = ["ObjectiveGradient", "Objective", "Constraints",
               "JacobianOperator", "JacobianTransposedOperator", "HessianOperator",
               "JacobianSparsity", "Jacobian", "HessianSparsity", "Hessian"]
  options["general"]["output_ignorelist"] = callbacks

  args = get_default_args()
  push!(args, "-I$include_dir")

  ctx = create_context(headers, args, options)
  build!(ctx)

  path = options["general"]["output_file_path"]
  wrappers = read(path, String)
  wrappers = replace(wrappers, r"^#.*\n?"m => "")
  for callback in callbacks
    wrappers = replace(wrappers, callback => "Ptr{Cvoid}")
  end
  write(path, wrappers)
  format_file(path, YASStyle())
  return nothing
end

# If we want to use the file as a script with `julia wrapper.jl`
if abspath(PROGRAM_FILE) == @__FILE__
  main()
end
