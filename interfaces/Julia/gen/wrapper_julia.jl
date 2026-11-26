# Script to parse Uno headers and generate Julia wrappers.
using Clang
using Clang.Generators
using JuliaFormatter

function main_julia()
  include_dir = joinpath(@__DIR__, "..", "..", "C")
  headers = [joinpath(include_dir, "Uno_C_API.h")]

  options_path = joinpath(@__DIR__, "uno.toml")
  options = load_options(options_path)

  callbacks = ["uno_objective_callback",
               "uno_objective_gradient_callback",
               "uno_constraints_callback",
               "uno_constraints_jacobian_callback",
               "uno_lagrangian_hessian_callback",
               "uno_constraints_jacobian_operator_callback",
               "uno_constraints_jacobian_transposed_operator_callback",
               "uno_lagrangian_hessian_operator_callback",
               "uno_notify_acceptable_iterate_callback",
               "uno_termination_callback",
               "uno_logger_stream_callback"]

  # the opaque handles uno_model and uno_solver are passed as Ptr{Cvoid}
  handles = ["uno_model", "uno_solver"]
  # the constants (enumerators and version macros) are defined in the prologue, as Cint
  constants = ["uno_optimization_sense", "uno_multiplier_sign_convention", "uno_base_indexing", "uno_triangular_part",
               "uno_option_type", "uno_options_status", "uno_model_status", "uno_solver_status", "uno_result_status",
               "uno_optimize_status", "uno_optimization_status", "uno_solution_status",
               "UNO_VERSION_MAJOR", "UNO_VERSION_MINOR", "UNO_VERSION_PATCH"]
  options["general"]["output_ignorelist"] = [callbacks; handles; constants; "uno_int"]

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
  for handle in handles
    wrappers = replace(wrappers, "Ptr{$handle}" => "Ptr{Cvoid}")
  end
  wrappers = replace(wrappers, "uno_int" => "Int32")
  wrappers = "# Copyright (c) 2026 Alexis Montoison and Charlie Vanaret\n# Licensed under the MIT license. See LICENSE file in the project directory for details.\n\n" * wrappers
  write(path, wrappers)
  format_file(path, YASStyle())
  return nothing
end

# If we want to use the file as a script with `julia wrapper_julia.jl`
if abspath(PROGRAM_FILE) == @__FILE__
  main_julia()
end