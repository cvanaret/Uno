#!/usr/bin/env julia
# Check that every option registered in uno/options/Options.cpp also has a
# default in DefaultOptions.cpp, is listed in uno_default.opt, and is documented
# in docs/options.md (and vice versa).

const ROOT = normpath(joinpath(@__DIR__, "..", ".."))

# Options intentionally absent from a given file.
const ALLOW_NO_DEFAULT = Set(["constraint_relaxation_strategy", "globalization_mechanism",
    "globalization_strategy", "inequality_handling_method", "option_file"])  # set by presets
const ALLOW_NO_OPTFILE = Set(["option_file"])  # command-line-only meta option

extract(path, re) = Set(m.captures[1] for line in eachline(path) for m in eachmatch(re, line))

types   = extract(joinpath(ROOT, "uno/options/Options.cpp"), r"\{\s*\"([a-zA-Z_0-9]+)\"\s*,\s*OptionType::")
default = extract(joinpath(ROOT, "uno/options/DefaultOptions.cpp"), r"options\.set_(?:string|double|integer|bool|unsigned)\(\"([a-zA-Z_0-9]+)\"")
optfile = extract(joinpath(ROOT, "uno_default.opt"), r"^#?([a-zA-Z_0-9]+)(?:\s|$)")
docs    = extract(joinpath(ROOT, "docs/options.md"), r"`([a-zA-Z_0-9]+)`")

issues = Dict{String,Vector{String}}()
flag(names, message) = foreach(name -> push!(get!(issues, name, String[]), message), names)

flag(setdiff(types, default, ALLOW_NO_DEFAULT), "missing a default in DefaultOptions.cpp")
flag(setdiff(types, optfile, ALLOW_NO_OPTFILE), "missing from uno_default.opt")
flag(setdiff(types, docs), "not documented in docs/options.md")
flag(setdiff(default, types), "set in DefaultOptions.cpp but not registered in Options.cpp")
flag(setdiff(optfile, types), "listed in uno_default.opt but not registered in Options.cpp")

if !isempty(issues)
    println("❌ Inconsistent options:\n")
    for name in sort(collect(keys(issues)))
        println("  ", name)
        foreach(message -> println("     - ", message), issues[name])
    end
    println(stderr, "\nWhen adding an option, update Options.cpp, DefaultOptions.cpp, uno_default.opt and docs/options.md.")
    exit(1)
end
println("✅ All $(length(types)) options are consistent.")
