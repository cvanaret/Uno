# Wrapping headers and generating wrappers

This directory contains scripts `wrapper_julia.jl`, `wrapper_fortran.jl` and `wrapper_rust.jl`
that can be used to automatically generate Julia, Fortran and Rust wrappers from the C headers of Uno.
This is done using [Clang.jl](https://github.com/JuliaInterop/Clang.jl).

# Usage

Activate and instantiate the project environment in this folder
to install the dependencies `Uno_jll.jl`, `Clang.jl` and `JuliaFormatter.jl`:
```julia
shell> cd interfaces/Julia/gen
shell> julia --project
julia> ]
(gen) pkg> instantiate
```

## Julia

Regenerate the Julia wrappers (`interfaces/Julia/src/libuno.jl`) with the following commands:
```julia
julia> include("wrapper_julia.jl")
julia> main_julia()
```

If you have already instantiated the environment, you can also run:
```bash
julia --project wrapper_julia.jl
```

Note that if new constants are added to `Uno_C_API.h`, a manual update of `prologue_julia.jl` is required.

## Fortran

Regenerate the Fortran wrappers (`interfaces/Fortran/src/uno_c.f90` and `interfaces/Fortran/src/uno_fortran.f90`) with the following commands:
```julia
julia> include("wrapper_fortran.jl")
julia> main_fortran()
```

If you have already instantiated the environment, you can also run:
```bash
julia --project wrapper_fortran.jl
```

Note that if new constants or callbacks are added to `Uno_C_API.h`, a manual update of `prologue_fortran.f90` is required.

## Rust

Regenerate the Rust FFI bindings (`interfaces/Rust/src/ffi.rs`) with the following commands:
```julia
julia> include("wrapper_rust.jl")
julia> main_rust()
```

If you have already instantiated the environment, you can also run:
```bash
julia --project wrapper_rust.jl
```

Note that if new constants or callbacks are added to `Uno_C_API.h`, a manual update of `prologue_rust.rs` is required.
