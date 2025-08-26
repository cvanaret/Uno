# Wrapping headers

This directory contains a script `wrapper.jl` that can be used to automatically generate wrappers from the C headers of Uno.
This is done using [Clang.jl](https://github.com/JuliaInterop/Clang.jl).

# Usage

Activate and instantiate the project environment in this folder
to install the dependencies `Uno_jll.jl`, `Clang.jl` and `JuliaFormatter.jl`:
```julia
shell> cd Uno.jl/gen
shell> julia --project
julia> ]
(gen) pkg> instantiate
```

Then, regenerate the Julia wrappers with the following commands:
```julia
julia> include("wrapper.jl")
julia> main()
```

If you have already instantiated the environment, you can also run:
```bash
julia --project wrapper.jl
```

Note that if new constants are added to `Uno_C_API.h`, a manual update of `prologue.jl` is required.