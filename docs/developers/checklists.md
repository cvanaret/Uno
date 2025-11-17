# Checklists

The purpose of this page is to collate a series of checklists for commonly performed changes to the source code of Uno.

In each case, copy the checklist into the description of the pull request.

## Making a release of Uno

In preparation for an Uno release, use the following checklist. The PR should have the title "Prep for vX.Y.Z.". These steps can be done in the same commit, or separately.

````
## Pre-release

- [ ] Change the version number in `Uno::current_version()` in `uno/Uno.cpp`
- [ ] Change the version number in `CITATION.cff` (twice)
- [ ] Change the version number in `CMakeLists.txt`
- [ ] Change the version number in `interfaces/C/Uno_C_API.h`
- [ ] Update the logo in `docs/figures/logo.png`
- [ ] The commit messages in this PR do not contain `[ci skip]`

## Post-release

- `Uno_jll.jl`
  - [ ] Update the [Yggdrasil tarballs](https://github.com/JuliaPackaging/Yggdrasil/blob/master/U/Uno/build_tarballs.jl)
  - [ ] Change the `Uno_jll` version number in `interfaces/Julia/Project.toml`
  - [ ] Change the `Uno_jll` version number in `interfaces/Julia/gen/Project.toml`
- [ ] Update the logo in [the Github settings](https://github.com/cvanaret/Uno/settings)
- [ ] Update the logo in `docs/figures/logo.png`
````

## Making a release of `UnoSolver.jl`

In preparation for an `UnoSolver.jl` release, use the following checklist. The title of the PR should start with "[UnoSolver.jl] Release".

````
## Pre-release

- [ ] Change the version number in `interfaces/Julia/Project.toml`
- [ ] Tag the register bot in a commit: `@JuliaRegistrator register subdir="interfaces/Julia"
````
