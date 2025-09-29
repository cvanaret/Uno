# Checklists

The purpose of this page is to collate a series of checklists for commonly performed changes to the source code of Uno.

In each case, copy the checklist into the description of the pull request.

## Making a release

In preparation for a release, use the following checklist. These steps can be done in the same commit, or separately. The last commit should have the message "Prep for vX.Y.Z."

````
## Pre-release

 - [ ] Change the version number in `Uno::current_version()` in `uno/Uno.cpp`
 - [ ] Change the version number in `CITATION.cff` (twice)
 - [ ] Change the version number in `CMakeLists.txt`
 - [ ] Change the version number in `bindings/C/Uno_C_API.h`
 - [ ] Change the version number in `bindings/Julia/Uno.jl/Project.toml`
 - [ ] Change the version number in `bindings/Julia/Uno.jl/gen/Project.toml`
 - [ ] The commit messages in this PR do not contain `[ci skip]`

````