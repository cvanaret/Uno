# Checklists

The purpose of this page is to collate a series of checklists for commonly performed changes to the source code of Uno.

In each case, copy the checklist into the description of the pull request.

## Making a release

In preparation for a release, use the following checklist. These steps can be done in the same commit, or separately. The last commit should have the message "Prep for vX.Y.Z."

````
## Pre-release

 - [ ] Change the version number in `Uno::current_version()` in `uno/Uno.cpp`
 - [ ] Change the version number in `CITATION.cff`
 - [ ] The commit messages in this PR do not contain `[ci skip]`

````
