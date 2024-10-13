# - Try to find SCOTCH
# Once done, this will define
#
#  SCOTCH_FOUND - system has SCOTCH
#  SCOTCH_INCLUDE_DIRS - the SCOTCH include directories
#  SCOTCH_LIBRARY - library to link to

include(LibFindMacros)

# Dependencies
# 

# Use pkg-config to get hints about paths
libfind_pkg_check_modules(SCOTCH_PKGCONF SCOTCH)


# Only check for the library since this is a link dependency
find_library(SCOTCH_LIBRARY
  NAMES libscotch
  PATHS ${SCOTCH_PKGCONF_LIBRARY_DIRS}
)

find_library(SCOTCH_ERR_LIBRARY
  NAMES libscotcherr
  PATHS ${SCOTCH_PKGCONF_LIBRARY_DIRS}
)

find_library(SCOTCH_ERR_EXIT_LIBRARY
  NAMES libscotcherrexit
  PATHS ${SCOTCH_PKGCONF_LIBRARY_DIRS}
)

find_library(SCOTCH_METIS_LIBRARY
  NAMES libscotchmetis
  PATHS ${SCOTCH_PKGCONF_LIBRARY_DIRS}
)

# Set the include dir variables and the libraries and let libfind_process do the rest.
# NOTE: Singular variables for this library, plural for libraries this this lib depends on.
set(SCOTCH_PROCESS_LIBS SCOTCH_LIBRARY)
libfind_process(SCOTCH)
