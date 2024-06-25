# - Try to find BLACS
# Once done, this will define
#
#  BLACS_FOUND - system has BLACS
#  BLACS_INCLUDE_DIRS - the BLACS include directories
#  BLACS_LIBRARY - library to link to

include(LibFindMacros)

# Dependencies
# 

# Use pkg-config to get hints about paths
libfind_pkg_check_modules(BLACS_PKGCONF BLACS)


# Only check for the library since this is a link dependency
find_library(BLACS_LIBRARY
  NAMES libblacs.a libblacs-openmpi.a
  PATHS ${BLACS_PKGCONF_LIBRARY_DIRS}
)

find_library(BLACS_INIT_LIBRARY
  NAMES libblacsCinit.a libblacsCinit-openmpi.a
  PATHS ${BLACS_PKGCONF_LIBRARY_DIRS}
)

# Set the include dir variables and the libraries and let libfind_process do the rest.
# NOTE: Singular variables for this library, plural for libraries this this lib depends on.
set(BLACS_PROCESS_LIBS BLACS_LIBRARY)
libfind_process(BLACS)
