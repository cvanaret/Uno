# - Try to find ScaLAPACK
# Once done, this will define
#
#  ScaLAPACK_FOUND - system has ScaLAPACK
#  ScaLAPACK_INCLUDE_DIRS - the ScaLAPACK include directories
#  ScaLAPACK_LIBRARY - library to link to

include(LibFindMacros)

# Dependencies
# 
#libfind_package(BLACS ScaLAPACK)
#libfind_package(ParMETIS ScaLAPACK)

# Use pkg-config to get hints about paths
libfind_pkg_check_modules(ScaLAPACK_PKGCONF ScaLAPACK)


# Only check for the library since this is a link dependency
find_library(ScaLAPACK_LIBRARY
  NAMES libscalapack libscalapack-openmpi
  PATHS ${ScaLAPACK_PKGCONF_LIBRARY_DIRS}
)

# Set the include dir variables and the libraries and let libfind_process do the rest.
# NOTE: Singular variables for this library, plural for libraries this this lib depends on.
set(ScaLAPACK_PROCESS_LIBS ScaLAPACK_LIBRARY)
libfind_process(ScaLAPACK)
