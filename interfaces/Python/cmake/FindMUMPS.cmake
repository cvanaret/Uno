# - Try to find MUMPS
# Once done, this will define
#
#  MUMPS_FOUND - system has MUMPS
#  MUMPS_INCLUDE_DIRS - the MUMPS include directories
#  MUMPS_LIBRARY - library to link to

include(LibFindMacros)

# Dependencies
# 
#libfind_package(MUMPS ScaLAPACK)

# Use pkg-config to get hints about paths
libfind_pkg_check_modules(MUMPS_PKGCONF MUMPS)

# Include dir
find_path(MUMPS_INCLUDE_DIR
  NAMES dmumps_c.h smumps_c.h
  PATHS ${MUMPS_PKGCONF_INCLUDE_DIRS}
)

# Finally the library itself
find_library(MUMPS_LIBRARY
  NAMES libdmumps
  PATHS ${MUMPS_PKGCONF_LIBRARY_DIRS}
)

find_library(MUMPS_COMMON_LIBRARY 
  NAMES libmumps_common
  PATHS ${MUMPS_PKGCONF_LIBRARY_DIRS}
)

find_library(MUMPS_PORD_LIBRARY
  NAMES libpord
  PATHS ${MUMPS_PKGCONF_LIBRARY_DIRS}
)

# libseq for sequential MUMPS
set(MUMPS_MPISEQ_DIR ${MUMPS_INCLUDE_DIR}/../libseq)

find_library(MUMPS_MPISEQ_LIBRARY
  NAMES libmpiseq
  PATHS ${MUMPS_MPISEQ_DIR}
)

# Set the include dir variables and the libraries and let libfind_process do the rest.
# NOTE: Singular variables for this library, plural for libraries this lib depends on.
set(MUMPS_PROCESS_INCLUDES MUMPS_INCLUDE_DIR)
set(MUMPS_PROCESS_LIBS "MUMPS_LIBRARY")
libfind_process(MUMPS)
