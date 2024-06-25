# - Try to find METIS
# Once done, this will define
#
#  METIS_FOUND - system has METIS
#  METIS_INCLUDE_DIRS - the METIS include directories
#  METIS_LIBRARY - library to link to

include(LibFindMacros)

# Dependencies
# 

# Use pkg-config to get hints about paths
libfind_pkg_check_modules(METIS_PKGCONF METIS)

# Include dir
find_path(METIS_INCLUDE_DIR
  NAMES metis.h
  PATHS ${METIS_PKGCONF_INCLUDE_DIRS}
)

# Finally the library itself
find_library(METIS_LIBRARY
  NAMES metis
  PATHS ${METIS_PKGCONF_LIBRARY_DIRS}
)

# Set the include dir variables and the libraries and let libfind_process do the rest.
# NOTE: Singular variables for this library, plural for libraries this this lib depends on.
set(METIS_PROCESS_INCLUDES METIS_INCLUDE_DIR)
set(METIS_PROCESS_LIBS METIS_LIBRARY)
libfind_process(METIS)
