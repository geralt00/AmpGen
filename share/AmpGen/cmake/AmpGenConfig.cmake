# - Config file for the AmpGen package
# It defines the following variables
#  AMPGEN_INCLUDE_DIRS - include directories for AmpGen
#  AMPGEN_LIBRARIES    - libraries to link against
#  AMPGEN_EXECUTABLE   - the bar executable

# Compute paths
get_filename_component(AMPGEN_CMAKE_DIR "${CMAKE_CURRENT_LIST_FILE}" PATH)
set(AMPGEN_INCLUDE_DIRS "/software/pc24403/AmpGen")
set(USE_OPENMP "TRUE")
# Our library dependencies (contains definitions for IMPORTED targets)
if(NOT TARGET AmpGen AND NOT AmpGen_BINARY_DIR)
  include("/software/pc24403/AmpGen/build/AmpGenTargets.cmake")
endif()

# These are IMPORTED targets created by FooBarTargets.cmake
set(AMPGEN_LIBRARIES AmpGen)

