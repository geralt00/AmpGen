#----------------------------------------------------------------
# Generated CMake target import file for configuration "Release".
#----------------------------------------------------------------

# Commands may need to know the format version.
set(CMAKE_IMPORT_FILE_VERSION 1)

# Import target "AmpGen::libAmpGen" for configuration "Release"
set_property(TARGET AmpGen::libAmpGen APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(AmpGen::libAmpGen PROPERTIES
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/lib/libAmpGen.so"
  IMPORTED_SONAME_RELEASE "libAmpGen.so"
  )

list(APPEND _IMPORT_CHECK_TARGETS AmpGen::libAmpGen )
list(APPEND _IMPORT_CHECK_FILES_FOR_AmpGen::libAmpGen "${_IMPORT_PREFIX}/lib/libAmpGen.so" )

# Commands beyond this point should not need to know the version.
set(CMAKE_IMPORT_FILE_VERSION)
