
#
# Module that checks whether GLPK is available and usable.
#
# Variables used by this module which you may want to set:
# GLPK_ROOT         Path list to search for GLPK
#
# Sets the follwing variable:
#
# GLPK_FOUND           True if GLPK available and usable.
# GLPK_INCLUDE_DIRS    Path to the GLPK include dirs.
# GLPK_LIBRARIES       Name to the GLPK library.
#

# look for header files, only at positions given by the user
find_path(GLPK_INCLUDE_DIR
  NAMES glpk.h
  PATHS ${GLPK_PREFIX} ${GLPK_ROOT}
  PATH_SUFFIXES "glpk" "include/glpk" "include" "SRC" "src"
  NO_DEFAULT_PATH
)

# look for header files, including default paths
find_path(GLPK_INCLUDE_DIR
  NAMES glpk.h
  PATH_SUFFIXES "glpk" "include/glpk" "include" "SRC" "src"
)

# look for library, only at positions given by the user
find_library(GLPK_LIBRARY
  NAMES "glpk" "libglpk"
  PATHS ${GLPK_PREFIX} ${GLPK_ROOT}
  PATH_SUFFIXES "lib" "lib32" "lib64" "libglpk"
  NO_DEFAULT_PATH
)

# look for library files, including default paths
find_library(GLPK_LIBRARY
  NAMES "glpk" "libglpk"
  PATH_SUFFIXES "lib" "lib32" "lib64" "libglpk"
)

# check version specific macros
include(CheckCSourceCompiles)
include(CMakePushCheckState)
cmake_push_check_state()

# we need if clauses here because variable is set variable-NOTFOUND
# if the searches above were not successful
# Without them CMake print errors like:
# "CMake Error: The following variables are used in this project, but they are set to NOTFOUND.
# Please set them or make sure they are set and tested correctly in the CMake files:"
#
if(GLPK_INCLUDE_DIR)
  set(CMAKE_REQUIRED_INCLUDES ${CMAKE_REQUIRED_INCLUDES} ${GLPK_INCLUDE_DIR})
endif(GLPK_INCLUDE_DIR)
if(GLPK_LIBRARY)
  set(CMAKE_REQUIRED_LIBRARIES ${CMAKE_REQUIRED_LIBRARIES} ${GLPK_LIBRARY})
endif(GLPK_LIBRARY)

# behave like a CMake module is supposed to behave
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(
  "GLPK"
  DEFAULT_MSG
  GLPK_INCLUDE_DIR
  GLPK_LIBRARY
)

mark_as_advanced(GLPK_INCLUDE_DIR GLPK_LIBRARY)

# if both headers and library are found, store results
if(GLPK_FOUND)
  set(GLPK_INCLUDE_DIRS ${GLPK_INCLUDE_DIR})
  set(GLPK_LIBRARIES    ${GLPK_LIBRARY})
  # log result
  file(APPEND ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeOutput.log
    "Determining location of GLPK succeeded:\n"
    "Include directory: ${GLPK_INCLUDE_DIRS}\n"
    "Library directory: ${GLPK_LIBRARIES}\n\n")
  set(GLPK_DUNE_COMPILE_FLAGS "-I${GLPK_INCLUDE_DIRS}"
    CACHE STRING "Compile flags used by DUNE when compiling GLPK programs")
  set(GLPK_DUNE_LIBRARIES ${GLPK_LIBRARIES}
    CACHE STRING "Libraries used by DUNE when linking GLPK programs")
else(GLPK_FOUND)
  # log errornous result
  file(APPEND ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeError.log
    "Determining location of GLPK failed:\n"
    "Include directory: ${GLPK_INCLUDE_DIRS}\n"
    "Library directory: ${GLPK_LIBRARIES}\n\n")
endif(GLPK_FOUND)

# set HAVE_GLPK for config.h
set(HAVE_GLPK ${GLPK_FOUND})

# register all GLPK related flags
if(GLPK_FOUND)
  dune_register_package_flags(COMPILE_DEFINITIONS "ENABLE_GLPK=1"
                              LIBRARIES "${GLPK_DUNE_LIBRARIES}"
                              INCLUDE_DIRS "${GLPK_INCLUDE_DIRS}")
endif()
