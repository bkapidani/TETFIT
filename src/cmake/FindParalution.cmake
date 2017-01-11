include(FindPackageHandleStandardArgs)

find_path(PARALUTION_INCLUDE_DIRS
          NAMES paralution.hpp
          PATHS /Users/matteo/Desktop/paralution/build/inc)

      find_library(PARALUTION_LIBRARIES
          NAMES paralution
          PATHS /Users/matteo/Desktop/paralution/build/lib)

FIND_PACKAGE_HANDLE_STANDARD_ARGS(PARALUTION DEFAULT_MSG PARALUTION_LIBRARIES PARALUTION_INCLUDE_DIRS)
