include(FindPackageHandleStandardArgs)

find_path(BLAZE_INCLUDE_DIRS
          NAMES blaze
          PATHS /usr/local/include /usr/include /mroot/blaze/include)

FIND_PACKAGE_HANDLE_STANDARD_ARGS(BLAZE DEFAULT_MSG BLAZE_INCLUDE_DIRS)
