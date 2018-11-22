include(FindPackageHandleStandardArgs)

find_path(BOOST_INCLUDE_DIRS
          NAMES boost
          PATHS /usr/local /usr/local/include /usr/include /mroot/boost)

FIND_PACKAGE_HANDLE_STANDARD_ARGS(BOOST DEFAULT_MSG GSL_INCLUDE_DIRS)
