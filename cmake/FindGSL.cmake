include(FindPackageHandleStandardArgs)

find_path(GSL_INCLUDE_DIRS
          NAMES gsl
          PATHS /usr/local/include /usr/include /mroot/gsl)

FIND_PACKAGE_HANDLE_STANDARD_ARGS(GSL DEFAULT_MSG GSL_INCLUDE_DIRS)
