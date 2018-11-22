include(FindPackageHandleStandardArgs)


find_path(SILO_INCLUDE_DIRS
          NAMES silo.h
          PATHS /usr/local/bin/visit/2.13.1/linux-x86_64/include/silo/include)#/mroot/silo/include /usr/include /usr/local/include)

find_library(SILO_LIBRARIES
             NAMES silo siloh5
             PATHS /usr/local/bin/visit/2.13.1/linux-x86_64/lib )#/mroot/silo/lib /usr/lib /usr/local/include)

FIND_PACKAGE_HANDLE_STANDARD_ARGS(SILO DEFAULT_MSG SILO_LIBRARIES SILO_INCLUDE_DIRS)
