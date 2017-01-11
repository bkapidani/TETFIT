find_path(PETSC_INCLUDE_DIRS
          NAMES petscksp.h
          PATHS /usr/local/include 
                /mroot/petsc/include)

find_library(PETSC_LIBRARIES
             NAMES petsc
             PATHS /usr/local/include
                   /mroot/petsc/lib)

FIND_PACKAGE_HANDLE_STANDARD_ARGS(PETSC DEFAULT_MSG PETSC_LIBRARIES PETSC_INCLUDE_DIRS)
