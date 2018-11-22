include(FindPackageHandleStandardArgs)

find_path(EIGEN_INCLUDE_DIRS
          NAMES Eigen
          PATHS /usr/local/include/eigen3 /usr/include/eigen3 /mroot/Eigen/include/eigen3)

FIND_PACKAGE_HANDLE_STANDARD_ARGS(EIGEN DEFAULT_MSG EIGEN_INCLUDE_DIRS)
