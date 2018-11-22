include(FindPackageHandleStandardArgs)

find_path(SPECTRA_INCLUDE_DIRS
          NAMES Spectra
          PATHS /usr/local/include/spectra /usr/include/spectra /mroot/spectra)

FIND_PACKAGE_HANDLE_STANDARD_ARGS(SPECTRA DEFAULT_MSG SPECTRA_INCLUDE_DIRS)
