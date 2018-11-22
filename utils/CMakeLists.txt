cmake_minimum_required (VERSION 3.0)

set(SOURCE
    ${SOURCE}
    # ${CMAKE_CURRENT_SOURCE_DIR}/dga_main.cpp
   ${CMAKE_CURRENT_SOURCE_DIR}/mapped_file.cpp
   ${CMAKE_CURRENT_SOURCE_DIR}/bessel.cpp
   ${CMAKE_CURRENT_SOURCE_DIR}/timecounter.cpp
    PARENT_SCOPE)
