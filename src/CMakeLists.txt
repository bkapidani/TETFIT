cmake_minimum_required (VERSION 3.0)
project (FDTD_UNIUD)

# Setup compilers
if ( APPLE )
    set (CMAKE_C_COMPILER "/myclang/bin/clang")
    set (CMAKE_CXX_COMPILER "/myclang/bin/clang++")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++1y -Wall -pedantic -O3 -I/mylibcxx/include")
    #set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++1y -Wall -pedantic -g -I/mylibcxx/include")

endif ( APPLE )

if ( CYGWIN )
	set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -std=gnu++14 -fdiagnostics-color=always \
-DENABLE_CYGWIN_MATH_FIX -fpermissive")
endif ( CYGWIN )

if ( ${CMAKE_SYSTEM_NAME} MATCHES "Linux" )
	set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -std=c++1y -fdiagnostics-color=always \
	-fpermissive")
	set(LINK_LIBS ${LINK_LIBS} -lpthread) # hack, must be fixed
endif( )

# Additional modules path for cmake
set (CMAKE_MODULE_PATH "${CMAKE_MODULE_PATH};${CMAKE_CURRENT_SOURCE_DIR}/cmake")

# Find EIGEN
find_package(EIGEN REQUIRED)
include_directories("${EIGEN_INCLUDE_DIRS}")

add_subdirectory(mesher)
add_subdirectory(fdtd)
add_subdirectory(main)

include_directories(${INCLUDE_DIRS})
include_directories(mesher fdtd)

add_executable(fdtd_uniud ${SOURCE})
target_link_libraries(fdtd_uniud "${LINK_LIBS}" "-lncurses")
