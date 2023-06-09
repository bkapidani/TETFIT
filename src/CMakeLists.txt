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
	set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -std=gnu++14 -fdiagnostics-color=always -Wno-narrowing -fpermissive")
endif ( CYGWIN )

if ( ${CMAKE_SYSTEM_NAME} MATCHES "Linux" )
	set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -std=c++1y -fdiagnostics-color=always \
	-Wno-narrowing -fpermissive")
	set(LINK_LIBS ${LINK_LIBS} -lpthread) # hack, must be fixed
endif( )

# Additional modules path for cmake
set (CMAKE_MODULE_PATH "${CMAKE_MODULE_PATH};${CMAKE_CURRENT_SOURCE_DIR}/cmake")

# Find EIGEN
find_package(EIGEN REQUIRED)
include_directories("${EIGEN_INCLUDE_DIRS}")

# Find GSL
find_package(GSL REQUIRED)
include_directories("${GSL_INCLUDE_DIRS}")

# Find SPECTRA
find_package(SPECTRA REQUIRED)
include_directories("${SPECTRA_INCLUDE_DIRS}")

# Find SILO
find_package(SILO REQUIRED)
include_directories("${SILO_INCLUDE_DIRS}")
set(LINK_LIBS "${LINK_LIBS}" "${SILO_LIBRARIES}")

# Find BOOST
find_package(BOOST REQUIRED)
include_directories("${BOOST_INCLUDE_DIRS}")

# add_subdirectory(mesher)
# add_subdirectory(fdtd)

include_directories(${INCLUDE_DIRS})
include_directories(FIT DGA utils agmg)

add_subdirectory(FIT)
add_subdirectory(DGA)
add_subdirectory(utils)
add_subdirectory(agmg)

add_executable(tetfit ./DGA/main/dga_main.cpp ${SOURCE})
# add_executable(fit_uniud ./FIT/main/fit_main.cpp ${SOURCE})
target_link_libraries(tetfit "${LINK_LIBS}" "-lncurses" "-lgsl" "-lblas" "-llapack")
# target_link_libraries(fit_uniud "${LINK_LIBS}" "-lncurses" "-lgsl" "-lblas" "-llapack")
