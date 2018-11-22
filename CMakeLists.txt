cmake_minimum_required (VERSION 3.0)
project (TETFIT3D)

# Setup compilers
if ( APPLE )
    set (CMAKE_C_COMPILER "/myclang/bin/clang")
    set (CMAKE_CXX_COMPILER "/myclang/bin/clang++")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++1y -Wall -pedantic -O3 -I/mylibcxx/include")
    #set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++1y -Wall -pedantic -g -I/mylibcxx/include")

endif ( APPLE )

if ( CYGWIN )
   set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -ansi -O3 -DNDEBUG -mavx -std=c++1y -fdiagnostics-color=always\
        -Wno-narrowing -Wall -fpermissive")
endif ( CYGWIN )

if ( ${CMAKE_SYSTEM_NAME} MATCHES "Linux" )
   set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -ansi -O3 -DNDEBUG -mavx2 -std=c++1y -fdiagnostics-color=always\
                                            -march=native -Wno-narrowing -Wall -fpermissive")
   #~ set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -std=c++1y -DNDEBUG -mavx -fdiagnostics-color=always\
                                            #~ -march=native -Wno-narrowing -Wall -fpermissive")
   set(LINK_LIBS ${LINK_LIBS} -lpthread) # hack, must be fixed
endif( )

# Additional modules path for cmake
set (CMAKE_MODULE_PATH "${CMAKE_MODULE_PATH};${CMAKE_CURRENT_SOURCE_DIR}/cmake")

#~ # Find EIGEN
#~ find_package(EIGEN REQUIRED)
#~ include_directories("${EIGEN_INCLUDE_DIRS}")

# Find BLAZE
find_package(BLAZE REQUIRED)
include_directories("${BLAZE_INCLUDE_DIRS}")

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
include_directories(tetfit utils)

add_subdirectory(utils)

add_executable(tetfit ./main/main.cpp ${SOURCE})
target_link_libraries(tetfit "${LINK_LIBS}" "-lncurses" "-lgsl" "-lblas" "-llapack")
install (TARGETS tetfit DESTINATION bin)
#~ install (FILES ${INCLUDE_DIRS} DESTINATION include)

