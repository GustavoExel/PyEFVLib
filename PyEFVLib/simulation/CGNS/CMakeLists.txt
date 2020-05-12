CMAKE_MINIMUM_REQUIRED (VERSION 3.14)

set (CMAKE_CXX_STANDARD 17)
# set (CMAKE_C_COMPILER $ENV{OPENMPI_DIR}/bin/mpicc)
# set (CMAKE_CXX_COMPILER $ENV{OPENMPI_DIR}/bin/mpicxx)

# include ($ENV{CMAKE_INCLUDES}/Configure.cmake.in)

# set (CMAKE_CXX_FLAGS_DEBUG ${CMAKE_CXX_FLAGS_DEBUG} "-g -DDEBUG" CACHE STRING "Debug compiler flags")
# set (CMAKE_CXX_FLAGS_RELEASE ${CMAKE_CXX_FLAGS_RELEASE} "-O3 -DNDEBUG" CACHE STRING "Release compiler flags")

# # add_definitions("-DTEST_INPUT_DIRECTORY=\"${CMAKE_SOURCE_DIR}/Zeta/Test/\"")

# set (CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/static_libs)
# set (CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/shared_libs)
set (CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/..)

##################################################################
# ABOUT THE PROJECT
##################################################################
project (CGNS)

##################################################################
# EXTERNAL LIBRARIES
##################################################################
set (CMAKE_MODULE_PATH $ENV{CMAKE_MODULES})
set (BUILD_TYPE "release")
set (CMAKE_BUILD_TYPE "Release")
set (LIBRARY_TYPE "shared")

##############
# BOOST
##############
set (BOOST_VERSION "1.70.0")
set (BOOST_ROOT $ENV{BOOST_DIR}/${BUILD_TYPE}/${LIBRARY_TYPE})
set (Components_Boost system filesystem chrono timer test_exec_monitor unit_test_framework)
set (Boost_USE_MULTITHREADED ON)
find_package (Boost ${BOOST_VERSION} EXACT COMPONENTS ${Components_Boost} REQUIRED)
if (Boost_FOUND)
    include_directories (SYSTEM ${Boost_INCLUDE_DIRS})
endif ()

##############
# CGNS
##############
set (CGNS_VERSION "3.4.0")
set (CGNS_DIR $ENV{CGNS_DIR}/${BUILD_TYPE}/${LIBRARY_TYPE})
find_package (CGNS REQUIRED)
if (CGNS_FOUND)
    include_directories (${CGNS_INCLUDE_DIR})
endif ()




add_executable(save save.cpp)
add_executable(create create.cpp)

target_link_libraries (save ${Boost_LIBRARIES})
target_link_libraries (save ${CGNS_LIBRARIES})

target_link_libraries (create ${Boost_LIBRARIES})
target_link_libraries (create ${CGNS_LIBRARIES})
