CMAKE_MINIMUM_REQUIRED (VERSION 3.14)

set (CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/static_libs)
set (CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/shared_libs)
set (CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

if (UNIX)
    set (CMAKE_BUILD_TYPE "Release")
    string (TOLOWER ${CMAKE_BUILD_TYPE} BUILD_TYPE)

    set (LIBRARY_TYPE "shared")

endif ()

project (CGNS)
set (CMAKE_MODULE_PATH $ENV{CMAKE_MODULES})

set (BOOST_ROOT $ENV{BOOST_DIR}/${BUILD_TYPE}/${LIBRARY_TYPE})

set (Components_Boost system filesystem unit_test_framework test_exec_monitor)
set (Boost_USE_MULTITHREADED ON)

find_package (Boost COMPONENTS ${Components_Boost} REQUIRED)
include_directories (SYSTEM ${Boost_INCLUDE_DIRS})

set (CGNS_DIR $ENV{CGNS_DIR}/${BUILD_TYPE}/${LIBRARY_TYPE})

find_package (CGNS COMPONENTS cgns_${LIBRARY_TYPE} REQUIRED)
include_directories (${CGNS_INCLUDE_DIR})

macro (add_executable _target)
    _add_executable (${_target} ${ARGN})
    target_link_libraries (${_target} ${Boost_LIBRARIES})
    target_link_libraries (${_target} ${CGNS_LIBRARIES})
endmacro ()


##################################################################
# MYPROJECT
################################################################## 

# add_library(CgnsOpener include/CgnsOpener.hpp source/CgnsOpener.cpp)
# target_include_directory(CgnsOpener ${CMAKE_SOURCE_DIR}/include)

add_executable(CgnsSaver main.cpp source/CgnsOpener.cpp)
# target_link_library(main PUBLIC CgnsOpener)
