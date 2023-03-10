######################################
# ---------------------------------- #
# -------- COMPILER VERSION -------- #
# ---------------------------------- #
######################################

cmake_minimum_required(VERSION 3.10)
project(dbg_compare VERSION 0.0.0)

message(STATUS "CMAKE_CXX_COMPILER_ID: ${CMAKE_CXX_COMPILER_ID}")
message(STATUS "CMAKE_CXX_COMPILER_VERSION: ${CMAKE_CXX_COMPILER_VERSION}")

# Compiler options.
#set(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -std=c++17")
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(CMAKE_CXX_FLAGS )
#set(ASAN_OPTIONS=check_initialization_order=1)
#set(ASAN_OPTIONS=detect_leaks=1)

set(CMAKE_CXX_FLAGS "-fsanitize=address -fexceptions")  # leak checking and address sanitization
#set(CMAKE_CXX_FLAGS "-fexceptions")

#add_definitions(-ggdb3 -O0 -Wall)       # Debugging + No optimization
add_definitions(-O3 -Wall)              # Much optimization

# Definitions needed to eliminate runtime dependency
# on the boost system library.
add_definitions(-DBOOST_SYSTEM_NO_DEPRECATED)
add_definitions(-DBOOST_ERROR_CODE_HEADER_ONLY)

# Some potentially helpful messages
message(STATUS "CMAKE_BINARY_DIR: ${CMAKE_BINARY_DIR}")
message(STATUS "CMAKE_SOURCE_DIR: ${CMAKE_SOURCE_DIR}")
message(STATUS "CMAKE_CURRENT_BINARY_DIR: ${CMAKE_CURRENT_BINARY_DIR}")
message(STATUS "CMAKE_CURRENT_SOURCE_DIR: ${CMAKE_CURRENT_SOURCE_DIR}")


#########################################
# ------------------------------------- #
# -------- SOURCES AND HEADERS -------- #
# ------------------------------------- #
#########################################

# Include header files
include_directories(
        "inc"
)

# Define our shared library sources. NOT test/executables.
set(SOURCES
        src/a.cpp
        )

project(dbg_compare)
add_library(dbg_compare STATIC ${SOURCES})

# To make sure the library is named Name.so,
# get rid of the "lib" prefix.
set_target_properties(dbg_compare PROPERTIES PREFIX "")

# Eliminate an extraneous -D during compilation.
set_target_properties(dbg_compare PROPERTIES DEFINE_SYMBOL "")


##############################################
# ------------------------------------------ #
# -------- LINKING EXTERNAL LIBRARY -------- #
# ---------------- htslib ------------------ #
# ------------------------------------------ #
##############################################

message(STATUS "Fetching htslib")

# https://github.com/samtools/htslib/releases/download/1.15.1/htslib-1.15.1.tar.bz2

find_package(ZLIB REQUIRED)

# Dependencies
find_library(ZLIB NAMES z)
if(${ZLIB} STREQUAL "ZLIB-NOTFOUND")
    message(WARNING "Couldn't find the 'z' library")
endif()

find_library(BZ2LIB bz2)
if(${BZ2LIB} STREQUAL "BZ2LIB-NOTFOUND")
    message(WARNING "Couldn't find the 'bz2' library")
endif()

find_library(CURLLIB curl)
if(${CURLLIB} STREQUAL "CURLLIB-NOTFOUND")
    message(WARNING "Couldn't find the 'curl' library")
endif()

include(FetchContent)

FetchContent_Declare(
        project_htslib
        URL      https://github.com/samtools/htslib/releases/download/1.16/htslib-1.16.tar.bz2
)

FetchContent_MakeAvailable(project_htslib)

# Ensure that main library has access to primary dependencies' and secondary dependencies' headers
include_directories(
        ${CMAKE_BINARY_DIR}/_deps/project_htslib-src
        ${CMAKE_BINARY_DIR}/_deps/project_htslib-src/htslib/
        ${CMAKE_BINARY_DIR}/_deps/project_htslib-src/cram/
        ${CMAKE_BINARY_DIR}/_deps/project_htslib-src/os/
        ${CMAKE_SOURCE_DIR}/external/
)

message(STATUS "htslib_SOURCE_DIR: ${project_htslib_SOURCE_DIR}")

set(HTS_DIR _deps/project_htslib-src)

# HTSLIB configure
add_custom_target(
        BUILD_HTS
        ALL
        #        OUTPUT ${PROJECT_SOURCE_DIR}/external/htslib/config.h # //TODO: at build time ${htslib_SOURCE_DIR} is not actually populated!! is an empty string
        WORKING_DIRECTORY ${HTS_DIR}
        COMMAND pwd
        COMMAND autoreconf -i
        #        COMMAND autoheader
        COMMAND ./configure --disable-lzma --disable-libcurl --disable-bz2 --disable-s3 --disable-gcs --without-libdeflate --disable-plugins
        COMMAND $(MAKE) print-config
        COMMAND $(MAKE) prefix=${CMAKE_SOURCE_DIR}/external/htslib/ install
)

#add_library(htslib ${HTSLIB_SOURCES})

add_library(htslib STATIC IMPORTED)
set_property(TARGET htslib
        PROPERTY IMPORTED_LOCATION ${CMAKE_SOURCE_DIR}/external/htslib/lib/libhts.a)


add_dependencies(htslib BUILD_HTS)
add_dependencies(dbg_compare htslib)


##############################################
# ------------------------------------------ #
# -------- LINKING EXTERNAL LIBRARY -------- #
# ---------------- bifrost ----------------- #
# ------------------------------------------ #
##############################################

message(STATUS "Fetching bifrost")

# Dependencies
find_library(ZLIB NAMES z)
if(${ZLIB} STREQUAL "ZLIB-NOTFOUND")
    message(WARNING "Couldn't find the 'z' library")
endif()

include(FetchContent)

message(STATUS "external: ${CMAKE_SOURCE_DIR}/external/")

FetchContent_Declare(
        project_bifrost
        URL      https://github.com/pmelsted/bifrost/archive/refs/tags/v1.2.0.tar.gz
)

FetchContent_MakeAvailable(project_bifrost)

set(BIFROST_DIR _deps/project_bifrost-src)

add_custom_target(
        BUILD_BIFROST
        ALL
        WORKING_DIRECTORY ${BIFROST_DIR}
        COMMAND pwd
        COMMAND cmake -DCMAKE_INSTALL_PREFIX=${CMAKE_SOURCE_DIR}/external/bifrost/ .
        COMMAND $(MAKE) install
)

message(STATUS "bifrost_SOURCE_DIR: ${project_bifrost_SOURCE_DIR}")

# Ensure that main library has access to primary dependencies' and secondary dependencies' headers
include_directories(
        ${CMAKE_BINARY_DIR}/external/bifrost/include/bifrost/
        ${CMAKE_BINARY_DIR}/external/bifrost/include/
)

add_library(bifrost STATIC IMPORTED)

set_property(TARGET bifrost
        PROPERTY IMPORTED_LOCATION ${CMAKE_SOURCE_DIR}/external/bifrost/lib/libbifrost.a)


#add_dependencies(BUILD_BIFROST project_bifrost)
add_dependencies(bifrost BUILD_BIFROST)
add_dependencies(dbg_compare bifrost)


##############################################
# ------------------------------------------ #
# -------- LINKING EXTERNAL LIBRARY -------- #
# --------------- cuttlefish --------------- #
# ------------------------------------------ #
##############################################

#message(STATUS "Fetching cuttlefish")
#
#include(FetchContent)
#
#message(STATUS "external: ${CMAKE_SOURCE_DIR}/external/")
#
#FetchContent_Declare(
#        project_cuttlefish
#        URL             https://github.com/COMBINE-lab/cuttlefish/archive/refs/tags/v2.2.0.tar.gz
#        CMAKE_ARGS      CMAKE_INSTALL_PREFIX=${CMAKE_SOURCE_DIR}/external/cuttlefish/ .
#        INSTALL_COMMAND $(MAKE) install
#)
#
#add_custom_target(
#        BUILD_CUTTLEFISH
#        ALL
#        WORKING_DIRECTORY ${CUTTLEFISH_DIR}
#        COMMAND pwd
#        COMMAND cmake -DCMAKE_INSTALL_PREFIX=${CMAKE_SOURCE_DIR}/external/cuttlefish/ .
#        COMMAND $(MAKE) install
#)
#
#add_custom_command(
#        TARGET BUILD_CUTTLEFISH POST_BUILD
#        COMMAND ${CMAKE_COMMAND} -E copy
#        ${CUTTLEFISH_DIR}/patches/kmc_patch.diff
#        ${CMAKE_SOURCE_DIR}/patches/kmc_patch.diff)
#
#message(STATUS "CUTTLEFISH_SOURCE_DIR: ${project_cuttlefish_SOURCE_DIR}")
#
## Ensure that main library has access to primary dependencies' and secondary dependencies' headers
#include_directories(
#        ${CUTTLEFISH_DIR}/include/
#        ${CUTTLEFISH_DIR}/include/kmc_api/
#        ${CMAKE_BINARY_DIR}/external/KMC-3.2.1/
#        ${CMAKE_BINARY_DIR}/external/KMC-3.2.1/kmc_api/
#        ${CMAKE_BINARY_DIR}/external/jemalloc-5.2.1/
#        ${CMAKE_BINARY_DIR}/external/include/
#)
#
#add_library(cuttlefish2 STATIC IMPORTED)
#
#set_property(TARGET cuttlefish2
#        PROPERTY IMPORTED_LOCATION ${CMAKE_SOURCE_DIR}/external/cuttlefish/lib/libcfcore_static.a)
#
#add_library(kmc STATIC IMPORTED)
#
#set_property(TARGET kmc
#        PROPERTY IMPORTED_LOCATION ${CMAKE_SOURCE_DIR}/external/cuttlefish/lib/libkmc_core.a)
#
#add_dependencies(cuttlefish2 BUILD_CUTTLEFISH)
##add_dependencies(cuttlefish2 project_cuttlefish)
#add_dependencies(dbg_compare cuttlefish2)


############################################
# ---------------------------------------- #
# -------- Generating executables -------- #
# ---------------------------------------- #
############################################

# -------- TESTS --------

set(TESTS
        test_bifrost
        test_cuttlefish
        )


foreach(FILENAME_PREFIX ${TESTS})
    add_executable(${FILENAME_PREFIX} src/test/${FILENAME_PREFIX}.cpp)
    target_link_libraries(${FILENAME_PREFIX}
            dbg_compare
            htslib
            bifrost
#            cuttlefish2
#            kmc
            )

endforeach()


# -------- EXECUTABLES --------

set_target_properties(dbg_compare PROPERTIES LINK_FLAGS "-static" )
SET(CMAKE_EXE_LINKER_FLAGS  "${CMAKE_EXE_LINKER_FLAGS} -static-libstdc++ -static-libgcc")

set(EXECUTABLES
        compare
        )

foreach(FILENAME_PREFIX ${EXECUTABLES})
    add_executable(${FILENAME_PREFIX} src/executable/${FILENAME_PREFIX}.cpp)
    target_link_libraries(${FILENAME_PREFIX}
            dbg_compare
            htslib
            bifrost
#            cuttlefish2
#            kmc
#            -static
            )

endforeach()

# -------- final steps --------

# Where to install
set(BINARY_INSTALL_DIR ${CMAKE_CURRENT_SOURCE_DIR}/bin)

# The shared library goes to the bin directory.
install(TARGETS
        dbg_compare
        DESTINATION ${BINARY_INSTALL_DIR})
