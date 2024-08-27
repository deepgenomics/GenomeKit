cmake_minimum_required(VERSION 3.21)
project(GenomeKit)

set(CMAKE_CXX_STANDARD 20)

FILE(GLOB SOURCES src/*.cpp src/*.h)
if(DEFINED ENV{IN_CLION})
    execute_process(
            COMMAND $ENV{CONDA_PREFIX}/bin/python -c "import sys; print('{}.{}'.format(sys.version_info.major, sys.version_info.minor), end='')"
            OUTPUT_VARIABLE PYTHON_VERSION
    )

    # Don't really need this lib (handled by setup.py), but required for CLion's search etc.
    # NOTE: Set env vars IN_CLION=1;CONDA_PREFIX=$HOME/conda/envs/genomekit_dev in CLion CMake preferences.
    # include numpy, Python.h, fmt
    include_directories(
            $ENV{CONDA_PREFIX}/lib/python${PYTHON_VERSION}/site-packages/numpy/core/include/
            $ENV{CONDA_PREFIX}/include/python${PYTHON_VERSION}/
    )
    find_package(ZLIB)
    find_package(fmt)
    add_library(gkdev ${SOURCES})
    target_link_libraries(gkdev fmt::fmt ZLIB::ZLIB)
    target_compile_definitions(gkdev PRIVATE -D_FILE_OFFSET_BITS=64 -DGKPY_LIBNAME=$ENV{CONDA_PREFIX}/lib -D_WANT_MAIN=1)
endif()

# We don't want any py_* files in gk because it brings in the python runtime and causes compile errors
set(NON_PY_SOURCES ${SOURCES})
list(FILTER NON_PY_SOURCES EXCLUDE REGEX "py_.*")

add_executable(main ${NON_PY_SOURCES})
target_compile_definitions(main PRIVATE -D_FILE_OFFSET_BITS=64 -D_FILE_ABS_PATH='\"${CMAKE_CURRENT_SOURCE_DIR}/src/main.cpp\"' -D_WANT_MAIN=1)
set_target_properties(main PROPERTIES LINKER_LANGUAGE CXX)
find_package(ZLIB)
find_package(fmt)
target_link_libraries(main fmt::fmt ZLIB::ZLIB)

include(CTest)
add_test(main unittestbuild/main)
set_tests_properties(main PROPERTIES WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
