cmake_minimum_required(VERSION 3.21)
project(GenomeKit LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

file(GLOB SOURCES CONFIGURE_DEPENDS src/*.cpp src/*.h)
file(GLOB GK_PY_EXTENSION_SOURCES CONFIGURE_DEPENDS src/*.cpp)
file(GLOB GK_HEADERS CONFIGURE_DEPENDS src/*.h)

find_package(Python COMPONENTS Interpreter Development.Module NumPy REQUIRED)
find_package(ZLIB REQUIRED)

if(NOT MSVC)
    find_program(CCACHE_PROGRAM ccache)
    if(CCACHE_PROGRAM)
        set(CMAKE_CXX_COMPILER_LAUNCHER "${CCACHE_PROGRAM}")
    endif()
endif()

Python_add_library(_cxx MODULE WITH_SOABI ${GK_PY_EXTENSION_SOURCES} ${GK_HEADERS})
target_include_directories(_cxx PRIVATE src "${CMAKE_CURRENT_BINARY_DIR}/gen")
target_link_libraries(_cxx PRIVATE Python::Module Python::NumPy ZLIB::ZLIB)
target_compile_definitions(_cxx PRIVATE GKPY_LIBNAME=genome_kit)
set_target_properties(_cxx PROPERTIES
    CXX_VISIBILITY_PRESET hidden
    VISIBILITY_INLINES_HIDDEN ON
)

if(MSVC)
    target_compile_definitions(_cxx PRIVATE _CRT_SECURE_NO_WARNINGS)
    target_compile_options(_cxx PRIVATE
        /permissive-
        /Zc:__cplusplus
        /Zc:strictStrings-
        /Zc:preprocessor
        /W3
        /EHsc
        /wd5033
        $<$<CONFIG:Debug>:/GS>
        $<$<CONFIG:Debug>:/Zi>
        $<$<CONFIG:Debug>:/Od>
        $<$<NOT:$<CONFIG:Debug>>:/GL>
        $<$<NOT:$<CONFIG:Debug>>:/Gy>
        $<$<NOT:$<CONFIG:Debug>>:/Oy>
        $<$<NOT:$<CONFIG:Debug>>:/Oi>
    )
    target_link_options(_cxx PRIVATE
        $<$<CONFIG:Debug>:/DEBUG>
        $<$<NOT:$<CONFIG:Debug>>:/LTCG>
        "/PDB:$<TARGET_FILE_DIR:_cxx>/_cxx.pdb"
    )
else()
    target_compile_definitions(_cxx PRIVATE _FILE_OFFSET_BITS=64)
    target_compile_options(_cxx PRIVATE
        -Wall
        -Wno-write-strings
        -Wno-invalid-offsetof
    )
    target_link_options(_cxx PRIVATE
        $<$<NOT:$<CONFIG:Debug>>:-Wl,-S>
        $<$<NOT:$<CONFIG:Debug>>:-Wl,-x>
    )
    if(APPLE)
        if(DEFINED ENV{MACOSX_DEPLOYMENT_TARGET})
            set(GENOMEKIT_MACOSX_DEPLOYMENT_TARGET "$ENV{MACOSX_DEPLOYMENT_TARGET}")
        else()
            set(GENOMEKIT_MACOSX_DEPLOYMENT_TARGET "10.15")
        endif()
        if(DEFINED ENV{CONDA_BUILD_SYSROOT})
            set(GENOMEKIT_MACOS_SYSROOT "$ENV{CONDA_BUILD_SYSROOT}")
        elseif(CMAKE_OSX_SYSROOT)
            set(GENOMEKIT_MACOS_SYSROOT "${CMAKE_OSX_SYSROOT}")
        else()
            execute_process(
                COMMAND xcrun --sdk macosx --show-sdk-path
                OUTPUT_VARIABLE GENOMEKIT_MACOS_SYSROOT
                OUTPUT_STRIP_TRAILING_WHITESPACE
                ERROR_QUIET
            )
        endif()
        target_compile_definitions(_cxx PRIVATE _LIBCPP_DISABLE_AVAILABILITY)
        target_compile_options(_cxx PRIVATE
            "-mmacosx-version-min=${GENOMEKIT_MACOSX_DEPLOYMENT_TARGET}"
            -stdlib=libc++
            -Wshorten-64-to-32
            -Wsign-compare
            -Wconditional-uninitialized
            -Wuninitialized
            -Wno-unknown-warning-option
        )
        target_link_options(_cxx PRIVATE
            "-mmacosx-version-min=${GENOMEKIT_MACOSX_DEPLOYMENT_TARGET}"
        )
        if(GENOMEKIT_MACOS_SYSROOT)
            target_compile_options(_cxx PRIVATE "-isysroot${GENOMEKIT_MACOS_SYSROOT}")
            target_link_options(_cxx PRIVATE "-isysroot${GENOMEKIT_MACOS_SYSROOT}")
        endif()
    endif()
endif()

install(TARGETS _cxx DESTINATION genome_kit)

if(DEFINED ENV{IN_CLION})
    execute_process(
            COMMAND $ENV{CONDA_PREFIX}/bin/python -c "import sys; print('{}.{}'.format(sys.version_info.major, sys.version_info.minor), end='')"
            OUTPUT_VARIABLE PYTHON_VERSION
    )

    # Don't really need this lib for packaging, but required for CLion's search etc.
    # NOTE: Set env vars IN_CLION=1;CONDA_PREFIX=$HOME/conda/envs/genomekit_dev in CLion CMake preferences.
    # include numpy, Python.h
    include_directories(
            $ENV{CONDA_PREFIX}/lib/python${PYTHON_VERSION}/site-packages/numpy/core/include/
            $ENV{CONDA_PREFIX}/include/python${PYTHON_VERSION}/
    )
    add_library(gkdev ${SOURCES})
    target_link_libraries(gkdev ZLIB::ZLIB)
    target_compile_definitions(gkdev PRIVATE -D_FILE_OFFSET_BITS=64 -DGKPY_LIBNAME=$ENV{CONDA_PREFIX}/lib -D_WANT_MAIN=1)
endif()

# We don't want any py_* files in gk because it brings in the python runtime and causes compile errors
set(NON_PY_SOURCES ${SOURCES})
list(FILTER NON_PY_SOURCES EXCLUDE REGEX "py_.*")

add_executable(main EXCLUDE_FROM_ALL ${NON_PY_SOURCES})
target_compile_definitions(main PRIVATE -D_FILE_OFFSET_BITS=64 -D_FILE_ABS_PATH='\"${CMAKE_CURRENT_SOURCE_DIR}/src/main.cpp\"' -D_WANT_MAIN=1)
set_target_properties(main PROPERTIES LINKER_LANGUAGE CXX)
target_link_libraries(main ZLIB::ZLIB)

include(CTest)
add_test(NAME main COMMAND main)
set_tests_properties(main PROPERTIES WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
