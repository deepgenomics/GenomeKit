# Copyright (C) 2016-2023 Deep Genomics Inc. All Rights Reserved.
import os
import platform
import sys
import sysconfig
from glob import glob
from pathlib import Path
from traceback import extract_stack

import shutil
from setuptools import Extension
from setuptools import setup, find_packages
from setuptools.command.build_ext import build_ext
from setuptools.command.egg_info import egg_info


COPYRIGHT_FILE = "COPYRIGHT.txt"
LICENSE_FILE = "LICENSE"

tests_require = [
    "twobitreader>=3.1",
]

version = "7.1.0"

# See https://stackoverflow.com/questions/9977889/how-to-include-license-file-in-setup-py-script/66443941#66443941
class egg_info_ex(egg_info):
    """Includes license file into `.egg-info` folder."""

    def run(self):
        # don't duplicate license into `.egg-info` when building a distribution
        if not self.distribution.have_run.get('install', True):
            # `install` command is in progress, copy license
            self.mkpath(self.egg_info)
            self.copy_file(COPYRIGHT_FILE, self.egg_info)
            self.copy_file(LICENSE_FILE, self.egg_info)

        egg_info.run(self)



libname = "genome_kit"

##############################################################
# Define the C Extension build config
##############################################################

debug = False  # Set to True for C++ debug symbols and extra memory/bounds checks
debug_info = debug
debug_objects = False  # Enable printing of Python C++ instances being constructed/destructed
toolset = "msvc" if platform.system() == "Windows" else "gcc"

gen_dir = os.path.join('build', 'gen')
include_dirs = [
    sys.prefix + "/include",
    gen_dir,
    ]
library_dirs = [sys.prefix + '/lib']
runtime_library_dirs = []

if not [x for x in extract_stack(limit=20) if 'load_setup_py_data' in x[2]]:
    # don't load numpy if we just need the version to fill the conda-build meta.yaml template
    import numpy as np
    include_dirs.append(np.get_include())

sources = glob("src/*.cpp")
headers = glob("src/*.h")
libraries = []

define_macros = [
    ("GKPY_LIBNAME", libname),
]

if debug_objects:
    define_macros += [
        ("GKPY_TRACE_MEM", None),  # ctor/dtor notifications from C objects
    ]

extra_compile_args = []
extra_link_args = []

if toolset == "gcc":
    # Using multiprocessing and ccache massively speeds up incremental builds
    ccache = shutil.which('ccache')
    if ccache:
        os.environ["CC"] = "{} {}".format(ccache, os.environ.get("CC", "gcc"))
        # TODO cannot use ccache in CXX since distutils hotpatches into LDSHARED:
        # https://github.com/python/cpython/blob/069306312addf87252e2dbf250fc7632fc8b7da3/Lib/distutils/unixccompiler.py#L191
        os.environ["LDSHARED"] = sysconfig.get_config_var("LDSHARED")
    else:
        print("WARNING: did not find ccache installed; files will be built from scratch every time")

    define_macros += [
        ("_FILE_OFFSET_BITS", 64),
    ]
    if debug:
        define_macros += [
            ("GK_DEBUG", None),  # Enable debug assertions etc
        ]

    libraries += [
        # python is not linked for conda's python
        # https://github.com/ContinuumIO/anaconda-issues/issues/9078#issuecomment-378321357
        "z",
    ]

    # GCC flags common to both debug and release modes
    extra_compile_args += [
        "-std=c++20",
        "-fvisibility=hidden",  # reduce symbols for code size/load times
        "-fvisibility-inlines-hidden",
        "-Wall",
        "-Wno-write-strings",  # char* in Python API
        "-Wno-invalid-offsetof",  # offsetof non-POD GK types for Python API
    ]

    if debug:
        opt_args = [
            "-O0",
            "-UNDEBUG",
        ]
    else:
        opt_args = [
            "-O3",
        ]

    if debug_info:
        opt_args += [
            "-g3",
        ]
    else:
        extra_link_args += [
            "-Wl,-S",
            "-Wl,-x",
        ]

    if platform.system() == "Darwin":
        # >=10.15 required for std::filesystem::remove
        osx_sdk = "-mmacosx-version-min={}".format(os.environ.get("MACOSX_DEPLOYMENT_TARGET", "10.15"))

        extra_compile_args += [
            osx_sdk,
            "-isysroot{}".format(os.environ.get("CONDA_BUILD_SYSROOT", "/Library/Developer/CommandLineTools/SDKs/MacOSX.sdk")),
            "-stdlib=libc++",
            "-Wshorten-64-to-32",  # catch implicit truncation as per MSVC
            "-Wsign-compare",  # match MSVC(/W3)/gcc(-Wall)
            "-Wconditional-uninitialized",  # gcc does better here but enable for safety
            "-Wuninitialized",
            "-Wno-unknown-warning-option",
        ]

        extra_link_args += [
            osx_sdk,
            "-isysroot{}".format(os.environ.get("CONDA_BUILD_SYSROOT", "/Library/Developer/CommandLineTools/SDKs/MacOSX.sdk")),
        ]
        define_macros += [
            # https://conda-forge.org/docs/maintainer/knowledge_base/#newer-c-features-with-old-sdk
            ("_LIBCPP_DISABLE_AVAILABILITY", None),
        ]

    extra_compile_args += opt_args
    extra_link_args += opt_args  # required for LTO


elif toolset == "msvc":

    os.environ["DISTUTILS_USE_SDK"] = "1"
    os.environ["MSSdk"] = "1"

    condalib_dir = sys.prefix + "/Library"
    condalib_inc = condalib_dir + "/include"
    condalib_lib = condalib_dir + "/lib"
    include_dirs.append(condalib_inc)
    library_dirs.append(condalib_lib)

    define_macros += [
        ("_CRT_SECURE_NO_WARNINGS", None),
    ]

    libraries += [
        "zlib",
    ]

    # VC flags common to both debug and release modes
    extra_compile_args += [
        "/std:c++20",
        "/permissive-",
        "/Zc:__cplusplus",
        "/Zc:strictStrings-",  # don"t let strings be written to by default
        # compatibility with __VA_OPT__
        # see https://devblogs.microsoft.com/cppblog/announcing-full-support-for-a-c-c-conformant-preprocessor-in-msvc/
        "/Zc:preprocessor",
        "/W3",  # Warning level 3
        "/EHsc",  # Enable C++ and structured exception handling (e.g. catch access violations)
        "/wd5033",  # Python usage of deprecated register keyword
    ]

    extra_link_args += [
        "/PDB:%s\\_cxx.pdb" % libname,  # Put it right beside the .pyd/.so file in "develop" mode
    ]

    if debug:
        # do NOT define _DEBUG; need non-debug runtime to match NDEBUG Python distribution
        # define_macros += [
        #      ("_DEBUG", None),
        # ]
        extra_compile_args += [
            "/GS",  # Enable buffer overrun checks
            "/Zi",  # Enable debug information .pdb
            "/Od",  # Disable optimizations
        ]
        extra_link_args += [
            "/DEBUG",
        ]
    else:
        extra_compile_args += [
            "/GL",  # Enable whole-program optimization
            "/Gy",  # Enable function-level linking
            "/Oy",  # Omit frame pointers
            "/Oi",  # Enable intrinsics
            #"/Zi",  # Enable debug information .pdb
        ]
        extra_link_args += [
            "/LTCG",  # Enable link-time code generation
            #"/DEBUG",
        ]


class NoCWarningsBuildExt(build_ext):
    def build_extensions(self):
        for x in ["-Wstrict-prototypes"]:
            try:
                self.compiler.compiler_so.remove(x)
            except (AttributeError, ValueError):
                continue
        build_ext.build_extensions(self)


extension = Extension(
    libname + "._cxx",
    sources=sources,
    depends=headers,
    include_dirs=include_dirs,
    define_macros=define_macros,
    library_dirs=library_dirs,
    libraries=libraries,
    runtime_library_dirs=runtime_library_dirs,
    extra_compile_args=extra_compile_args,
    extra_link_args=extra_link_args,
    )

##############################################################


if __name__ == "__main__":
    setup(
        author="Deep Genomics",
        author_email="info@deepgenomics.com",
        python_requires=">=3.9, <4",
        classifiers=[
            "Development Status :: 5 - Production/Stable",
            "License :: OSI Approved :: Apache Software License",
            "Programming Language :: Python :: 3.9",
            "Programming Language :: Python :: 3.10",
            "Programming Language :: Python :: 3.11",
            "Programming Language :: Python :: 3.12",
        ],
        description="GenomeKit is a Python library for fast and easy access to genomic resources such as sequence, data tracks, and annotations.",
        long_description=(Path(__file__).parent / "README.md").read_text(),
        long_description_content_type='text/markdown',
        install_requires=[
            "appdirs",
            "numpy",
            "google-cloud-storage",
            "boto3",
            "tqdm",
            "importlib-metadata",
            "typing-extensions",
        ],
        license="Apache License 2.0",
        license_files=(COPYRIGHT_FILE, LICENSE_FILE,),
        name="genomekit",
        packages=find_packages(include=["genome_kit"]),
        project_urls={
            "Documentation": "https://deepgenomics.github.io/GenomeKit"
        },
        cmdclass={
            'build_ext': NoCWarningsBuildExt,
            'egg_info': egg_info_ex
        },
        ext_modules=[extension],
        test_suite="tests",
        tests_require=tests_require,
        url=f"https://github.com/deepgenomics/GenomeKit",
        version=version,
        zip_safe=False,
    )
