import os
import platform
import shutil
import sys
import sysconfig
from glob import glob
from traceback import extract_stack

from setuptools import Extension
from setuptools._distutils import ccompiler
from setuptools.command.build_ext import build_ext


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
        # cannot use ccache in CXX since distutils hotpatches into LDSHARED:
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
# TODO convert to meson or scikit-build

# monkey-patch for parallel compilation
# taken from http://stackoverflow.com/questions/11013851/speeding-up-build-process-with-distutils

PARALLEL_JOBS = 4  # number of parallel compilations


def gcc_parallel_ccompile(self,
                          sources,
                          output_dir=None,
                          macros=None,
                          include_dirs=None,
                          debug=0,
                          extra_preargs=None,
                          extra_postargs=None,
                          depends=None):
    # those lines are copied from distutils.ccompiler.CCompiler directly
    macros, objects, extra_postargs, pp_opts, build = self._setup_compile(output_dir, macros, include_dirs, sources,
                                                                          depends, extra_postargs)
    cc_args = self._get_cc_args(pp_opts, debug, extra_preargs)
    # parallel code
    import multiprocessing.pool

    def _single_compile(obj):
        try:
            src, ext = build[obj]
        except KeyError:
            return
        self._compile(obj, src, ext, cc_args, extra_postargs, pp_opts)

    # convert to list, imap is evaluated on-demand
    list(multiprocessing.pool.ThreadPool(PARALLEL_JOBS).imap(_single_compile, objects))
    return objects


def windows_parallel_ccompile(self,
                              sources,
                              output_dir=None,
                              macros=None,
                              include_dirs=None,
                              debug=0,
                              extra_preargs=None,
                              extra_postargs=None,
                              depends=None):
    if not self.initialized:
        self.initialize()
    compile_info = self._setup_compile(output_dir, macros, include_dirs, sources, depends, extra_postargs)
    macros, objects, extra_postargs, pp_opts, build = compile_info

    from setuptools._distutils.errors import CompileError
    from setuptools._distutils.errors import DistutilsExecError

    compile_opts = extra_preargs or []
    compile_opts.append("/c")
    compile_opts.extend(self.compile_options_debug if debug else self.compile_options)

    def _compile_obj(obj):
        try:
            src, ext = build[obj]
        except KeyError:
            return
        if debug:
            # pass the full pathname to MSVC in debug mode,
            # this allows the debugger to find the source file
            # without asking the user to browse for it
            src = os.path.abspath(src)

        if ext in self._c_extensions:
            input_opt = "/Tc" + src
        elif ext in self._cpp_extensions:
            input_opt = "/Tp" + src
        elif ext in self._rc_extensions:
            # compile .RC to .RES file
            input_opt = src
            output_opt = "/fo" + obj
            try:
                self.spawn([self.rc] + pp_opts + [output_opt] + [input_opt])
            except DistutilsExecError as msg:
                raise CompileError(msg)
            return
        elif ext in self._mc_extensions:
            # Compile .MC to .RC file to .RES file.
            #   * "-h dir" specifies the directory for the
            #     generated include file
            #   * "-r dir" specifies the target directory of the
            #     generated RC file and the binary message resource
            #     it includes
            #
            # For now (since there are no options to change this),
            # we use the source-directory for the include file and
            # the build directory for the RC file and message
            # resources. This works at least for win32all.
            h_dir = os.path.dirname(src)
            rc_dir = os.path.dirname(obj)
            try:
                # first compile .MC to .RC and .H file
                self.spawn([self.mc] + ["-h", h_dir, "-r", rc_dir] + [src])
                base, _ = os.path.splitext(os.path.basename(src))
                rc_file = os.path.join(rc_dir, base + ".rc")
                # then compile .RC to .RES file
                self.spawn([self.rc] + ["/fo" + obj] + [rc_file])

            except DistutilsExecError as msg:
                raise CompileError(msg)
            return
        else:
            # how to handle this file?
            raise CompileError("Don't know how to compile %s to %s" % (src, obj))

        output_opt = "/Fo" + obj
        try:
            self.spawn([self.cc] + compile_opts + pp_opts + [input_opt, output_opt] + extra_postargs)
        except DistutilsExecError as msg:
            raise CompileError(msg)

    import multiprocessing.pool
    list(multiprocessing.pool.ThreadPool(PARALLEL_JOBS).imap(_compile_obj, objects))

    return objects


if sys.platform == "win32":
    import setuptools._distutils._msvccompiler
    setuptools._distutils._msvccompiler.MSVCCompiler.compile = windows_parallel_ccompile
else:
    ccompiler.CCompiler.compile = gcc_parallel_ccompile
