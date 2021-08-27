import codecs
import configparser
import os
import sys

import numpy as np
from Cython.Build import cythonize
from setuptools import Extension, setup


class PawpyBuildError(Exception):
    pass


with codecs.open("README.md", "r", encoding="utf8") as fh:
    long_description = fh.read()

DEBUG = True

srcfiles = [
    "density",
    "gaunt",
    "linalg",
    "projector",
    "pseudoprojector",
    "quadrature",
    "radial",
    "reader",
    "sbt",
    "utils",
    "momentum",
]

# SEARCH FOR AND READ CONFIGURATION FILE
config = configparser.ConfigParser()
user_cfg_file = os.path.expanduser("~/.pawpyseed-site.cfg")
if os.path.isfile("site.cfg"):
    config.read("site.cfg")
elif os.path.isfile(user_cfg_file):
    config.read(user_cfg_file)
else:
    config.read_file(open("site.cfg.default"))

# SET COMPILER AND LINKER IF SET IN CONFIG FILE
if "compiler_name" in config["compiler"]:
    os.environ["CC"] = config["compiler"]["compiler_name"]
    if not "linker_name" in config["compiler"]:
        os.environ["LDSHARED"] = config["compiler"]["compiler_name"] + " -shared"
if "linker_name" in config["compiler"]:
    os.environ["LDSHARED"] = config["compiler"]["linker_name"]

# SET PARALLELIZATION AND INTERFACE OPTIONS
sdl = config["mkl"].getboolean("sdl")
omp_loops = config["threading"].getboolean("omp_loops")
threaded_mkl = config["threading"].getboolean("threaded_mkl")

if sys.platform == "darwin":
    # platform_link_args = ['-lmkl_avx512']
    platform_link_args = []
    sdl_platform_link_args = []
    os.environ[
        "CPATH"
    ] = "/Applications/Xcode.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX.sdk/usr/include"
else:
    platform_link_args = ["-Wl,--no-as-needed", "-lmkl_def"]
    sdl_platform_link_args = ["-Wl,--no-as-needed"]


if sdl:
    link_args = sdl_platform_link_args + "-lmkl_rt -liomp5 -lpthread -lm -ldl".split()
else:
    interfacelib = "-lmkl_intel_lp64"
    if threaded_mkl:
        threadlib = "-lmkl_intel_thread"
        omplib = "-liomp5"
    else:
        threadlib = "-lmkl_sequential"
        omplib = ""
    link_args = "{} {} -lmkl_core {} -lpthread -lm -ldl".format(
        interfacelib,
        threadlib,
        omplib,
    )
    link_args = platform_link_args + link_args.split()

# SET OTHER COMPILER ARGS
extra_args = "-std=c11 -fPIC -Wall".split()
if omp_loops:
    extra_args.append("-fopenmp")

# ADD ADDITIONAL MKL LIBRARIES IF FOUND IN CONFIG
lib_dirs = []
inc_dirs = ["pawpyseed/core", np.get_include()]
if "root" in config["mkl"]:
    root_dirs = config["mkl"]["root"].split(":")
    for r in root_dirs:
        lib_dirs.append(os.path.join(r, "lib/intel64"))
        lib_dirs.append(os.path.join(r, "lib"))
        inc_dirs.append(os.path.join(r, "include"))
if "extra_libs" in config["compiler"]:
    extra_dirs = config["compiler"]["extra_libs"].split(":")
    for d in extra_dirs:
        lib_dirs.append(d)

# SET UP COMPILE/LINK ARGS AND THREADING
cfiles = [f + ".c" for f in srcfiles]
ext_files = cfiles
ext_files = ["pawpyseed/core/" + f for f in ext_files]
if DEBUG:
    inc_dirs.append("pawpyseed/core/tests")
rt_lib_dirs = lib_dirs[:]

if not DEBUG:
    extra_args += ["-g0", "-O2"]
else:
    extra_args += ["-g"]

extensions = [
    Extension(
        "pawpyseed.core.pawpyc",
        ext_files + ["pawpyseed/core/pawpyc.pyx"],
        define_macros=[
            ("MKL_Complex16", "double complex"),
            ("MKL_Complex8", "float complex"),
        ],
        library_dirs=lib_dirs,
        extra_link_args=extra_args + link_args,
        extra_compile_args=extra_args,
        runtime_library_dirs=rt_lib_dirs,
        include_dirs=inc_dirs,
    )
]

if DEBUG:
    extensions.append(
        Extension(
            "pawpyseed.core.tests.testc",
            ["pawpyseed/core/tests/testc.pyx", "pawpyseed/core/tests/tests.c"]
            + ext_files,
            define_macros=[
                ("MKL_Complex16", "double complex"),
                ("MKL_Complex8", "float complex"),
            ],
            library_dirs=lib_dirs,
            extra_link_args=extra_args + link_args,
            extra_compile_args=extra_args,
            runtime_library_dirs=rt_lib_dirs,
            include_dirs=inc_dirs,
        )
    )

packages = ["pawpyseed", "pawpyseed.core", "pawpyseed.analysis"]
if DEBUG:
    packages.append("pawpyseed.core.tests")

setup(
    name="pawpyseed",
    version="0.7.0",
    description="Parallel C/Python package for numerical analysis of PAW DFT wavefunctions",
    long_description=long_description,
    long_description_content_type="text/markdown",
    author="Kyle Bystrom",
    author_email="kylebystrom@gmail.com",
    maintainer="Kyle Bystrom",
    maintainer_email="kylebystrom@gmail.com",
    license="BSD",
    setup_requires=["mkl-devel", "numpy>1.14", "Cython>=0.29.21"],
    install_requires=[
        "mkl-devel",
        "numpy>=1.14",
        "scipy>=1.0",
        "pymatgen>=2018.2.13",
        "sympy>=1.1.1",
        "matplotlib>=0.2.5",
    ],
    python_requires=">=3.6",
    packages=packages,
    data_files=[("", ["LICENSE", "README.md"])],
    url="https://github.com/kylebystrom/pawpyseed",
    classifiers=(
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: BSD License",
    ),
    ext_modules=cythonize(
        extensions,
        include_path=[os.path.join(os.path.abspath(__file__), "pawpyseed/core")],
    ),
    zip_safe=False,
)
