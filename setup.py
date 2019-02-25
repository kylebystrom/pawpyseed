# coding: utf-8

from setuptools import setup, Extension
import os, subprocess
import shutil
import sys
import codecs
from Cython.Build import cythonize
import numpy as np 

class PawpyBuildError(Exception):
	pass

with codecs.open('README.md', 'r', encoding='utf8') as fh:
	long_description = fh.read()

DEBUG = True

reqs = "numpy>=1.14,scipy>=1.0,pymatgen>=2018.2.13,sympy>=1.1.1,matplotlib>=0.2.5".split(',')

srcfiles = ['density', 'gaunt', 'linalg', 'projector', 'pseudoprojector', 'quadrature',\
			'radial', 'reader', 'sbt', 'utils']

cfiles = [f+'.c' for f in srcfiles]
#hfiles = [f+'.h' for f in srcfiles]
ext_files = cfiles
ext_files = ['pawpyseed/core/' + f for f in ext_files]
lib_dirs = ['/global/home/users/kbystrom/.local/lib']
inc_dirs = ['/global/home/users/kbystrom/.local/include', 'pawpyseed/core', np.get_include()]
if DEBUG:
	inc_dirs.append('pawpyseed/core/tests')
#if 'MKLROOT' in os.environ:
#	MKLROOT = os.environ['MKLROOT']
#	lib_dirs.append('%s/lib/intel64_lin' % MKLROOT)
#	inc_dirs.append('%s/include' % MKLROOT)
rt_lib_dirs = lib_dirs[:]
#if 'C_INCLUDE_PATH' in os.environ:
#	inc_dirs += os.environ['C_INCLUDE_PATH'].split(':')
#if 'LD_LIBRARY_PATH' in os.environ:
#	lib_dirs += os.environ['LD_LIBRARY_PATH'].split(':')
#if 'LIBRARY_PATH' in os.environ:
#	rt_lib_dirs += os.environ['LIBRARY_PATH'].split(':')
extra_args = '-std=c11 -lmkl_rt -fopenmp -fPIC -Wall'.split()
extra_args = '-std=c11 -fPIC -Wall'.split()
if not DEBUG:
	extra_args += ['-g0', '-O2']
if DEBUG:
	extra_args += ['-g']
link_args = '-lmkl_sequential -lmkl_intel_lp64 -lmkl_core -lpthread -lm -ldl'.split()
link_args = '-lmkl_rt -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -lm -ldl'.split()

extensions = [Extension('pawpyseed.core.pawpyc', ext_files + ['pawpyseed/core/pawpyc.pyx'],
	define_macros=[('MKL_Complex16', 'double complex'), ('MKL_Complex8', 'float complex')],
	library_dirs=lib_dirs,
	extra_link_args=link_args,
	extra_compile_args=extra_args,
	runtime_library_dirs=lib_dirs,
	include_dirs=inc_dirs)]
	#depends=[os.path.join('pawpyseed/core', '*.h'), os.path.join('pawpyseed/core', '*.pxd')])]
if DEBUG:
	extensions.append(Extension('pawpyseed.core.tests.testc',
		['pawpyseed/core/tests/testc.pyx', 'pawpyseed/core/tests/tests.c'] + ext_files,
		define_macros=[('MKL_Complex16', 'double complex'), ('MKL_Complex8', 'float complex')],
		library_dirs=lib_dirs,
		extra_link_args=extra_args + link_args,
		extra_compile_args=extra_args,
		runtime_library_dirs=lib_dirs,
		include_dirs=inc_dirs))
		#depends=[os.path.join('pawpyseed/core/tests', '*.h'), os.path.join('pawpyseed/core/tests', '*.pxd')]))

packages = ['pawpyseed', 'pawpyseed.core', 'pawpyseed.analysis']
if DEBUG:
	packages.append('pawpyseed.core.tests')

setup(name='pawpyseed',
	version='0.3.1',
	description='Parallel C/Python package for numerical analysis of PAW DFT wavefunctions',
	long_description=long_description,
	long_description_content_type='text/markdown',
	author='Kyle Bystrom',
	author_email='kylebystrom@berkeley.edu',
	license='BSD',
	install_requires=reqs,
	packages=packages,
	#package_data={'pawpyseed.core': cfiles+hfiles},
	data_files=[('', ['LICENSE', 'README.md'])],
	#scripts=['scripts/pawpy'],
	url="https://github.com/kylebystrom/pawpyseed",
	classifiers=(
		"Programming Language :: Python :: 3",
		"License :: OSI Approved :: BSD License"
	),
	ext_modules=cythonize(extensions,
		include_path=[os.path.join(os.path.abspath(__file__), 'pawpyseed/core')]),
	zip_safe=False
)

