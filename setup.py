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

#os.environ["CC"] = "icc"
#os.environ["CXX"] = "icc"

srcfiles = ['density', 'gaunt', 'linalg', 'projector', 'pseudoprojector', 'quadrature',\
			'radial', 'reader', 'sbt', 'tests', 'utils']
cfiles = [f+'.c' for f in srcfiles]
hfiles = [f+'.h' for f in srcfiles]
ext_files = cfiles + ['pawpyc.pyx']
ext_files = ['pawpyseed/core/' + f for f in ext_files]

extensions = [Extension('pawpy', ext_files,
	define_macros=[('MKL_Complex16', 'double complex'), ('MKL_Complex8', 'float complex')],
	library_dirs=['/opt/intel/compilers_and_libraries_2018/linux/mkl/lib/intel64_lin'],
	extra_link_args='-std=c11 -lmkl_rt -fopenmp -lpthread -ldl -lm -O3 -fPIC -Wall'.split(),
	runtime_library_dirs=['/opt/intel/compilers_and_libraries_2018/linux/mkl/lib/intel64_lin'],
	include_dirs=["/opt/intel/compilers_and_libraries_2018/linux/mkl/include", "pawpyseed/core", np.get_include()])]

setup(name='pawpyseed',
	version='0.2.0',
	description='Parallel C/Python package for numerical analysis of PAW DFT wavefunctions',
	long_description=long_description,
	long_description_content_type='text/markdown',
	author='Kyle Bystrom',
	author_email='kylebystrom@berkeley.edu',
	license='BSD',
	packages=['pawpyseed', 'pawpyseed.core', 'pawpyseed.analysis'],
	package_data={'pawpyseed.core': cfiles+hfiles+['Makefile', 'pawpy.so']},
	data_files=[('', ['LICENSE.txt', 'README.md'])],
	scripts=['scripts/pawpy'],
	url="https://github.com/kylebystrom/pawpyseed",
	classifiers=(
		"Programming Language :: Python :: 3",
		"License :: OSI Approved :: BSD License"
	),

	ext_modules=cythonize(extensions),
	zip_safe=False
)

if len(sys.argv) > 1 and sys.argv[1] == 'build':
	currdir = os.getcwd()
	os.chdir(os.path.join(currdir,'build/lib/pawpyseed/core'))
	if not "PAWPYCC" in os.environ:
		if subprocess.call("which icc".split()) == 0:
			os.environ["PAWPYCC"] = "icc"
		elif subprocess.call("which gcc".split()) == 0:
			os.environ["PAWPYCC"] = "gcc"
		else:
			raise PawpyBuildError("Can't find icc or gcc compiler!")

	status = subprocess.call('make pawpyinst'.split())
	if status != 0:
		raise PawpyBuildError("Can't compile pawpy.so! Check the C error output for details.")
	os.chdir(currdir)

