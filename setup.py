from distutils.core import setup, Extension
import os, subprocess
import shutil
import sys


os.environ["CC"] = "icc"
os.environ["CXX"] = "icc"

srcfiles = ['density', 'gaunt', 'mkl_linalg', 'gsl_linalg', 'projector', 'pseudoprojector', 'quadrature',\
			'radial', 'reader', 'sbt', 'tests', 'utils']
cfiles = [f+'.c' for f in srcfiles]
hfiles = [f+'.h' for f in srcfiles]

setup(name='pawpyseed',
	version='0.0.1',
	description='Parallelized DFT wavefunction analysis utilities',
	author='Kyle Bystrom',
	author_email='kylebystrom@berkeley.edu',
	packages=['pawpyseed', 'pawpyseed.core', 'pawpyseed.analysis'],
	package_data={'pawpyseed.core': cfiles+hfiles+['Makefile', 'pawpy.so']},
	scripts=['scripts/pawpy']
)

if len(sys.argv) > 1 and sys.argv[1] == 'build':
	currdir = os.getcwd()
	os.chdir(os.path.join(currdir,'build/lib/pawpyseed/core'))
	shutil.copyfile('mkl_linalg.c', 'linalg.c')
	shutil.copyfile('mkl_linalg.h', 'linalg.h')
	subprocess.call('make pawpy_icc'.split())
	os.chdir(currdir)

