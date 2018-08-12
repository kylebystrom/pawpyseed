from distutils.core import setup, Extension
import os, subprocess
import shutil
import sys

class PawpyBuildError(Exception):
	pass

with open('README.md', 'r') as fh:
	long_description = fh.read()

os.environ["CC"] = "icc"
os.environ["CXX"] = "icc"

srcfiles = ['density', 'gaunt', 'linalg', 'projector', 'pseudoprojector', 'quadrature',\
			'radial', 'reader', 'sbt', 'tests', 'utils']
cfiles = [f+'.c' for f in srcfiles]
hfiles = [f+'.h' for f in srcfiles]

setup(name='pawpyseed',
	version='0.0.1',
	description='Parallelized DFT wavefunction analysis utilities',
	long_description=long_description,
	long_description_content_type='text/markdown',
	author='Kyle Bystrom',
	author_email='kylebystrom@berkeley.edu',
	packages=['pawpyseed', 'pawpyseed.core', 'pawpyseed.analysis'],
	package_data={'pawpyseed.core': cfiles+hfiles+['Makefile', 'pawpy.so']},
	scripts=['scripts/pawpy'],
	url="https://github.com/kylebystrom/pawpyseed",
	classifiers=(
		"Programming Language :: Python :: 3",
		"License :: OSI Approved :: BSD License"
	),
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

	status = subprocess.call('make pawpy'.split())
	if status != 0:
		raise PawpyBuildError("Can't compile pawpy.so! Check the C error output for details.")
	os.chdir(currdir)

