from distutils.core import setup

cargs = ['-std=c11', '-lmkl_rt', '-fopenmp', '-lpthread', '-ldl', '-lm', '-O3', '-fPIC', '-Wall']

setup(name='PAWpySeed',
	version='0.0.1',
	description='Parallelized DFT wavefunction analysis utilities',
	author='Kyle Bystrom',
	author_email='kylebystrom@berkeley.edu',
	packages=['pawpyseed'],
	ext_modules=[Extension('pawpy', ['pseudoprojector.c', 'projector.c', 'fft.c', 'quadrature.c', \
                'reader.c', 'utils.c', 'radial.c', 'tests.c', 'sbt.c'],
                library_dirs=['${MKLROOT}/lib/intel64'],
                include_dirs=['${MKLROOT}/include'],
                extra_compile_args=cargs,
                extra_link_args=cargs)]
	)