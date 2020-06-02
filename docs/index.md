# pawpyseed

<https://kylebystrom.github.io/pawpyseed/> <br/>

PAWpySeed is a parallelized Python and C tool for reading and
analyzing the optimized band structure and wave functions
of VASP DFT calculations. The code is written for the PAW
formalism developed by P.E. Blochl and implemented
in VASP.

## Installation

TO compile the C code, pawpyseed needs to link with the Intel Math Kernel Library
(MKL). You can customize how this is done via the config file (see "The customizable way").
However, for most users, the easy way described below is adequate.

### The easy way

First, install `mkl-devel` in your local installation:

```
pip install mkl-devel --user
```

Next, find your Python installation directory (call it `INST_DIR`) and add the `include` and `lib` to your
C compiler paths to link with MKL:
```
export C_INCLUDE_PATH=$INST_DIR/include:$C_INCLUDE_PATH
export LIBRARY_PATH=$INST_DIR/lib:$LIBRARY_PATH
export LD_LIBRARY_PATH=$INST_DIR/lib:$LD_LIBRARY_PATH
```

Then run
```
python setup.py build
python setup.py install
```

OR with `pip`.

```
pip install pawpyseed
```

This has been tested on Scientific Linux 7 and Linux Mint 18,
but should work for systems that have the appropriate
packages and environment variables defined as described below.

#### Simple site.cfg approach

After installing MKL locally with `pip`,
copy `site.cfg.default` from the pawpyseed repository to `~/.pawpyseed-site.cfg`.
Open it and set `root=<your home directory>/.local` under the `[mkl]` heading
and uncomment it by removing the `#`.

IF you get linking issues at runtime using this method, try setting
`sdl=True` in the config file and then setting the environment
variable `MKL_THREADING_LAYER=sequential`.

#### Intel easy way

If you have `icc` and prefer to use it to compile pawpyseed's C libraries,
set `compiler=icc` in `~/.pawpyseed-site.cfg` and then
set `root` to your MKL installation directory. Set `sdl=True`. You MKL distribution must have also
been compiled with the Intel compiler. Run the `setup.py` script or `pip install pawpyseed`.

### The customizable way

To make your own config file,
make a file in your home directory `~` named `.pawpyseed-site.cfg`. This file allows you
to customize configuration settings for pawpyseed. It is read using the Python
configparser module <https://docs.python.org/3/library/configparser.html>. If you don't
want to learn about the configuration file options, skip to "The Easy Way" to get this done with quickly.
The options (with their defaults) are as follows:

```
[compiler]
# Name of the compiler to be used. By default, let's the setup script choose this.
compiler_name = <no default> # Examples: icc or gcc
# Name of the linker to be used. By default, let's the setup script choose this.
# NOTE: Don't forget the -shared tag for icc or gcc!
# NOTE: If compiler_name is set and linker_name is not,
#       linker_name is set to compiler_name with the -shared tag appended.
linker_name = <no default> # Examples: icc -shared or gcc -shared

[mkl]
# MKL installation directory. <root>/lib or <root>/lib/intel64 must contain
# the MKL shared object libraries, while <root>/include must containn the MKL headers.
# Many systems have the environment variable
# MKLROOT set for this, so you can set this to the output of:
#       echo $MKLROOT
# If you pip install mkl-devel, this goes to /usr or /usr/local.
# If you pip install mkl-devel --user, this goes to ~/.local
# (or the equivalents for your system).
root = <no default>
# Don't change this. Might be used in the future but currently has no effect.
interface32 = True
# Whether to compile with the single dynamic library libmkl_rt.so.
# If you are having installation/performance problems, try setting sdl=True,
# pip installing mkl-devel with the --user option, set root=~/.local,
# and set the environment variable MKL_THREADING_LAYER=sequential
# when using pawpyseed.
sdl = False

[threading]
# Whether to use omp_loops. If False, pawpyseed doesn't do any parallelization
# on its own, though its calls to MKL, BLAS, LAPACK, etc might be
# threaded by MKL.
omp_loops = True
# Whether to use threading for MKL. OVERRIDDEN if sdl=True, in which case
# MKL is threaded by default and you can run sequentially by
# setting the environment variable MKL_THREADING_LAYER=sequential
# when using pawpyseed.
# NOTE: ONLY SET THIS TRUE IF COMPILING WITH icc (OR gcc with omp_loops=False),
# and only if you are very concerned about speed. This is because gcc
# cannot do nested parallelism with MKL
threaded_mkl = False

# NOTE: Do not link to icc-compiled MKL libraries when compiling the C
# extension with gcc or vice versa. Use like type compilers.
```

### Dependencies

All dependencies indicate the minimum version tested.
PAWpySeed might work fine with earlier versions, but
use of older versions will not be officially supported.

Python requirements:
```
Python>=3.5
numpy>=1.14
scipy>=1.0
pymatgen>=2018.2.13
sympy>=1.1.1
matplotlib>=0.2.5
```

C requirements:
```
icc >= 16.0.4 OR gcc >= 4.8.5
Intel Math Kernel Library >= 11.3.4
```
If you don't want to `pip install` Intel MKL,
it is available for free installation on a variety of platforms.
Most computing clusters will have Intel MKL, and you can install it
on your desktop (or any system to which you have root access) by following
the relevant instructions at the following URL:
<https://software.intel.com/en-us/articles/free-ipsxe-tools-and-libraries>.

## Theory and Input

### PAW

The projector augmented wave (PAW) method is a technique
used in plane wave density functional theory to simplify
the description of the wavefunctions near the nuclei
of a system. The strong Coulombic forces near an atomic
nucleus creates quickly oscillating wavefunctions that are
not well described by plane waves without prohibitively
large basis sets, so a "pseudopotential" is introduced
near the atomic nuclei which results in smooth 
"pseudo wavefunctions" well described by plane waves. The
full wavefunctions can be recovered by a linear transform
of the pseudowavefunctions. The PAW method requires
three sets of functions: projector functions, onto which
pseudowavefunctions are projected to probe their character;
full partial waves, which describe atomic valence states
derived from the true potential; and pseudo partial waves,
which are derived from the full partial waves and
pseudopotential.

### Files

The projector functions and partial waves are unique
to each element and stored in the POTCAR file
used in a VASP calculation. The pseudowavefunction
is the part of the wavefunction optimized during a DFT
calculation and is stored in the WAVECAR output file
in VASP. PAWpySeed parses both files to retrieve
all parts of the full Kohn Sham wavefunctions.

## The Code

The main purpose of PAWpySeed is to evaluate overlap
operators between Kohn-Sham wavefunctions from different
structures, which is not done by standard plane-wave DFT codes.
Such functionality can be useful for analyzing the composition
of defect levels in solids, which is main application for which
the code is currently focused.

### Implementation

* Python Interface
* Computationally intensive tasks in C
* Parallelized with openmp

### Current Functionality

* Read pseudowavefunctions
* Read projectors and partial waves from VASP POTCAR
* Evaluate overlap operators between bands,
including when bands belong to different structures
with the same lattice
* Project point defect levels onto bulk valence
and conduction bands
* Convenient pycdt interface
* Perturbation-extrapolation correction for point defect calculations
* Read noncollinear pseudo wavefunctions and construct all-electron wavefunctions (no overlap operator evaluation for noncollinear data)
* Calculate plane-wave matrix elements between bands.

## Acknowledgments

The code in PAWpySeed is based on a several algorithms and codes, which are enumerated
here.

1. **PAW**:
The PAW method was developed by P. E. Blochl in 1994. His paper deriving the method
was helpful to me in deriving the extensions to the formalism needed to develop
this code.
    * P. E. Blochl. Projector augmented-wave method. Phys. Rev. B, 50:17953, 1994.
2. **VASP**:
PAWpySeed is primarily built to read and process the output of VASP calculations.
PAWpySeed reads PAW wavefunctions and calculate overlap operators using algorithms
derived from VASP and other plane-wave codes, so the following citations are necessary.
The last citation is specifically for the PAW method and potentials. See the VASP website
at <http://community.hartree.stfc.ac.uk/wiki/site/admin/vasp.html> for information
on citing specific functionals.
    * G. Kresse and J. Hafner. Ab initio molecular dynamics for liquid metals. Phys. Rev. B, 47:558, 1993.
    * G. Kresse and J. Hafner. Ab initio molecular-dynamics simulation of the liquid-metal-amorphous-semiconductor transition in germanium. Phys. Rev. B, 49:14251, 1994.
    * G. Kresse and J. Furthmüller. Efficiency of ab-initio total energy calculations for metals and semiconductors using a plane-wave basis set. Comput. Mat. Sci., 6:15, 1996.
    * G. Kresse and J. Furthmüller. Efficient iterative schemes for ab initio total-energy calculations using a plane-wave basis set. Phys. Rev. B, 54:11169, 1996.
    * G. Kresse and D. Joubert. From ultrasoft pseudopotentials to the projector augmented-wave method. Phys. Rev. B, 59:1758, 1999.
3. **NUMSBT**:
NUMSBT is a code written by J. D. Talman, which implements an algorithm that
calculates the spherical Bessel transform (SBT) in O(NlogN) time.
PAWpySeed
employs the high-k transform algorithm implemented in NUMSBT
to calculate the overlap operators between
overlapping augmentation spheres that have different positions or elements.
It is also used to filter out high-frequency components from AE partial
waves, which allows projections from pseudowavefunctions to AE partial
waves to be performed in real space, which is a vital component of the
code.
NUMSBT is distributed under the Standard CPC License, and the algorithm is
developed in the following paper:
    * Talman, J. Computer Physics Communications 2009, 180, 332 –338.
4. **WaveTrans**:
reader.c and reader.h, which read WAVECAR files from VASP output,
are based on the Fortran program, WaveTrans, written by
R. M. Feenstra and M. Widom from the Dept. of Physics at Carnegie
Mellon University. To see the original work,
please visit <https://www.andrew.cmu.edu/user/feenstra/wavetrans/>.
5. **Doxygen**:
Doxygen is a documentation generator from which I built the docs for PAWpySeed.
It is an excellent tool that allows for clean, up-to-date documentaton
that is easy to make and navigate. Check it out at <http://www.stack.nl/~dimitri/doxygen/>.

## Citing

If you find pawpyseed useful, please encourage its development by
citing the following paper in your research:

```
@ARTICLE{2019arXiv190411572B,
       author = { {Bystrom}, Kyle and {Broberg}, Danny and {Dwaraknath}, Shyam and
         {Persson}, Kristin A. and {Asta}, Mark},
        title = "{Pawpyseed: Perturbation-extrapolation band shifting corrections for point defect calculations}",
      journal = {arXiv e-prints},
     keywords = {Condensed Matter - Materials Science},
         year = "2019",
        month = "Apr",
          eid = {arXiv:1904.11572},
        pages = {arXiv:1904.11572},
archivePrefix = {arXiv},
       eprint = {1904.11572},
 primaryClass = {cond-mat.mtrl-sci},
       adsurl = {https://ui.adsabs.harvard.edu/abs/2019arXiv190411572B},
      adsnote = {Provided by the SAO/NASA Astrophysics Data System}
}

```

You can view the paper here: <https://arxiv.org/abs/1904.11572>

## Questions and Comments

Find a bug? Areas of code unclearly documented? Other questions? Feel free to contact
Kyle Bystrom at kylebystrom@gmail.com with the subject "pawpyseed: (Your Topic)"
AND/OR create an issue on the Github page at <https://github.com/kylebystrom/pawpyseed>.
