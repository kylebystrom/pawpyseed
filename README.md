# PAWpySeed

<https://kylebystrom.github.io/pawpyseed/> <br/>

**Pawpyseed has converted from ctypes to Cython for its C interface!!! This should
make the code more portable and easier to use.**

**WARNING: PAWpySeed is still in early development. Documentation is
incomplete, and some features are not yet thoroughly tested.
The evaluation of overlap operators is tested, but a standard test
suite is not yet published, and some features still require more thorough testing.
Also, the setup.py script is only tested on a couple systems, and same
with the pip command.**

PAWpySeed is a parallelized Python and C tool for reading and
analyzing the optimized band structure and wave functions
of VASP DFT calculations. The code is written for the PAW
formalism developed by P.E. Blochl and implemented
in VASP.

## Installation

### The Easy Way

First, install `mkl-devel`:

```
pip install mkl-devel (--user)
```

If you are on a computing cluster with "module" options, or your
system has different MKL builds (particularly one is built with icc),
you need to be careful how you link your libraries. Find the directory
in which MKL was installed (this should correspond to the location of
your Python packages, e.g. if your Python packages are in
`~/.local/lib/python3.5/site-packages`, your MKL is in `~/.local/lib`
and `~/.local/include`).

Installation can also be performed
by cloning this repository and running the `setup.py` script
in the root directory of the repository.

```
python setup.py build
python setup.py install
```

OR you can install PAWpySeed with `pip`.

```
pip install pawpyseed
```

This has been tested on Scientific Linux 7 and Linux Mint 18,
but should work for systems that have the appropriate
packages and environment variables defined as described below.

### Advanced Options



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
"pseudowavefunctions" well described by plane waves. The
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

### Future Functionality

* Localize orbitals with SCDM-k
* Atomic Hartree Fock and GGA DFT
database for use in charge corrections
and other applications
* Read noncollinear pseudowavefunctions
* Convert PAW wavefunctions to NC wavefunctions
(for use in GW calculations)
* Perturbative charge corrections
* Read pseudopotential, atomic charge
density, and other POTCAR data
* Perform general operator
expectation values on full wavefunctions

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


## Questions and Comments

Find a bug? Areas of code unclearly documented? Other questions? Feel free to contact
Kyle Bystrom at kylebystrom@berkeley.edu with the subject "pawpyseed: <Topic>".
