# PAWpySeed

PAWpySeed is a parallelized Python and C tool for reading and
analyzing the optimized band structure and wave functions
of VASP DFT calculations. The code is written for the PAW
formalism developed by P.E. Blochl and implemented
in VASP.

## Installation

In the future, PAWpySeed will be installable with pip.
For the time being, installation should be performed
by cloning this repository and running the `setup.py` script
in the root directory of the repository.

```
python setup.py build
python setup.py install
```

The `build` command, in addition to the standard distutils setup,
compiles the C code in the `pawpyseed.core` module into a shared
object in the core module, `pawpy.so`. Currently, compiling the
C code requires the Intel C compiler, `icc`, as well as the Intel
Math Kernel Library. When the `build` command is run, an environment
variable `MKLROOT` must be present and point to the Math Kernel Library.

### Dependencies

All dependencies indicate the minimum version tested.
PAWpySeed might work fine with earlier versions, but
errors related to use of earlier versions are not supported
and will not be addressed.

Python requirements:
```
python>=3.5
numpy>=1.14
scipy>=1.0
pymatgen>=2018.2.13
```

C requirements:
```
icc >= 16.0.4 OR gcc >= 4.8.5
Intel Math Kernel Library >= 11.3.4
```
Intel MKL is available for free installation on a variety of platforms.
Most computing clusters will have Intel MKL, and you can install it
on your desktop (or any system to which you have root access) by following
the relevant instructions at the following URL:
<https://software.intel.com/en-us/articles/free-ipsxe-tools-and-libraries>.
After Intel MKL is installed, add the following line to your .bashrc
to link MKL (NOTE: this might not be the exact directory that MKL is in,
you need to check that first):
```
export MKLROOT=/opt/intel/compilers_and_libraries_2018/linux/mkl
export LD_LIBRARY_PATH=$MKLROOT/lib/intel64_lin:$LD_LIBRARY_PATH
export C_INCLUDE_PATH=$MKLROOT/include:$C_INCLUDE_PATH
```
The last line is optional but might be useful for future PAWpySeed builds
and other programs which make use of MKL.
The setup.py file will now take care of C compilation.

Optional Python dependencies:
```
sympy>=1.1.1
matplotlib>=0.2.5
```

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

## Examples

PAWpySeed is designed so that useful data can be generated
and formatted in a few lines of code. Example 1 illustrates
using the Wavefunction class utilities to calculate the overlap
operators of the pseudowavefunctions of several point defect structures
with the those of the bulk supercell of the same size (structures
generated with PyCDT https://bitbucket.org/mbkumar/pycdt). Example 2
shows how using the PawpyData subclasses can make some data generation
and formatting even easier.

NOTE: Having trouble writing a script to get the data you want?
Check out pawpyseed.analysis.defect_composition, which contains a few
examples of using PAWpySeed in the form of PawpyData subclasses.
The PawpyData.makeit function takes a generator from the Wavefunction.setup_multiple_projections function and uses it to calculate different values, depending on the subclass. These
subclasses cn serve as examples for how to calculate required values from a Wavefunction.

### Example 1: Overlap operators of pseudowavefunctions

Consider calculating the overlaps of the pseudowavefunctions of several point defect
structures with the bulk supercell. Specifically, a few different charge states
of the boron substitutional, phosphorous substitutional, and silicon vacancy in silicon.
We will use PAWpySeed to calculate these overlaps and use them to find the proportion
of the defect bands that project onto the conduction and valence bands of the bulk.

```
from pawpyseed.core.wavefunction import *
import yaml # for storing our data

dat = {} # dictionary of conduction/valence character values
# Directories of the VASP output for the defect structures
def_lst = ['Pcharge_0', 'Pcharge_1', 'charge_-1', 'charge_0', 'charge_1', 'charge_2', 'charge_-2', 'Bcharge_-1', 'Bcharge_0']
# Initialize the bulk from the bulk directory
basis = Wavefunction.from_directory('bulk')

# loop over the defect directories
for wf_dir in def_lst:
	wf = Wavefunction.from_directory(wf_dir)
	dat[wf_dir] = {}
	# loop over bands near the band gap
	for i in range(250, 262):
		# if pseudo is true, the overlap operators of the pseudowavefunctions
		# is evaluated, rather than of the all electron wavefunctions,
		# which is much faster but less quantitatively informative
		# v + c = 1 for pseudo = True
		# v is the valence band character and c is the conduction band character
		v, c = wf.proportion_conduction(i, basis, pseudo=True, spinpol=True)
		dat[wf_dir][i] = (v, c)
	print ('FINISHED DEFECT %s' % wf_dir)
	# Wavefunction objects use C code and memory, make sure to free it!
	wf.free_all()

basis.free_all()
f = open('res.yaml', 'w')
yaml.dump(dat, f)
f.close()
```

Now 'res.yaml' will contain the valence and conduction band characters for the selected bands
and structures in a yaml-formatted dictionary.

### Example 2: AE overlap operators using PawpyData subclasses

To calculate the AE overlap operators, it is first necessary to run the setup
routines which calculate the overlap operators of the projector functions
with the PS wavefunctions and with each other. This is demonstrated in the
following script:

```
from pawpyseed.core.wavefunction import *

basis = Wavefunction.from_directory('bulk')
wf = Wavefunction.from_directory('defect')

wf.setup_projection(basis, setup_basis=True)

v, c = wf.proportion_condition(253, basis, pseudo=False, spinpol=True)
```

The following example does the same thing as Example 1, except the all electron
overlap operators are calculated, and only four lines of code are used!

```
from pawpyseed.core.wavefunction import *
from pawpyseed.analysis.defect_composition import *

generator = Wavefunction.setup_multiple_projections(*pycdt_dirs('.'))
bcs = BulkCharacter.makeit(generator)
```

This script results in a list `bcs` of `BulkCharacter` objects, which each contain bulk valence
and conduction band character for 20 bands above and below the Fermi level for each defect, the
defect crystal structures as pymatgen Structure objects, and density of states data for each
defect, which can be useful for plotting and analysis later.

## Questions and Comments

Find a bug? Areas of code unclearly documented? Other questions? Feel free to contact
Kyle Bystrom at kylebystrom@berkeley.edu with the subject "pawpyseed: <Topic>".
