# PAWpySeed

PAWpySeed is a parallelized Python and C tool for reading and
analyzing the optimized band structure and wave functions
of VASP DFT calculations. The code is written for the PAW
formalism developed by P.E. Blochl and implemented
in VASP.

## PAW

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

## Files

The projector functions and partial waves are unique
to each element and stored in the POTCAR file
used in a VASP calculation. The pseudowavefunction
is the part of the wavefunction optimized during a DFT
calculation and is stored in the WAVECAR output file
in VASP. PAWpySeed parses both files to retrieve
all parts of the full Kohn Sham wavefunctions.

## Implementation

* Python Interface
* Computationally intensive tasks in C
* Parallelized with openmp

## Functionality

* Read pseudowavefunctions
* Read projectors and partial waves
* Read pseudopotential, atomic charge
density, and other POTCAR data
* Perform projections and other operator
expectation values on full wavefunctions
* Localize orbitals with SCDM-k