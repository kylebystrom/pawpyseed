# AE overlap operators for wavefunctions in different structures: The core of pawpyseed

In defect physics and other problems involving perturbed crystal structures, it is
potentially useful to describe the Kohn-Sham (KS) states of one structure in the
basis of the KS states of another structure. This is easy to do with pseudowavefunctions,
but since pseudowavefunctions are not orthonormal, the results are imprecise and qualitative.
This projection is difficult for the AE wavefunctions in the PAW formalism, as the augmentation
regions that define the wavefunction near the nuclei overlap with each other
and the plane-wave pseudowavefunction in nontrivial ways.

But pawpyseed can handle this! In fact, it's what pawpyseed is built for.
To do AE state projections, we need to initialize a `Projector` class from two `Wavefunction`
objects: a basis (or bulk, for defect problems) `Wavefunction`, and a perturbed (or defect)
`Wavefunction`.

```
from pawpyseed.core.projector import * # also imports the wavefunction module

wf = Wavefunction.from_directory('defect', setup_projectors=False)
basis = Wavefunction.from_directory('bulk', setup_projectors=False)
pr = Projector(wf, basis)
```

We choose `setup_projectors=False` for this code snippet because the `Projector`
class does it's own setup, which combines all the projector functions and partial
waves for the structures in `basis` and `wf`. **Keep in mind that the crystal
structures for `basis` and `wf` can be different but must have the
same lattice, k-points, and energy cutoff**.

Now that we have intialized a `Projector`, we can call some useful functions,
like `defect_band_analysis`, which projects the 40 bands in `wf` closest
to the Fermi level onto all of the bands in `basis`, and then sums
the magnitudes of these projections to give proportion valence
and conduction band character values for the bands in `wf`.

```
res = pr.defect_band_analysis(num_below_ef=20,
		num_above_ef=20, pseudo = False, spinpol = False)
```


How is this different than just using pseudowavefunctions? With AE wavefunctions,
the KS states of `basis` form an orthonormal\*
basis set for the space of periodic functions in the lattice. This gives
the results quantitative significance. In the case of defect problems,
if the valence and conduction proportions of a band sum to less
than 1, it might indicate that the state has some deep level character!
(It might also mean some high-energy states outside the basis set are
mixed in.)

When we're done, we simply free the data.

\* Theoretically, this basis set is complete when all states are included.
Of course, we always calculate a finite number of KS eigenstates,
so this basis will not be complete.
