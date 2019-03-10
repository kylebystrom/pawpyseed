# Desymmetrizing structures

Let's say you turned symmetry on for a set of calculations to save time,
and now you want to use pawpyseed to analyze the data. The problem is,
to do wavefunction projections, you need your two DFT calculations
to have the same set of k-points. Not a problem! As long as your
VASP calculations all had the same KPOINTS file, pawpyseed can
build Wavefunction objects with all the necessary k-points.

Let's say you have a 2x2x2 k-point mesh, but your `basis`
wavefunction has a lot of symmetry and therefore only
1 k-point. Due to these symmetry operations, it is possible
to extrapolate the wavefunctions at the other k-points.

(You can skip this paragraph if you just want to learn
how to tell pawpyseed to do desymmetrization for you.)
To "desymmetrize" a Wavefunction object manually, you can
simply call the `Wavefunction.desymmetrized_copy()`
function. This function optionally accepts a set of k-points
and weights, which you can use to specify an exact order
of the k-points and bypass the k-point-finding routine
that the code uses. This is important if you have an
already desymmetrized `Wavefunction` object that you want
to project onto/be projected onto, because the k-points
need to be in the same order with the same weights
for the projections to work.

To tell pawpyseed to do the desymmetrization when setting
up a `Projector` object, just pass `unsym_wf=True`
and/or `unsym_basis=True` to the initializer as
necessary. 