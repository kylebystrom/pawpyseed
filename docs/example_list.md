# Examples and Tutorials

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

List of Examples:
1. [Overlap operators of pseudowavefunctions]({{ site.baseurl }}/examples/example1.md)
2. [AE overlap operators using PawpyData subclasses]({{ site.baseurl }}/examples/example2.md)