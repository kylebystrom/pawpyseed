# AE overlap operators using PawpyData subclasses

To calculate the AE overlap operators for KS states from two different
structures, there is some tricky and computationally intensive math
that goes into setting up the problem. Intializing a `Wavefunction` object
does not automatically do this work. Instead, you need to initialize
a `Projector` class. This is demonstrated in the
following script:

```
from pawpyseed.core.projector import *

basis = Wavefunction.from_directory('bulk')
wf = Wavefunction.from_directory('defect')

pr = Projector(wf, basis)

v, c = wf.proportion_condition(253, basis, pseudo=False, spinpol=True)
```

The following example does the same thing as Example 1, except the all electron
overlap operators are calculated, and only four lines of code are used!
(pycdt_dirs is a function which searches for a directory named bulk
and returns it as well as a list of all other directories containing
complete VASP output.)

```
from pawpyseed.core.wavefunction import *
from pawpyseed.analysis.defect_composition import *

generator = Projector.setup_multiple_projections(*pycdt_dirs('.'))
bcs = BulkCharacter.makeit(generator)
```

This script results in a list `bcs` of `BulkCharacter` objects, which each contain bulk valence
and conduction band character for 20 bands above and below the Fermi level for each defect, the
defect crystal structures as pymatgen Structure objects, and density of states data for each
defect, which can be useful for plotting and analysis later.