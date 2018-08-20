# AE overlap operators using PawpyData subclasses

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