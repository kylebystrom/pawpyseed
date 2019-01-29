# Overlap operators of pseudowavefunctions

Consider calculating the overlaps of the pseudowavefunctions of several point defect
structures with the bulk supercell. Specifically, a few different charge states
of the boron substitutional, phosphorous substitutional, and silicon vacancy in silicon.
We will use PAWpySeed to calculate these overlaps and use them to find the proportion
of the defect bands that project onto the conduction and valence bands of the bulk.

```
from pawpyseed.core.projector import Projector, Wavefunction
import yaml # for storing our data

dat = {} # dictionary of conduction/valence character values
# Directories of the VASP output for the defect structures
def_lst = ['Pcharge_0', 'Pcharge_1', 'charge_-1', 'charge_0', 'charge_1', 'charge_2', 'charge_-2', 'Bcharge_-1', 'Bcharge_0']
# Initialize the bulk from the bulk directory
basis = Wavefunction.from_directory('bulk')

# loop over the defect directories
for wf_dir in def_lst:
    wf = Wavefunction.from_directory(wf_dir)
    pr = Projector(wf, basis, pseudo=True)
    dat[wf_dir] = {}
    # loop over bands near the band gap
    for i in range(250, 262):
        # if pseudo is true, the overlap operators of the pseudowavefunctions
        # is evaluated, rather than of the all electron wavefunctions,
        # which is much faster but less quantitatively informative
        # v + c = 1 for pseudo = True
        # v is the valence band character and c is the conduction band character
        v, c = pr.proportion_conduction(i, spinpol=True)
        dat[wf_dir][i] = (v, c)
    print ('FINISHED DEFECT %s' % wf_dir)

f = open('res.yaml', 'w')
yaml.dump(dat, f)
f.close()
```

Now 'res.yaml' will contain the valence and conduction band characters for the selected bands
and structures in a yaml-formatted dictionary.
