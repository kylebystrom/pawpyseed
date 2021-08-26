# coding: utf-8

from pymatgen.io.vasp.inputs import Poscar, Potcar

from pawpyseed.wavefunction import *

if __name__ == "__main__":
    posb = Poscar.from_file("bulk/CONTCAR").structure
    posd = Poscar.from_file("Pcharge_0/CONTCAR").structure
    potb = Potcar.from_file("bulk/POTCAR")
    potd = Potcar.from_file("Pcharge_0/POTCAR")
    pwf1 = PseudoWavefunction("bulk/WAVECAR", "bulk/vasprun.xml")
    pwf2 = PseudoWavefunction("Pcharge_0/WAVECAR", "Pcharge_0/vasprun.xml")

    wf1 = Wavefunction(posb, pwf1, CoreRegion(potb), (120, 120, 120))
    wf2 = Wavefunction(posd, pwf2, CoreRegion(potd), (120, 120, 120))
    wf2.setup_projection(wf1)
    print(wf2.defect_band_analysis(wf1, spinpol=True))
    # for i in range(253,254):
    #       wf2.single_band_projection(i, wf1)

    wf1.free_all()
    wf2.free_all()
