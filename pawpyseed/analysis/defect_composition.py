import yaml
import matplotlib.pyplot as plt 
import numpy as np 
import os, subprocess
from pymatgen.io.vasp.outputs import Vasprun
from pymatgen.io.vasp.inputs import Poscar
from pymatgen import Spin

"""
f = open('ddres.yaml', 'r')
data = yaml.load(f)
f.close()

bs = []
vs = []
cs = []

for defect in data:
	bs = []
	vs = []
	cs = []
	dat = data[defect]
	for b in dat:
		bs.append(b)
		#bs.append(2*b)
		#bs.append(2*b+1)
		vs.append(dat[b][0][0])
		#vs.append(dat[b][0][1])
		cs.append(dat[b][1][0])
		#cs.append(dat[b][1][1])

	fig, ax1 = plt.subplots()
	ax1.set_xlabel('band')
	ax1.set_ylabel('valence', color='b')
	ax1.bar(bs, vs, color='b')
	ax1.set_ylim(0,1)
	ax2 = ax1.twinx()
	ax2.set_ylabel('conduction', color='r')
	ax2.bar(bs, cs, color='r')
	ax2.set_ylim(0,1)
	ax2.invert_yaxis()
	plt.title(defect)
	plt.savefig('../../../Si_stuff/new/'+defect+'.png')
"""

class PawpyData:

	def __init__(self, dos, structure, data):
		
		self.energies = energies
		if type(dos) == list:
			self.energies = dos[0]
			self.densities = dos[1]
			self.efermi = dos[2]
		else:
			self.energies = dos.energies
			self.densities = (dos.densities[Spin.up] + dos.densities[Spin.down]) / 2
			self.efermi = dos.efermi
		self.structure = structure
		self.data = data

	def as_dict(self):

		data = {}
		data['structure'] = Poscar(self.structure).get_string()
		data['energies'] = self.energies
		data['densities'] = self.densities
		data['efermi'] = self.efermi
		data['data'] = self.data
		return data

	def write_yaml(self, filename):

		data = self.as_dict()
		f = open(filename, 'w')
		yaml.dump(data, f)
		f.close()

	@classmethod
	def from_dict(cls, data):
		return cls([data['energies'], data['densities'], data['efermi']],
			data['structure'], data['data'])

	@classmethod
	def from_yaml(cls, filename):
		f = open(filename, 'r')
		data = yaml.load(f)
		return cls.from_dict(data)



class BulkCharacter(PawpyData):

	def __init__(self, dos, structure, band_dict):

		if type(dos) == list:
			self.energies = dos[0]
			self.densities = dos[1]
			self.efermi = dos[2]
		else:
			self.energies = dos.energies
			self.densities = (dos.densities[Spin.up] + dos.densities[Spin.down]) / 2
			self.efermi = dos.efermi
		self.structure = structure
		self.data = band_dict

	def plot(self, name):
		bs = []
		vs = []
		cs = []
		for b in self.data:
			bs.append(b)
			vs.append(self.data[b][0][0])
			vs.append(self.data[b][0][1])
			cs.append(self.data[b][1][0])
			cs.append(self.data[b][1][1])

		bs = np.array(bs) - np.mean(bs)
		cs = np.array(cs)
		vs = np.array(vs)
		fig, (ax1, ax3) = plt.subplots(2, 1, gridspec_kw = {'height_ratios':[3, 1]},
			figsize=[6.4,6.4])
		ax1.set_xlabel('band')
		ax1.set_ylabel('valence', color='b')
		ax1.bar(bs-0.2, vs[::2], width=0.4, color='b')
		ax1.bar(bs+0.2, vs[1::2], width=0.4, color='b')
		ax1.set_ylim(0,1)
		ax2 = ax1.twinx()
		ax2.set_ylabel('conduction', color='r')
		ax2.bar(bs-0.2, cs[::2], width=0.4, color='r')
		ax2.bar(bs+0.2, cs[1::2], width=0.4, color='r')
		ax2.set_ylim(0,1)
		ax2.invert_yaxis()
		plt.title(name + ' band character')
		#plt.savefig('BAND_'+name)

		ax3.plot(self.energies - self.efermi, self.densities)
		ax3.set_xlabel('Energy (eV)')
		ax3.set_ylabel('Totla DOS')
		ax3.set_xlim(-2,2)
		ax3.set_ylim(0,max(self.densities))
		#plt.savefig('DOS_'+name)
		plt.savefig(name)

	@staticmethod
	def makeit(generator):
		#Example: 
		#>>> def_lst = ['charge_1', 'charge_0', 'charge_-1']
		#>>> generator = Projector.setup_multiple_protections('bulk', def_lst)
		#>>> objs = BulkComposition.makeit()

		bcs = {}

		for wf_dir, wf in generator:
			vr = Vasprun(os.path.join(wf_dir, 'vasprun.xml'))
			dos = vr.tdos
			data = wf.defect_band_analysis(spinpol = True)
			bcs[wf_dir] = BulkCharacter(dos, wf.structure, data)

		return bcs



class BasisExpansion(PawpyData):

	def __init__(self, dos, structure, projs):

		self.energies = dos.energies
		if type(dos) == list:
			self.energies = dos[0]
			self.densities = dos[1]
			self.efermi = dos[2]
		else:
			self.energies = dos.energies
			self.densities = (dos.densities[Spin.up] + dos.densities[Spin.down]) / 2
			self.efermi = dos.efermi
		self.structure = structure
		self.data = projs


	@staticmethod
	def makeit(generator):
		#Example: 
		#>>> def_lst = ['charge_1', 'charge_0', 'charge_-1']
		#>>> generator = Projector.setup_multiple_protections('bulk', def_lst)
		#OR
		#>>> generator = Projector.setup_multiple_projections(*pycdt_dirs('.'))
		#
		#>>> objs = BasisExpansion.makeit()

		bes = {}

		for wf_dir, wf in generator:

			vr = Vasprun(os.path.join(wf_dir, 'vasprun.xml'))
			dos = vr.tdos
			basis = wf.basis
			expansion = np.zeros((wf.nband, basis.nband * basis.nwk * basis.nspin),
				dtype=np.complex128)
			for b in range(wf.nband):
				expansion[b,:] = wf.single_band_projection(b)
			bes[wf_dir] = BasisExpansion(dos, wf.structure, expansion)

		return bes

def pycdt_dirs(top_dir):

	bulk = os.path.join(top_dir, 'bulk')
	wfdirs = []
	for root, dirs, files in os.walk(top_dir):
		if 'bulk' in root:
			continue
		if 'OUTCAR' in files:
			wfdirs.append(root)

	return bulk, wfdirs

