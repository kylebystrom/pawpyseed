# coding: utf-8

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
		if self.energy_levels:
			data['energy_levels'] = self.energy_levels
		return data

	def write_yaml(self, filename):

		data = self.as_dict()
		f = open(filename, 'w')
		yaml.dump(data, f)
		f.close()

	@classmethod
	def from_dict(cls, data):
		if 'energy_levels' in data:
			return cls([data['energies'], data['densities'], data['efermi']],
				data['structure'], data['data'], data['energy_levels'])
		return cls([data['energies'], data['densities'], data['efermi']],
			data['structure'], data['data'])

	@classmethod
	def from_yaml(cls, filename):
		f = open(filename, 'r')
		data = yaml.load(f.read().encode('utf-8'))
		return cls.from_dict(data)



class BulkCharacter(PawpyData):

	def __init__(self, dos, structure, band_dict, energy_levels = None,
		vbm = None, cbm = None):

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
		self.energy_levels = energy_levels
		self.vbm = vbm
		self.cbm = cbm
		if vbm:
			self.efermi = self.vbm

	def plot(self, name, bandgap = None, spinpol = False, title=None):
		if title == None:
			title = name
		bs = []
		vs = []
		cs = []
		for b in self.data:
			bs.append(b)
			vs.append(self.data[b][0][0])
			vs.append(self.data[b][0][1])
			cs.append(self.data[b][1][0])
			cs.append(self.data[b][1][1])

		if bandgap == None and self.vbm != None and self.cbm != None:
			bandgap = self.cbm - self.vbm

		bs = np.array(bs) - np.mean(bs)
		cs = np.array(cs)
		vs = np.array(vs)
		if self.energy_levels == None:
			fig, (ax1, ax3) = plt.subplots(2, 1, gridspec_kw = {'height_ratios':[3, 1]},
				figsize=[6.4,6.4])
		else:
			fig, (ax1, ax3) = plt.subplots(2, 1, gridspec_kw = {'height_ratios':[1, 1]},
				figsize=[6.4,8])
		ax1.set_xlabel('band', fontsize=18)
		ax1.set_ylabel('valence', color='b', fontsize=18)
		ax1.bar(bs-0.2, vs[::2], width=0.4, color='b')
		ax1.bar(bs+0.2, vs[1::2], width=0.4, color='b')
		ax1.set_ylim(0,1)
		ax2 = ax1.twinx()
		ax2.set_ylabel('conduction', color='r', fontsize=18)
		ax2.bar(bs-0.2, cs[::2], width=0.4, color='r')
		ax2.bar(bs+0.2, cs[1::2], width=0.4, color='r')
		ax2.set_ylim(0,1)
		ax2.invert_yaxis()
		plt.title(title + ' band character', fontsize=20)
		#plt.savefig('BAND_'+name)

		if self.energy_levels == None:
			ax3.plot(self.energies - self.efermi, self.densities)
			ax3.set_xlabel('Energy (eV)', fontsize=18)
			ax3.set_ylabel('Total DOS', fontsize=18)
			ax3.set_xlim(-2,2)
			ax3.set_ylim(0,max(self.densities))
		elif type(list(self.energy_levels.values())[0]) == int:
			print('here')
			bs = list(self.energy_levels.keys())
			bmean = np.mean(bs)
			#print(self.energy_levels)
			for b in self.energy_levels:
				ax3.plot([b-0.4-bmean,b+0.4-bmean],
					[self.energy_levels[b] - self.efermi] * 2,
					color = 'r')
			if bandgap != None:
				bmin = min(bs) - bmean - 0.5
				bmax = max(bs) - bmean + 0.5
				ax3.plot([bmin,bmax], [0,0], color='black')
				ax3.plot([bmin,bmax], [bandgap,bandgap], color='black')
			#	bs.append(b)
			#	es.append(self.energy_levels[b])
			#bs = np.array(bs) - np.mean(bs)
			#es = np.array(es)
			#ax3.bar(bs, es - self.efermi)
			ax3.set_xlabel('band', fontsize=18)
			ax3.set_ylabel('Energy (eV)', fontsize=18)
		else:
			bs = list(self.energy_levels.keys())
			bmean = np.mean(bs)
			#print(self.energy_levels)
			cmap = plt.get_cmap('plasma')
			for b in self.energy_levels:
				i = 0
				delta = 0.8 / len(self.energy_levels[b])
				for en, occ in self.energy_levels[b]:
					color = cmap(1-occ)
					disp = i * delta - 0.4
					span = [b-bmean+disp, b-bmean+disp+delta]
					ax3.plot(span,
						[en - self.efermi] * 2,
						color = color)
					i += 1
			if self.vbm != None and self.cbm != None:
				bmin = min(bs) - bmean - 0.5
				bmax = max(bs) - bmean + 0.5
				ax3.plot([bmin, bmax], [0,0], color='black')
				ax3.plot([bmin, bmax], [self.cbm-self.vbm,self.cbm-self.vbm], color='black')
			elif bandgap != None:
				bmin = min(bs) - bmean - 0.5
				bmax = max(bs) - bmean + 0.5
				ax3.plot([bmin,bmax], [0,0], color='black')
				ax3.plot([bmin,bmax], [bandgap,bandgap], color='black')
			ax3.set_xlabel('band', fontsize=18)
			ax3.set_ylabel('Energy (eV)', fontsize=18)
		plt.savefig(name.replace(' ', '_'))

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
			data, energy_levels = wf.defect_band_analysis(num_above_ef=5, num_below_ef=5,
				spinpol = True, return_energies=True)
			bcs[wf_dir] = BulkCharacter(dos, wf.structure, data, energy_levels)

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
		if 'bulk' in root or 'dielectric' in root:
			continue
		if 'OUTCAR' in files:
			wfdirs.append(root)

	return bulk, wfdirs

