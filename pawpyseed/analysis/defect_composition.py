# coding: utf-8

## File presenting the PawpyData class
# and its children for analysis of pawpyseed
# output primarily from core.projector.Projector

import yaml
import matplotlib.pyplot as plt 
import numpy as np 
import os, subprocess
from pymatgen.io.vasp.outputs import Vasprun
from pymatgen.io.vasp.inputs import Poscar
from pymatgen import Spin

class PawpyData:

	def __init__(self, structure, data, dos = None, vbm = None, cbm = None):
		"""
		Arguments:
			structure (pymatgen.core.structure.Structure): crystal structure
			data: Whatever data is stored
			dos (pymatgen.electronic_structure.dos.DOS or list, None): A pymatgen
				density of states or a list containing: 1) energies, 2) density
				of states values at those energies, 3) the Fermi level.
			vbm (float, None): valence band maximum
			cbm (float, None): conduction band minimum
		"""
		if type(structure) == str:
			structure = Poscar.from_string(structure).structure
		if dos is None:
			self.energies = None
			self.densities = None
			self.efermi = None
		elif type(dos) == list:
			self.energies = dos[0]
			self.densities = dos[1]
			self.efermi = dos[2]
		else:
			self.energies = dos.energies
			self.densities = (dos.densities[Spin.up] + dos.densities[Spin.down]) / 2
			self.efermi = dos.efermi
		self.structure = structure
		self.data = data
		self.cbm = cbm
		self.vbm = vbm
		if (not (cbm is None)) and (not (vbm is None)):
			self.bandgap = max(0, cbm - vbm)
		else:
			self.bandgap = None

	def set_band_properties(vbm, cbm):
		"""
		Set the VBM, CBM, and band gap.

		Arguments:
			vbm (float): valence band maximum
			cbm (float): conduction band minimum
		"""
		self.bandgap = max(0, cbm - vbm)
		self.cbm = cbm
		self.vbm = vbm

	def as_dict(self):
		"""
		Returns a representation of the PawpyData
		as a dictionary.
		"""
		data = {}
		data['structure'] = Poscar(self.structure).get_string()
		data['energies'] = self.energies
		data['densities'] = self.densities
		data['data'] = self.data
		data['vbm'] = self.vbm
		data['cbm'] = self.cbm
		data['efermi'] = self.efermi
		return data

	def write_yaml(self, filename):
		"""
		Write the PawpyData as a yaml file
		called filename
		"""
		data = self.as_dict()
		f = open(filename, 'w')
		yaml.dump(data, f)
		f.close()

	@classmethod
	def from_dict(cls, data):
		"""
		Takes the dictionary--data--and
		returns a PawpyData instance.
		"""
		return cls(data['structure'], data['data'],
			dos = [data['energies'], data['densities'], data['efermi']],
			vbm = data['vbm'], cbm = data['cbm'])

	@classmethod
	def from_yaml(cls, filename):
		"""
		Reads a PawpyData instance from
		a file called filename.
		"""
		f = open(filename, 'r')
		data = yaml.load(f.read().encode('utf-8'), Loader=yaml.Loader)
		return cls.from_dict(data)



class BulkCharacter(PawpyData):
	"""
	The data member for the BulkCharacter is dictionary of form
	{band_index : (v, c)}, where v and c are the valence and
	conduction character of the band, respectively.
	"""

	def __init__(self, structure, data, energy_levels = None,
		dos = None, vbm = None, cbm = None, metadata=None):
		"""
		Arguments:
			structure (pymatgen.core.structure.Structure): crystal structure
			energy_levels (list of list): list of energy levels at the different
				k-points and spins in each band. Returned by Projector.defect_band_analysis
			data: {band_index : (v, c)}, where v and c are the valence and
				conduction character of the band, respectively.
				Returned by Projector.defect_band_analysis
			dos (pymatgen.electronic_structure.dos.DOS or list, None): A pymatgen
				density of states or a list containing: 1) energies, 2) density
				of states values at those energies, 3) the Fermi level.
			vbm (float, None): valence band maximum
			cbm (float, None): conduction band minimum
		"""

		self.energy_levels = energy_levels
		self.kws = None
		self.nspin = None
		if metadata is not None:
			if 'kws' in metadata:
				self.kws = metadata['kws']
			if 'nspin' in metadata:
				self.nspin = metadata['nspin']
		self.metadata = metadata
		super(BulkCharacter, self).__init__(structure, data, dos, vbm, cbm)

	def as_dict(self):
		data = super(BulkCharacter, self).as_dict()
		data['energy_levels'] = self.energy_levels
		data['metadata'] = {}
		if self.nspin is not None:
			data['metadata']['nspin'] = self.nspin
		if self.kws is not None:
			data['metadata']['kws'] = self.kws
		return data

	@classmethod
	def from_dict(cls, data):
		"""
		Takes the dictionary--data--and
		returns a PawpyData instance.
		"""
		if 'metadata' in data:
			metadata = data['metadata']
		else:
			metadata = None
		return cls(data['structure'], data['data'], energy_levels = data['energy_levels'],
				dos = [data['energies'], data['densities'], data['efermi']],
				vbm = data['vbm'], cbm = data['cbm'], metadata = metadata)

	def plot(self, name, title=None, spinpol = False):
		"""
		Plot the bulk character data. If the energy levels are available,
		those are plotted for reference. Otherwise, if the dos is available,
		it is plotted for reference.

		Arguments:
			name (str): Name of the file to which the plot is saved.
				Spaces will be replaced with underscores.
			title (str, None): Title of the plot. Defaults to name,
				so it is recommended to set this.
			spinpol (bool, False): Show the spin polarized valence
				and conduction character. Only works if the VASP
				calculation was spin polarized
		"""
		kws = self.kws if (self.kws is not None) else np.array([])
		if spinpol:
			if self.nspin == 2:
				spin = 2
			else:
				raise ValueError('Must verify that nspin==2 to use spinpol')
		else:
			spin = 1
			if self.nspin == 2:
				kws = np.append(kws, kws)
			elif self.nspin == None and len(kws) == len(list(self.energy_levels.values())[0])//2:
				kws = np.append(kws, kws)

		if (kws is not None) and (len(kws) != len(list(self.energy_levels.values())[0]))//spin:
			raise ValueError('kws is the wrong shape')

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

		bandgap = self.bandgap

		bs = np.array(bs) - np.mean(bs)
		xticks = bs[:]
		cs = np.array(cs)
		vs = np.array(vs)
		if self.energy_levels == None:
			if self.energies is None or self.efermi is None or self.densities is None:
				fig, ax1 = plt.subplots()
			else:
				fig, (ax1, ax3) = plt.subplots(2, 1, gridspec_kw = {'height_ratios': [3, 1]},
					figsize=[6.4,6.4])
		else:
			fig, (ax1, ax3) = plt.subplots(2, 1, gridspec_kw = {'height_ratios': [1, 1]},
				figsize=[6.4,8])
		ax1.set_xlabel('band', fontsize=18)
		ax1.set_xticks(xticks)
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
			if not (self.energies is None or self.efermi is None or self.densities is None):
				ax3.plot(self.energies - self.efermi, self.densities)
				ax3.set_xlabel('Energy (eV)', fontsize=18)
				ax3.set_ylabel('Total DOS', fontsize=18)
				ax3.set_xlim(-2,2)
				ax3.set_ylim(0,max(self.densities))
		else:
			bs = list(self.energy_levels.keys())
			bmean = np.mean(bs)
			cmap = plt.get_cmap('rainbow')
			for b in self.energy_levels:
				i = 0
				delta = 0.8 / spin# len(self.energy_levels[b])
				enlists = []
				occlists = []
				for s in range(spin):
					length = len(self.energy_levels[b])
					interval = (s*length//spin, (s+1)*length//spin)
					enlists.append([self.energy_levels[b][t][0] for t in range(*interval)])
					occlists.append([self.energy_levels[b][t][1] for t in range(*interval)])
				for i, endat in enumerate(zip(enlists, occlists)):
					enlist, occlist = endat
					color = cmap(1-np.mean(occlist))
					disp = i * delta - 0.4
					span = [b-bmean+disp, b-bmean+disp+delta]
					if self.kws is None:
						en = np.mean(enlist)
					else:
						en = np.dot(enlist, kws)
					errs = (en-min(enlist), max(enlist)-en)
					ax3.bar(b-bmean+disp+delta/2, max(enlist)-min(enlist),
						width=delta, bottom=min(enlist)-self.efermi, color='0.8')
					ax3.plot(span, [en - self.efermi] * 2,
						color = color)
				#for en, occ in self.energy_levels[b]:
				#	color = cmap(1-occ)
				#	disp = i * delta - 0.4
				#	span = [b-bmean+disp, b-bmean+disp+delta]
				#	ax3.plot(span,
				#		[en - self.efermi] * 2,
				#		color = color)
				#	i += 1
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
			ax3.set_xticks(xticks)
			ax3.set_ylabel('Energy (eV)', fontsize=18)
		plt.tight_layout()
		plt.savefig(name.replace(' ', '_'))

	@staticmethod
	def makeit(generator):
		#Example: 
		#>>> def_lst = ['charge_1', 'charge_0', 'charge_-1']
		#>>> generator = Projector.setup_multiple_protections('bulk', def_lst)
		#>>> objs = BulkComposition.makeit(generator)

		bcs = {}

		for wf_dir, pr in generator:
			vr = Vasprun(os.path.join(wf_dir, 'vasprun.xml'))
			bg, cbm, vbm, _ = vr.eigenvalue_band_properties
			dos = vr.tdos
			data, energy_levels = pr.defect_band_analysis(num_above_ef=5, num_below_ef=5,
				spinpol = True, return_energies=True)
			bcs[wf_dir] = BulkCharacter(pr.wf.structure, data,
				energy_levels = energy_levels, dos = dos, vbm = vbm, cbm = cbm,
				metadata = {'nspin': pr.wf.nspin, 'kws': pr.wf.kws})

		return bcs


class BasisExpansion(PawpyData):
	"""
	The data member for the BasisExpansion is 2D array of shape
	(wf.nband, basis.nband * basis.nwk * basis.nspin).
	Each item is the projection of a band of wf onto a band of
	basis for a given k-point index and spin index.
	"""

	@staticmethod
	def makeit(generator):
		#Example: 
		#>>> def_lst = ['charge_1', 'charge_0', 'charge_-1']
		#>>> generator = Projector.setup_multiple_protections('bulk', def_lst)
		#OR
		#>>> generator = Projector.setup_multiple_projections(*pycdt_dirs('.'))
		#
		#>>> objs = BasisExpansion.makeit(generator)

		bes = {}

		for wf_dir, pr in generator:

			vr = Vasprun(os.path.join(wf_dir, 'vasprun.xml'))
			bg, cbm, vbm, _ = vr.eigenvalue_band_properties
			dos = vr.tdos
			basis = pr.basis
			expansion = np.zeros((pr.wf.nband, basis.nband * basis.nwk * basis.nspin),
				dtype=np.complex128)
			for b in range(pr.wf.nband):
				expansion[b,:] = pr.single_band_projection(b)
			bes[wf_dir] = BasisExpansion(pr.wf.structure, expansion, dos=dos,
				vbm = vbm, cbm = cbm)

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

