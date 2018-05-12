import yaml
import matplotlib.pyplot as plt 
import numpy as np 

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

class PawpyData:

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

	@staticmethod
	def from_dict(data):

		dat.energies = data['energies']
		dat.densities = data['densities']
		dat.efermi = data['efermi']
		dat.structure = data['structure']
		dat.data = data['data']

		return dat

	@staticmethod
	def read_yaml(filename):
		f = open(filename, 'r')
		data = yaml.load(f)
		return self.read_dict(data)



class BulkCharacter(PawpyData):

	def __init__(self, band_dict, dos, structure):

		self.energies = dos.energies
		self.densities = (dos.densities[Spin.up] + dos.densities[Spin.down]) / 2
		self.efermi = dos.efermi
		self.structure = structure
		self.data = band_dict

	@staticmethod
	def makeit(generator):

		bcs = {}

		for wf_dir, basis, wf in generator:
			vr = Vasprun(os.path.join(wf_dir, 'vasprun.xml'))
			dos = vr.tdos
			data = wf.defect_band_analysis(basis, spinpol = True)
			bcs[wf_dir] = BulkCharacter(data, dos, vr.structure)

		return bcs



class BasisExpansion(PawpyData):

	def __init__(self, projs, dos, structure):

		self.energies = dos.energies
		self.densities = dos.densities
		self.efermi = dos.efermi
		self.structure = structure
		self.data = projs


	@staticmethod
	def makeit(generator):

		for wf_dir, basis, wf in generator:

			expansion = np.zeros((wf.nband, basis.nband), dtype=np.float64)
			for b in range(wf.nband):
				expansion[b,:] = single_band_projection(b, basis)


