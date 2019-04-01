from pawpyseed.core.projector import Projector, Wavefunction
import numpy as np
from pymatgen.analysis.defects.corrections import DefectCorrection
from pawpyseed.analysis.defect_composition import BulkCharacter

class PerturbationCorrection(DefectCorrection):

	def __init__(self):
		self.metadata = {"defect_levels": [], "potalign": 0}

	def get_correction(self, d, filename):

		defect_path = d.parameters["path"]
		mb = True
		defect = BulkCharacter.from_yaml(filename).data

		potalign = d.parameters["potalign"]
		vbm = d.parameters["vbm"] + potalign
		cbm = d.parameters["cbm"] + potalign
	
		bulk_vbm = d.parameters["vbm"]
		bulk_cbm = d.parameters["cbm"]
		hybrid_vbm = d.parameters["hybrid_vbm"]
		hybrid_cbm = d.parameters["hybrid_cbm"]

		nband = d.parameters['nband']
		nwk = d.parameters['nwk']
		nspin = d.parameters['nspin']

		ens = []
		occs = []
		for b in range(nband):
			enset = []
			for item in d.parameters['eigenvalues'].values():
				for k in range(nwk):
					enset.append(item[k][b][0])
					occs.append(item[k][b][1])
			ens.append(enset)

		ens = np.array(ens)
		occs = np.array(occs)

		#GET projection amounts
		proj_amounts = {}
		corr = 0
	
		weights = d.parameters['kptweights']
		print ("cbmshift, vbmshift", hybrid_cbm - bulk_cbm, hybrid_vbm - bulk_vbm)
	
		num_vbm = 0	
		for band in defect.keys():
			proj_amounts[band] = {}
			for spin in range(nspin):
				band_occ = np.sum(occs[band*nwk*nspin+spin*nwk : band*nwk*nspin+(spin+1)*nwk] * weights)
				if nspin == 1:
					band_occ *= 2
				proj_amounts[band][spin] = (defect[band][0][spin], defect[band][1][spin])
				corr_term = (proj_amounts[band][spin][0] * (hybrid_vbm - bulk_vbm) \
					+ proj_amounts[band][spin][1] * (hybrid_cbm - bulk_cbm))
				en = np.sum(ens[band][spin*nwk:(spin+1)*nwk] * weights)
				new_en = en + corr_term
				if en < hybrid_vbm + potalign:
					num_vbm += band_occ
					print('HOLE IN VB', band, spin, band_occ*(hybrid_vbm-bulk_vbm))
					corr += band_occ * (hybrid_vbm - bulk_vbm)
				elif en > hybrid_cbm + potalign:
					print('ELEC IN CB', band, spin, band_occ*(hybrid_cbm-bulk_cbm))
					corr += band_occ * (hybrid_cbm - bulk_cbm)
				else:
					print('STATE IN GAP', band, spin, corr_term, band_occ)
					corr += band_occ * corr_term

		return corr, proj_amounts, num_vbm

class DelocalizedStatePerturbationCorrection(DefectCorrection):

	def __init__(self):
		self.metadata = {"defect_levels": [], "potalign": 0}

	def get_correction(self, d, filename):
		from pycdt.utils.plotter import SingleParticlePlotter
		from pycdt.utils.parse_calculations import SingleDefectParser

		defect_path = d.parameters["path"]

		local_bands = d.parameters["defect_ks_delocal_data"]["localized_band_indices"]
		contain_nums = d.parameters["defect_ks_delocal_data"]["contain_nums"]
		local_dat = {}
		for spin in local_bands:
			local_dat[spin] = dict(zip(local_bands[spin], contain_nums[spin]))

		mb = True
		defect = BulkCharacter.from_yaml(filename).data

		potalign = d.parameters["potalign"]
		vbm = d.parameters["vbm"] + potalign
		cbm = d.parameters["cbm"] + potalign
	
		bulk_vbm = d.parameters["vbm"]
		bulk_cbm = d.parameters["cbm"]
		hybrid_vbm = d.parameters["hybrid_vbm"]
		hybrid_cbm = d.parameters["hybrid_cbm"]

		nband = d.parameters['nband']
		nwk = d.parameters['nwk']
		nspin = d.parameters['nspin']

		ens = []
		occs = []
		for b in range(nband):
			enset = []
			for item in d.parameters['eigenvalues'].values():
				for k in range(nwk):
					enset.append(item[k][b][0])
					occs.append(item[k][b][1])
			ens.append(enset)

		ens = np.array(ens)
		occs = np.array(occs)

		#GET projection amounts
		proj_amounts = {}
		corr = 0
	
		weights = d.parameters['kptweights']
		print ("cbmshift, vbmshift", hybrid_cbm - bulk_cbm, hybrid_vbm - bulk_vbm)
	
		num_vbm = 0	
		for band in defect.keys():
			proj_amounts[band] = {}
			for spin in range(nspin):
				if band in local_dat[spin]:
					print('LOCALIZED BAND', band, spin, local_dat[spin][band])
					x = local_dat[spin][band]
					frac = min(1, max(1 - x, 0))
					print(frac)
				else:
					frac = 1
				band_occ = np.sum(occs[band*nwk*nspin+spin*nwk : band*nwk*nspin+(spin+1)*nwk] * weights)
				if nspin == 1:
					band_occ *= 2
				proj_amounts[band][spin] = (defect[band][0][spin], defect[band][1][spin])
				corr_term = (proj_amounts[band][spin][0] * (hybrid_vbm - bulk_vbm) \
					+ proj_amounts[band][spin][1] * (hybrid_cbm - bulk_cbm))
				corr_term *= frac
				en = np.sum(ens[band][spin*nwk:(spin+1)*nwk] * weights)
				new_en = en + corr_term
				if en < hybrid_vbm + potalign:
					num_vbm += band_occ
					print('HOLE IN VB', band, spin, band_occ*(hybrid_vbm-bulk_vbm))
					corr += band_occ * (hybrid_vbm - bulk_vbm)
				elif en > hybrid_cbm + potalign:
					print('ELEC IN CB', band, spin, band_occ*(hybrid_cbm-bulk_cbm))
					corr += band_occ * (hybrid_cbm - bulk_cbm)
				else:
					print('STATE IN GAP', band, spin, corr_term, band_occ)
					corr += band_occ * corr_term

		return corr, proj_amounts, num_vbm

