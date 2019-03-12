# coding: utf-8

## @package pawpyseed.core.symmetry
# Utilities related to symmetry of the crystal structure,
# namely finding symmetrically identical k-points and the
# space group operators that map between them.

from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.core.operations import SymmOp
import numpy as np

def get_symmops(structure, symprec):
	"""
	Helper function to get the symmetry operations of the structure
	in the reciprocal lattice fractional coordinates.

	Args:
		structure (pymatgen.core.structure.Structure)
		symprec (number): symmetry precision for pymatgen SpacegroupAnalyzer
	"""

	sga = SpacegroupAnalyzer(structure, symprec * max(structure.lattice.abc))
	symmops = sga.get_symmetry_operations(cartesian = True)
	lattice = structure.lattice.matrix
	invlattice = structure.lattice.inv_matrix
	newops = []
	for op in symmops:
		newrot = np.dot(lattice, op.rotation_matrix)
		newrot = np.dot(newrot, invlattice)
		newtrans = np.dot(op.translation_vector, invlattice)
		newops.append(SymmOp.from_rotation_and_translation(
			newrot, newtrans))
	return newops

def get_nosym_kpoints(kpts, structure, init_kpts = None, symprec=1e-4,
	gen_trsym = True, fil_trsym = True):
	"""
	Starting with a set of k-points (kpts), finds all of the k-points
	that are symmetrically identical to a k-point in kpts by symmetry
	transformations of a crystal (structure).

	Args:
		kpts (np.ndarray shape=(n,3))
	"""

	allkpts = [] if init_kpts == None else [kpt for kpt in init_kpts]
	orig_kptnums = []
	op_nums = []
	symmops = get_symmops(structure, symprec)
	trs = []
	for i, op in enumerate(symmops):
		for k, kpt in enumerate(kpts):
			newkpt = np.dot(op.rotation_matrix, kpt)
			newkpt -= np.around(newkpt)
			newkpt[ abs(newkpt + 0.5) < 1e-5 ] = 0.5
			#if ((newkpt > 0.5+1e-6) + (newkpt < -0.5+1e-6)).any():
			#	continue
			if fil_trsym:
				if newkpt[2] < -1e-6 or \
					(abs(newkpt[2]) < 1e-6 and newkpt[1] < -1e-6) or \
					(abs(newkpt[2]) < 1e-6 and abs(newkpt[1]) < 1e-6 and newkpt[0] < -1e-6):
					continue
			unique = True
			for nkpt in allkpts:
				diff = (newkpt - nkpt) % 1
				oppdiff = 1 - diff
				tst = (np.abs(diff) < 1e-4) + (np.abs(oppdiff) < 1e-4)
				if ( tst.all() ):
					unique = False
					break
			if unique:
				allkpts.append(newkpt)
				orig_kptnums.append(k)
				op_nums.append(i)
				trs.append(0)
	if gen_trsym:
		for i, op in enumerate(symmops):
			for k, kpt in enumerate(kpts):
				newkpt = np.dot(op.rotation_matrix, kpt) * -1
				newkpt -= np.around(newkpt)
				newkpt[ abs(newkpt + 0.5) < 1e-5 ] = 0.5
				#if ((newkpt > 0.5+1e-6) + (newkpt < -0.5+1e-6)).any():
				#	continue
				if fil_trsym:
					if newkpt[2] < -1e-10 or \
						(abs(newkpt[2]) < 1e-6 and newkpt[1] < -1e-6) or \
						(abs(newkpt[2]) < 1e-6 and abs(newkpt[1]) < 1e-6 and newkpt[0] < -1e-6):
						continue
				unique = True
				for nkpt in allkpts:
					diff = (newkpt - nkpt) % 1
					oppdiff = 1 - diff
					tst = (np.abs(diff) < 1e-4) + (np.abs(oppdiff) < 1e-4)
					if ( tst.all() ):
						unique = False
						break
				if unique:
					allkpts.append(newkpt)
					orig_kptnums.append(k)
					op_nums.append(i)
					trs.append(1)
	return np.array(allkpts), orig_kptnums, op_nums, symmops, trs

def get_kpt_mapping(allkpts, kpts, structure, symprec=1e-4, gen_trsym = True):
	symmops = get_symmops(structure, symprec)
	orig_kptnums = []
	op_nums = []
	trs = []
	for nkpt in allkpts:
		match = False
		for i, op in enumerate(symmops):
			for k, kpt in enumerate(kpts):
				newkpt = np.dot(op.rotation_matrix, kpt)
				#if ((newkpt > 0.5+1e-6) + (newkpt < -0.5+1e-6)).any():
				#	continue
				diff = (newkpt - nkpt) % 1
				oppdiff = 1 - diff
				tst = (np.abs(diff) < 1e-4) + (np.abs(oppdiff) < 1e-4)
				if tst.all():
					match = True
					orig_kptnums.append(k)
					op_nums.append(i)
					trs.append(0)
					break
			if match:
				break
		if match:
			continue
		for i, op in enumerate(symmops):
			for k, kpt in enumerate(kpts):
				newkpt = np.dot(op.rotation_matrix, kpt) * -1
				#if ((newkpt > 0.5+1e-6) + (newkpt < -0.5+1e-6)).any():
				#	continue
				diff = (newkpt - nkpt) % 1
				oppdiff = 1 - diff
				tst = (np.abs(diff) < 1e-4) + (np.abs(oppdiff) < 1e-4)
				if tst.all():
					match = True
					orig_kptnums.append(k)
					op_nums.append(i)
					trs.append(1)
					break
			if match:
				break
		if not match:
			raise PAWpyError("Could not find kpoint mapping to %s" % str(nkpt))
	return orig_kptnums, op_nums, symmops, trs
