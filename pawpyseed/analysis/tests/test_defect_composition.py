# coding: utf-8

import unittest
import os, subprocess, sys
import time
import scipy

import numpy as np
from numpy.testing import assert_almost_equal

from scipy.special import lpmn, sph_harm
from nose import SkipTest
from nose.tools import nottest

from pawpyseed.core.projector import Wavefunction, Projector

from pawpyseed.analysis.defect_composition import *

MODULE_DIR = os.path.dirname(os.path.abspath(__file__))

class PawpyTestError(Exception):
	"""
	Class for handling errors that occur during execution
	of Python functions in pawpyseed
	"""
	def __init__(self, msg):
		self.msg = msg

class TestBulkCharacter:

	def setup(self):
		os.chdir(os.path.join(MODULE_DIR, '../../../test_files'))
		generator = Projector.setup_multiple_projections('.', ['.', '.'])
		self.bcs = BulkCharacter.makeit(generator)

	def teardown(self):
		pass

	def test_plot(self):
		self.bcs['.'].plot('tst')
		self.bcs['.'].plot('tst', spinpol=True)

	def test_read_write(self):
		self.bcs['.'].write_yaml('test_file.yaml') 
		tst = BulkCharacter.from_yaml('test_file.yaml')
		assert type(tst.structure) == type(self.bcs['.'].structure)
		assert tst.structure == self.bcs['.'].structure
		assert (tst.energies == self.bcs['.'].energies).all()
		assert (tst.densities == self.bcs['.'].densities).all()
		assert tst.efermi == self.bcs['.'].efermi
		assert tst.data == self.bcs['.'].data
		assert tst.vbm == self.bcs['.'].vbm
		assert tst.cbm == self.bcs['.'].cbm
		assert tst.energy_levels == self.bcs['.'].energy_levels
		assert tst.nspin == 2
		assert tst.kws[0] == 0.5
		assert tst.kws[1] == 0.5

class TestBasisExpansion:

	def setup(self):
		os.chdir(os.path.join(MODULE_DIR, '../../../test_files'))
		generator = Projector.setup_multiple_projections('.', ['.', '.'])
		self.bes = BasisExpansion.makeit(generator)

	def teardown(self):
		pass

	def test_read_write(self):
		self.bes['.'].write_yaml('test_file.yaml') 
		tst = BasisExpansion.from_yaml('test_file.yaml')
		assert type(tst.structure) == type(self.bes['.'].structure)
		assert tst.structure == self.bes['.'].structure
		assert (tst.energies == self.bes['.'].energies).all()
		assert (tst.densities == self.bes['.'].densities).all()
		assert tst.efermi == self.bes['.'].efermi
		assert (tst.data == self.bes['.'].data).all()
		assert tst.vbm == self.bes['.'].vbm
		assert tst.cbm == self.bes['.'].cbm
