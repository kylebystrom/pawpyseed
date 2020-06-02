# coding: utf-8

import unittest
import os, subprocess, sys
import time
import scipy

import numpy as np
from numpy.testing import assert_almost_equal, assert_equal,\
						assert_raises, assert_array_almost_equal

from scipy.special import lpmn, sph_harm
from nose import SkipTest
from nose.tools import nottest
from nose.plugins.skip import Skip

from pawpyseed.core import pawpyc
from pawpyseed.core.tests import testc

class PawpyTestError(Exception):
	"""
	Class for handling errors that occur during execution
	of Python functions in pawpyseed
	"""
	def __init__(self, msg):
		self.msg = msg

MODULE_DIR = os.path.dirname(os.path.abspath(__file__))

from pawpyseed.core.wavefunction import Wavefunction
from pawpyseed.core.noncollinear import NCLWavefunction
from pawpyseed.core.momentum import MomentumMatrix

class TestMomentumMatrix:

	def setup(self):
		self.currdir = os.getcwd()
		os.chdir(os.path.join(MODULE_DIR, '../../../test_files'))
		self.initialize_wf_and_mm()

	def initialize_wf_and_mm(self):
		SIZE = 60
		self.wf = Wavefunction.from_directory('.')
		vol = self.wf.structure.volume
		self.realspace_wf = self.wf.get_state_realspace(0,0,0,
					dim=(SIZE,SIZE,SIZE), remove_phase=True)
		#self.realspace_chg = np.abs(self.realspace_wf)**2
		self.realspace_chg = self.wf.get_state_realspace_density(0,0,0,
									dim=(SIZE,SIZE,SIZE))
		self.recipspace_wf = np.fft.fftn(self.realspace_wf) / SIZE**3 * np.sqrt(vol)
		self.recipspace_chg = np.fft.fftn(self.realspace_chg) / SIZE**3 * vol
		self.mm_real = MomentumMatrix(self.wf, encut=3000)
		self.mm_direct = MomentumMatrix(self.wf)
		self.mm_direct2 = MomentumMatrix(self.wf, encut=self.wf.encut)

		self.ncl_wf = NCLWavefunction.from_directory('noncollinear')
		self.ncl_realspace_wf = self.ncl_wf.get_state_realspace(0,0,0,
					dim=(SIZE,SIZE,SIZE), remove_phase=True)
		self.ncl_recipspace_wf = (np.fft.fftn(self.ncl_realspace_wf[0]) / SIZE**3 * np.sqrt(vol),\
							np.fft.fftn(self.ncl_realspace_wf[1]) / SIZE**3 * np.sqrt(vol))

	def teardown(self):
		os.chdir(self.currdir)

	def test_ncl_transform(self):
		chg = np.sum(np.abs(self.ncl_recipspace_wf[0])**2) +\
				np.sum(np.abs(self.ncl_recipspace_wf[1])**2)
		assert_array_almost_equal(chg, 1, 3)

	def test_encut_insensitivity(self):
		res = self.mm_direct.get_momentum_matrix_elems(0,0,0,0,0,0)
		res2 = self.mm_direct2.get_momentum_matrix_elems(0,0,0,0,0,0)
		assert_almost_equal(res[0], 1, 4)
		assert_almost_equal(res2[0], 1, 4)
		assert_almost_equal(res[:6], res2[:6], 7)

	def test_get_momentum_matrix_elems(self):
		res = self.mm_direct.get_momentum_matrix_elems(0,0,0,0,0,0)
		grid = self.mm_direct.momentum_grid
		for i in range(grid.shape[0]):
			if (np.abs(grid[i]) < 3).all():
				#print(grid[i], res[i], self.recipspace_chg[grid[i][0],grid[i][1],grid[i][2]])
				assert_almost_equal(res[i],
					self.recipspace_chg[grid[i][0],grid[i][1],grid[i][2]],
					3)
		with assert_raises(ValueError):
			self.mm_direct.get_momentum_matrix_elems(0,0,0,0,-1,0)

	def test_get_reciprocal_fullfw(self):
		res = self.mm_real.get_reciprocal_fullfw(0,0,0)
		print("check size", np.sum(np.abs(res)**2))
		grid = self.mm_real.momentum_grid
		for i in range(grid.shape[0]):
			if (np.abs(grid[i]) < 2).all():
				#print(grid[i], res[i], self.recipspace_wf[grid[i][0],grid[i][1],grid[i][2]])
				assert_almost_equal(res[i],
					self.recipspace_wf[grid[i][0],grid[i][1],grid[i][2]],
					3)
		with assert_raises(ValueError):
			self.mm_real.get_reciprocal_fullfw(50,0,0)

	def test_g_from_wf(self):
		grid = self.mm_real.momentum_grid
		for i in range(grid.shape[0]):
			if (np.abs(grid[i]) < 2).all():
				#print(grid[i], self.mm_real.g_from_wf(0,0,0,0,0,0,grid[i]), self.recipspace_chg[grid[i][0],grid[i][1],grid[i][2]])
				assert_almost_equal(self.mm_real.g_from_wf(0,0,0,0,0,0,grid[i]),
					self.recipspace_chg[grid[i][0],grid[i][1],grid[i][2]],
					3)
		with assert_raises(ValueError):
			self.mm_real.g_from_wf(100,0,0,0,0,0,[0,0,0])
		
