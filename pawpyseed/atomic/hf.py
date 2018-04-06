import numpy as np 
from ctypes import *

MODULE_DIR = os.path.dirname(os.path.abspath(__file__))
HFC = CDLL(os.path.join(MODULE_DIR, "hfc.so"))

HFC.get_E.restype = c_double
HFC.get_occs.restype = POINTER(c_double)
HFC.get_Ps.restype = POINTER(c_double)
HFC.setup_H.restype = POINTER(None)
HFC.setup.restype = POINTER(None)

if __name__ == '__main__':
	grid = np.exp(10**-5 * 1.02**np.arange(0,600,1,np.float64))
	hwf = HFC.setup_H(600, 92, 3, numpy_to_cdouble(grid))
	print(HFC.get_Ps(c_void_p(hwf)))
