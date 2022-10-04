from slfit import slfit

import numpy as np

ratiofile = "ratio_4_SLfit_test_Q49.txt"

ratio = np.loadtxt(ratiofile)

slfit (ratio, 50e-3, 5e-3, 600.13, 4.5, "1D3Z_f.pdb")