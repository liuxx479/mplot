from scipy import *
from pylab import *
import os

import numpy as np
# the purpose of this code is simply organize all the errors into 1 array

fisher_mat=np.ndarray(shape=(100,3))
i=0
for iseed in np.arange(0,100000,1000):
	f='/direct/astro+astronfs01/workarea/jia/fit_noise_seed%s/err_g6.976_ngal15_bins15.ls'%(iseed)
	fisher_mat[i]=np.genfromtxt(f)
	i+=1