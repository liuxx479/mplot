import numpy as np
from scipy import *

s_arr = (0.0, 0.5, 1.0, 1.5, 2.0)
bins_arr =  (10, 20, 30, 40, 50, 60, 70, 80, 90, 100)

fit_dir = '/direct/astro+u/jia/magbias/peaks/'

#### put each chisqcontour fitting file, powell fitting into 1 single file ####
filename = fit_dir+'chiscontour_true.ls'

chisqells = np.ndarray(shape=(0,7))
for bins in bins_arr:
	f=abs(np.genfromtxt(fit_dir+'chisqcontour_true_bins%s.ls'%(bins))[1,1:])
	f=np.append(f,bins)
	chisqells=np.append(chisqells,[f],axis=0)

chisqells.shape
savetxt(filename,chisqells)



#### brute force fitting file ##########

file_cube = (fit_dir+'true_ells_brute.ls', fit_dir+'bias_ells_brute.ls')

# bias
mat_ells = np.ndarray(shape=(0,21))
for s in s_arr:
	for bins in bins_arr:
		file_ellss = fit_dir+'bias_ells_brute_s%s_bins%s.ls'%(s,bins)
		ells = np.genfromtxt(file_ellss)
		mat_ells = np.append(mat_ells,[ells],axis=0)

np.savetxt(file_cube[1],mat_ells)
  
# true
mat_ells = np.ndarray(shape=(0,21))
for bins in bins_arr:
	file_ellss = fit_dir+'true_ells_brute_bins%s.ls'%(bins)
	ells = np.genfromtxt(file_ellss)
	mat_ells = np.append(mat_ells,[ells],axis=0)

np.savetxt(file_cube[0],mat_ells)

