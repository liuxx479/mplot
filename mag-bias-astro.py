import numpy as np
from scipy import *
from scipy import interpolate
from pylab import *
import pyfits
import scipy.ndimage as snd
import os #math
from findpeak import findpeak
from scipy import optimize
import sys


########## constants & input ######

#s, bins = np.array(sys.argv[1:3],dtype=float)
s=float(sys.argv[1])
bins=int(sys.argv[2])

s_arr = (0.0, 0.5, 1.0, 1.5, 2.0)
bins_arr =  (10, 20, 30, 40, 50, 60, 70, 80, 90, 100)
initguess = (0.27, -1.1, 0.76) # initial guesses for fitter
sigma = 6.976
#sigma = 9.865

low, high = -0.02, 0.19 # -0.15 , 0.25 kappa range
kmapnum=1000 # number of maps per cosmology
initguess = (0.27, -1.1, 0.76) # initial guesses for fitter

kmap_dir = '/direct/astro+u/jia/kmap/' # where all the map directories locate
peak_dir = '/direct/astro+u/jia/magbias/peaks/' # folder to save peak files
fit_dir = '/direct/astro+u/jia/magbias/fit/' # folder to save fit files
plot_dir = '/direct/astro+u/jia/magbias/plot/' # folder to save plots

#kmap_dir = '/Users/jia/Documents/weaklensing/map2/' # where all the map directories locate
#peak_dir = '/Users/jia/Documents/weaklensing/magbias/peaks/' # folder to save peak files
#plot_dir = '/Users/jia/Documents/weaklensing/magbias/plot/' # folder to save plots
#fit_dir = '/Users/jia/Documents/weaklensing/magbias/fit/' # folder to save fit files

cosmo_all = os.listdir(kmap_dir) #a list of 7 cosmologies, suppose all the map folders are in kmap_dir

########### equations & functions that takes (bins, s, sigma) ########
#bin_size = (high - low) / bins
#hist_bins = np.linspace (low+bin_size/2, high-bin_size/2, bins) # center of each bin

smooth = lambda kmap, sigma: snd.filters.gaussian_filter(kmap,sigma)
bias = lambda kmap, sigma, s: snd.filters.gaussian_filter(kmap*(1+(5*s-2)*kmap),sigma)/snd.filters.gaussian_filter(1+(5*s-2)*kmap,sigma)

def peaklist(fitsfile, s=s, bins=bins, sigma=sigma, control=0):
	'''Grab kmap from fits file and find peaks for bias & true
	Input: 	fits file, 
		sigma for smoothing the map
		s for bias
	Output:	array of peaks in smoothed true map,
		array of peaks in smoothed bias map
	Procedure: First look into the peakfile direcotry to see if the peak count file exist, if yes then grab from there. If no, then apply bias to true map, count peaks for both maps, and save into a ascii file in the peak_dir
	'''
	if control==0:
		peakfile_true = peak_dir+fitsfile[-93:-4]+'_true_g%s.ls'%(sigma) ## only dep on sigma not s
		if os.path.isfile(peakfile_true):
			peak_true = np.genfromtxt(peakfile_true)
		else:
			hdulist = pyfits.open(fitsfile)#, ignore_missing_end = True)
			kmap = np.array(hdulist[0].data)
			kmap_true = smooth (kmap, sigma)
			peak_true = findpeak(kmap_true)
			np.savetxt(peakfile_true, peak_true)
			print 'counting true peaks for',fitsfile[-93:-4]
		return peak_true
	if control==1:
		peakfile_true = peak_dir+fitsfile[-93:-4]+'_true_g%s.ls'%(sigma) ## only dep on sigma not s
		peakfile_bias = peak_dir+fitsfile[-93:-4]+'_bias_s%s_g%s.ls'%(s,sigma)
		filestate1,filestate2 = 0,0
		if os.path.isfile(peakfile_true):
			peak_true = np.genfromtxt(peakfile_true)
			filestate1 = True
		if os.path.isfile(peakfile_bias):
			peak_bias = np.genfromtxt(peakfile_bias)
			filestate2 = True
		if filestate2 == False:
			hdulist = pyfits.open(fitsfile)#, ignore_missing_end = True)
			kmap = np.array(hdulist[0].data)
			if filestate1 == False: # no true peaks
				kmap_true = smooth (kmap, sigma)
				peak_true = findpeak(kmap_true)
				np.savetxt(peakfile_true, peak_true)
				print 'counting true peaks for',fitsfile[-93:-4]
			kmap_bias = bias (kmap, sigma, s)
			peak_bias = findpeak(kmap_bias)
			print 'counting bias peaks for',fitsfile[-93:-4]
			np.savetxt(peakfile_bias, peak_bias)
		return peak_true, peak_bias

def batch_kmap(kmap_dir,cosmo_num, s=s, bins=bins, sigma=sigma, control=0):
	'''Batch operate on 1000 maps (with specific s & theta_g(sigma)), or how many fits files that are in the cosmology directory
	Input:	kmap directory
		specific cosmology (int 0-6)
	Output:	2 matrices (true, bias)
		shape = (number of maps * bins)
	Procedure: First look into the peakfile direcotry to see if the peak matrix file exist, if yes then grab from there. If no, then create an empty matrix, prevoke peaklist routine to generate individual map peak count. Make a histogram out of the peakfile, and append each map's peak count to peak_matrix_true & bias.
	'''
	cosmo = cosmo_all[cosmo_num]
	kmaps = os.listdir(kmap_dir+cosmo)
	if control == 0:
		file_matrix_true = fit_dir+'peak_matrix_true_'+str(kmapnum)+cosmo+'_g%s_bins%s.ls'%(sigma,bins)
		# detect see if matrix file already exist, if yes, then grab from there, if not then calculate from peak_lists
		if os.path.isfile(file_matrix_true):
			peak_matrix_true = np.genfromtxt(file_matrix_true)
		else:
			peak_matrix_true = np.ndarray (shape = (len(kmaps), bins), dtype = float32)
			for n in range(len(kmaps)):
				fitsfile = kmap_dir+cosmo+'/'+kmaps[n]
				peak_true, peak_bias = peaklist(fitsfile, s=s, bins=bins, sigma=sigma, control=0)#!! to change!!
				peak_matrix_true[n] = np.histogram(peak_true, range=(low,high), bins=bins)[0]
				print 'Now making histogram for fitsfile #'+str(n)
				#		print 'n, len_true-bias, sum_diff_matrix[n]', n, len(peak_true)-len(peak_bias), sum(peak_matrix_true[n]-peak_matrix_bias[n])
			np.savetxt(file_matrix_true,peak_matrix_true)
		return peak_matrix_true
	if control == 1:
		file_matrix_true = fit_dir+'peak_matrix_true_'+str(kmapnum)+cosmo+'_g%s_bins%s.ls'%(sigma,bins)
		file_matrix_bias = fit_dir+'peak_matrix_bias_'+str(kmapnum)+cosmo+'_s%s_g%s_bins%s.ls'%(s,sigma,bins)
		# detect see if matrix file already exist, if yes, then grab from there, if not then calculate from peak_lists
		if os.path.isfile(file_matrix_true) and os.path.isfile(file_matrix_bias):
			peak_matrix_true = np.genfromtxt(file_matrix_true)
			peak_matrix_bias = np.genfromtxt(file_matrix_bias)
		else:
			peak_matrix_true = np.ndarray (shape = (len(kmaps), bins), dtype = float32)
			peak_matrix_bias = np.ndarray (shape = (len(kmaps), bins), dtype = float32)
			for n in range(len(kmaps)):
				fitsfile = kmap_dir+cosmo+'/'+kmaps[n]
				peak_true, peak_bias = peaklist(fitsfile, s=s, bins=bins, sigma=sigma,control=1)
				peak_matrix_true[n] = np.histogram(peak_true, range=(low,high), bins=bins)[0]
				peak_matrix_bias[n] = np.histogram(peak_bias, range=(low,high), bins=bins)[0]
				print 'Now making histogram for fitsfile #'+str(n)
				#		print 'n, len_true-bias, sum_diff_matrix[n]', n, len(peak_true)-len(peak_bias), sum(peak_matrix_true[n]-peak_matrix_bias[n])
			np.savetxt(file_matrix_true,peak_matrix_true)
			np.savetxt(file_matrix_bias,peak_matrix_bias)
		return peak_matrix_true, peak_matrix_bias

################# set parameters ########################
params = np.ndarray(shape=(len(cosmo_all),3))
for i in range(params.shape[0]):
	iOm = cosmo_all[i][12:17]
	iw = cosmo_all[i][27:33]
	isi = cosmo_all[i][44:49]
	params[i] = iOm, iw, isi

Om = np.array([0.23, 0.26, 0.29]) # Omega Mass
w  = np.array([-1.2, -1.0, -0.8]) # W
si = np.array([0.750, 0.798, 0.85]) # Sigma8

hi_Om = int(np.where(params[:,0]==np.amax(params[:,0]))[0])
lo_Om = int(np.where(params[:,0]==np.amin(params[:,0]))[0])
hi_w  = int(np.where(params[:,1]==np.amax(params[:,1]))[0])
lo_w  = int(np.where(params[:,1]==np.amin(params[:,1]))[0])
hi_si = int(np.where(params[:,2]==np.amax(params[:,2]))[0])
lo_si = int(np.where(params[:,2]==np.amin(params[:,2]))[0])
contr = int(np.delete(np.arange(len(cosmo_all)),[hi_Om,hi_si,hi_w,lo_Om,lo_si,lo_w]))

cosmo_str = ['']*7
cosmo_str[hi_Om]='hi_Om'
cosmo_str[hi_si]='hi_si'
cosmo_str[hi_w] ='hi_w'
cosmo_str[lo_Om]='lo_Om'
cosmo_str[lo_si]='lo_si'
cosmo_str[lo_w] ='lo_w'
cosmo_str[contr]='contr'

##### create true_mat & bias_mat for certain (bins, s, sigma) ####
##### in preparation for chi-sq fitting ####

true_mat = np.ndarray(shape=(len(cosmo_all), kmapnum, bins))
#bias_mat = np.ndarray(shape=(len(cosmo_all), kmapnum, bins))

for i in np.arange(len(cosmo_all)):
#	print 'cosmology ' + str(i) + ' of ' + str(len(cosmo_all))
	if i == contr:
		true_mat[i], bias_mat = batch_kmap(kmap_dir,i,control=1) # here we create the actual files, most computing time
	else:
		true_mat[i] = batch_kmap(kmap_dir,i,control=0)

fidu_mat = average(true_mat,axis=1)
#fidu_bias_mat = average(bias_mat,axis=1)


########### extrapolate Ni_arr for certain bins ##############

interpOm = lambda ibin: interpolate.UnivariateSpline(Om,fidu_mat[[lo_Om,contr,hi_Om], ibin], k=2, s=0)
interpw  = lambda ibin: interpolate.UnivariateSpline(w, fidu_mat[[lo_w ,contr,hi_w ], ibin], k=2, s=0)
interpsi = lambda ibin: interpolate.UnivariateSpline(si,fidu_mat[[lo_si,contr,hi_si], ibin], k=2, s=0)

Ni = lambda ibin, iOm, iw, isi: interpOm(ibin)(iOm)+interpw(ibin)(iw)+ interpsi(ibin)(isi)-2.0*fidu_mat[contr,ibin]

Ni_arr = lambda iOm, iw, isi: np.array(map(Ni, np.arange(bins),array([iOm]*bins),array([iw]*bins),array([isi]*bins)))


################ chi-square fitting ##############################
cov_mat = np.cov(true_mat[contr],rowvar=0)
cov_inv = np.mat(cov_mat).I

def chisq (initguess, ibmap):
	iOm, iw, isi = initguess
	del_Ni = np.mat(Ni_arr(iOm, iw, isi) - ibmap)
	chisquare = float64(del_Ni * cov_inv * del_Ni.T)
#	print chisquare, initguess
	return chisquare

def minimizechisq (ibmap, initguess=initguess):
	a = optimize.minimize(chisq,initguess,args=(ibmap,),method='L-BFGS-B')
#	print'%.2f\t%.4f\t%.4f\t%.3f'%(a.fun,a.x[0],a.x[1],a.x[2])
	return a.fun, a.x[0],a.x[1],a.x[2]

######## find fit to each map, prepare for chi-sq contour #######
def fitdistrib (s=s, bins=bins, sigma=sigma):
	'''Input: s, bins
	Output: 2 99*4 arrays (true fit, bias fit)
	each row is chisquare of the fit, and fitted 3 parameters
	'''	
	
	true_fit = np.ndarray(shape=(true_mat.shape[1],4))
	bias_fit = np.ndarray(shape=(true_mat.shape[1],4))
	
	true_fit_filename = fit_dir+'true_fit_g%s_bins%s.ls'%(sigma,bins)
	bias_fit_filename = fit_dir+'bias_fit_s%s_g%s_bins%s.ls'%(s,sigma,bins)
	
	if os.path.isfile(true_fit_filename) and os.path.isfile(bias_fit_filename):
		true_fit = np.genfromtxt(true_fit_filename)
		bias_fit = np.genfromtxt(bias_fit_filename)
	else:
		if os.path.isfile(true_fit_filename):
			true_fit = np.genfromtxt(true_fit_filename)
			n = 0
			while n < true_fit.shape[0]:
				ibmap_bias = bias_mat[n] #changed from bias_mat[contr,n]
				b = minimizechisq (ibmap_bias)
				bias_fit[n] = b
				print 'fitting to map #',n,'bias',b
				n+=1		
			savetxt(bias_fit_filename,bias_fit)
		else:
			n = 0
			while n < true_fit.shape[0]:
				ibmap_true = true_mat[contr,n]
				ibmap_bias = bias_mat[n] #changed from bias_mat[contr,n]
				a = minimizechisq (ibmap_true)
				b = minimizechisq (ibmap_bias)
				true_fit[n] = a
				bias_fit[n] = b
				print 'fitting to map #',n,'true',a,'bias',b
				n+=1		
			savetxt(true_fit_filename,true_fit)
			savetxt(bias_fit_filename,bias_fit)
	return true_fit, bias_fit

true_fit, bias_fit = fitdistrib (s=s, bins=bins, sigma=sigma)


############Method 2 - chi-sq contour use grid search first, then fit an ellipsoid to it ##########

def chisqgrid(chisq0=3.53,Ns=60):
	lower = np.sort(true_fit,axis=0)[3,1:]
	upper = np.sort(true_fit,axis=0)[-3,1:]
	filename = fit_dir+'chisqgrid_true_g%s_bins%s.ls'%(sigma,bins)
	if os.path.isfile(filename):
		chisqgrid = genfromtxt(filename).reshape(Ns,Ns,Ns)
	else:
		print 'chisqgrid running'
		chisqgrid = np.ndarray(shape=(Ns,Ns,Ns))
		i = 0
		for iOm in linspace(lower[0],upper[0],Ns):
			j = 0
			for iw in linspace(lower[1],upper[1],Ns): 
				k = 0
				for isi in linspace(lower[2],upper[2],Ns):
					ichisq = chisq((iOm,iw,isi),fidu_mat[contr])
					chisqgrid[i,j,k] = ichisq
					k += 1
				j += 1
			i += 1
		chisqgrid_reshape = chisqgrid.reshape(-1)
		savetxt(filename,chisqgrid_reshape)
	return chisqgrid,lower,upper

print 'Now searching in chisq grid, Ns = 60'
chisqgrid,lower,upper = chisqgrid(chisq0=3.53,Ns=60)

#### added 2/13/2012 Fisher matrix ####

hi_params = (hi_Om, hi_w, hi_si)
lo_params = (lo_Om, lo_w, lo_si)
paramvalue = [Om, w, si]

savetxt(fit_dir+'covmat_g%s_bins%s.ls'%(sigma,bins),cov_mat)


errfile=fit_dir+'err_g%s_bins%s.ls'%(sigma,bins)
if os.path.isfile(errfile):
	err = genfromtxt(errfile)
else:
	hi_arr = fidu_mat[[hi_si,hi_w,hi_Om]]
	lo_arr = fidu_mat[[lo_si,lo_w,lo_Om]]
	contr_arr = fidu_mat[contr]

	devhi = np.array([0.052,0.2,0.03])
	devlo = np.array([0.048,0.2,0.03])

	Fisher = np.ndarray(shape=(3,3))
	for x in range(3):
		for y in range(3):
			#Na = np.mat(arr[x]/dev[x])
			#Nb = np.mat(arr[y]/dev[y])
			Nahi = np.mat((hi_arr[x]-contr_arr)/devhi[x])
			Nbhi = np.mat((hi_arr[y]-contr_arr)/devhi[y])
			Nalo = -np.mat((lo_arr[x]-contr_arr)/devlo[x])
			Nblo = -np.mat((lo_arr[y]-contr_arr)/devlo[y])
			Na = 0.5* (Nahi+Nalo)
			Nb = 0.5* (Nbhi+Nblo)
			M = Na.T*Nb+Nb.T*Na
			Fisher[x,y] = 0.5 * trace(cov_inv*M)

	err = sqrt(mat(Fisher)**-1)[[0,1,2],[0,1,2]]/sqrt(20000./3.46**2)
	savetxt(errfile,err)

xiuyuan = np.array([0.0078,0.036,.0057])
print 'jia err',err


#def fishererr (params):
	#Fisher = np.ndarray(shape=(3,3))
	#for x in range(3):
		#for y in range(3):
			#a = params[x]
			#b = params[y]
			#M = (np.mat(fidu_mat[a]-fidu_mat[contr]).T*np.mat(fidu_mat[b]-fidu_mat[contr]) + np.mat(fidu_mat[b]-fidu_mat[contr]).T*np.mat(fidu_mat[a]-fidu_mat[contr]))/((paramvalue[x][-1]-paramvalue[x][1])*(paramvalue[y][-1]-paramvalue[y][1]))
			#Fisher[x,y] = 0.5 * trace(cov_inv*M)
	#Fisher = np.mat(Fisher)
	#F = sqrt(Fisher.I)
	#errOm, errw, errsi = F[0,0], F[1,1], F[2,2]
	#return errOm, errw, errsi

#errfw_file = fit_dir+'errfw_g%s_bins%s.ls'%(sigma,bins)
#errbw_file = fit_dir+'errbw_g%s_bins%s.ls'%(sigma,bins)

#if os.path.isfile(errfw_file):
	#errfw = genfromtxt(errfw_file)
#else:
	#errfw = fishererr (hi_params)
	#savetxt(errfw_file,errfw)
	
#if os.path.isfile(errbw_file):
	#errbw = genfromtxt(errbw_file)
#else:
	#errbw = fishererr (lo_params)
	#savetxt(errbw_file,errbw)
print 'Done-Done-Done-Done'