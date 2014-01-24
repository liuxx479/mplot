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

#inline
#python s bins ngal interpmethod
########## constants & input ######

#s, bins = np.array(sys.argv[1:3],dtype=float)
s = float(sys.argv[1])
bins = int(sys.argv[2])
ngal_arcmin = float(sys.argv[3])
interpmethod = str(sys.argv[4]) #spline,fw,or bw
#ngal=ngal_arcmin/(2048.0/3.46/60.0)**2 # galaxies per pixel
ngal=12.0*60**2*ngal_arcmin/2048.0**2

z=1.0
gamma=0.15+0.035*z
noise_width=sqrt((gamma*gamma)/ngal)

initguess = (0.27, -1.1, 0.76) # initial guesses for fitter
sigma = 6.976
#sigma = 9.865

low, high = -0.02, 0.19 # -0.15 , 0.25 kappa range
kmapnum=1000 # number of maps per cosmology
fidu_params = np.array([0.26, -1.0, 0.798])

iseed=0

#analytical_dir = '/direct/astro+astronfs01/workarea/jia/fit_noise_sim45_lowk/' # folder to save fit files
analytical_dir = '/direct/astro+astronfs01/workarea/jia/fit_noise_sim45/'
peak_dir = '/direct/astro+astronfs01/workarea/jia/peaks_noise_sim45/'
peaks_mat_dir = '/direct/astro+astronfs01/workarea/jia/peaks_mat_noise_sim45/'
kmap_sim45_dir = '/direct/astro+u/jia/kmap_sim45/'

kmap_dir = '/direct/astro+u/jia/kmap/' # where all the map directories locate
peaks_mat_fidu_dir = '/direct/astro+u/jia/magbias/peaks_mat_noise/'
peak_dir_sim5 = '/direct/astro+u/jia/magbias/peaks_noise/'
cosmo_all = os.listdir(kmap_dir) #a list of 7 cosmologies, suppose all the map folders are in kmap_dir

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

########### equations & functions that takes (bins, s, sigma) ########

smooth = lambda kmap, sigma: snd.filters.gaussian_filter(kmap,sigma)

biasnoisesmooth = lambda kmap, noise, sigma, s: snd.filters.gaussian_filter(kmap*(1+(5*s-2)*kmap)+noise,sigma)/snd.filters.gaussian_filter(1+(5*s-2)*kmap,sigma)


def peaklist_noise(fitsfile, s=s, bins=bins, peak_dir=peak_dir, sigma=sigma, bias=0):
	if bias:
		peakfile = peak_dir+fitsfile[-93:-4]+'_bias_s%s_g%s_ngal%s.ls'%(s,sigma,int(ngal_arcmin))
		 
	else:
		peakfile = peak_dir+fitsfile[-93:-4]+'_true_g%s_ngal%s.ls'%(sigma,int(ngal_arcmin))

	if os.path.isfile(peakfile):
		peaks = np.genfromtxt(peakfile)
	else:
		hdulist = pyfits.open(fitsfile)#, ignore_missing_end = True)
		kmap = np.array(hdulist[0].data)
		noise = np.array(np.random.normal(0.0, noise_width, kmap.shape),dtype=float32)
		if bias:
			kmap_smooth = biasnoisesmooth(kmap,noise,sigma,s)
		else:
			kmap_noise = kmap + noise
			kmap_smooth = smooth (kmap_noise, sigma)
		peaks = findpeak(kmap_smooth)
		np.savetxt(peakfile, peaks)
		print 'counting peaks for bias=',bias,fitsfile[-93:-4]
	return peaks


def batch_kmap_noise(kmap_dir = kmap_sim45_dir, peaks_mat_dir = peaks_mat_dir, cosmo_num = contr, s=s, bins=bins, sigma=sigma, bias=0):
	cosmo = cosmo_all[cosmo_num]
	kmaps = os.listdir(kmap_dir+cosmo)
	if bias:
		file_matrix = peaks_mat_dir+'peak_matrix_bias_'+str(kmapnum)+cosmo+'_s%s_g%s_ngal%s_bins%s.ls'%(s,sigma,int(ngal_arcmin),bins)
	else:
		file_matrix = peaks_mat_dir+'peak_matrix_true_'+str(kmapnum)+cosmo+'_g%s_ngal%s_bins%s.ls'%(sigma,int(ngal_arcmin),bins)
	if os.path.isfile(file_matrix):
		peak_matrix = np.genfromtxt(file_matrix)
	else:
		peak_matrix = np.ndarray (shape = (len(kmaps), bins), dtype = float32)
		for n in range(len(kmaps)):
			seed(iseed+n) # not sure if I should seed here
			fitsfile = kmap_dir+cosmo+'/'+kmaps[n]
			if peaks_mat_dir == peaks_mat_fidu_dir:
				peaks = peaklist_noise(fitsfile, s=s, bins=bins, peak_dir=peak_dir_sim5, sigma=sigma, bias=bias)
			else:				
				peaks = peaklist_noise(fitsfile, s=s, bins=bins, sigma=sigma, bias=bias)
			peak_matrix[n] = np.histogram(peaks, range=(low,high), bins=bins)[0]
			print 'Now making histogram for fitsfile # bias=',str(bias),str(n)
		np.savetxt(file_matrix, peak_matrix)
	return peak_matrix
	


##### create true_mat & bias_mat for certain (bins, s, sigma) ####
##### in preparation for chi-sq fitting ####

###### here we get sim45 maps ####
true_mat_noise = batch_kmap_noise(kmap_dir = kmap_sim45_dir,cosmo_num = contr, bias=0)
bias_mat_noise = batch_kmap_noise(kmap_dir = kmap_sim45_dir,cosmo_num = contr, bias=1)

###### here we get model fidu map that's already calculated ###
true_mat_noise_all = np.ndarray(shape=(len(cosmo_all), kmapnum, bins))
for i in np.arange(len(cosmo_all)):
	#	print 'cosmology ' + str(i) + ' of ' + str(len(cosmo_all))
	true_mat_noise_all[i] = batch_kmap_noise(kmap_dir = kmap_dir, peaks_mat_dir = peaks_mat_fidu_dir, cosmo_num = i, bias=0)

fidu_mat = average(true_mat_noise_all,axis=1)
fidu_mat_sim45 = average(true_mat_noise,axis=0)


################ chi-square fitting ##############################
#cov_mat = np.cov(true_mat[contr],rowvar=0)
num=bins ### important change: use only up to kappa=0.085 when it's bins/2

cov_mat = np.cov(true_mat_noise_all[contr][:,:num],rowvar=0)
cov_inv = np.mat(cov_mat).I

if interpmethod=='fw':
	dNOm=fidu_mat[hi_Om]-fidu_mat[contr]
	dNw =fidu_mat[hi_w ]-fidu_mat[contr]
	dNsi=fidu_mat[hi_si]-fidu_mat[contr]
	dp = np.array([0.03,0.2,0.052])
elif interpmethod=='bw':
	dNOm=fidu_mat[contr]-fidu_mat[lo_Om]
	dNw =fidu_mat[contr]-fidu_mat[lo_w ]
	dNsi=fidu_mat[contr]-fidu_mat[lo_si]
	dp = np.array([0.03,0.2,0.048])
dNdOm=dNOm/dp[0]
dNdw =dNw/dp[1] 
dNdsi=dNsi/dp[2]
X=np.mat([dNdOm,dNdw,dNdsi])[:,:num]


def analytical_fits(peak_mat, s=s, bins=bins, sigma=sigma, bias=1):
	if bias:
		fit_filename = analytical_dir+'bias_analytical_s%s_g%s_ngal%s_bins%s.ls'%(s,sigma,int(ngal_arcmin),bins)
	else:
		fit_filename = analytical_dir+'true_analytical_g%s_ngal%s_bins%s.ls'%(sigma,int(ngal_arcmin),bins)
	if os.path.isfile(fit_filename):
		fit_analytical = np.genfromtxt(fit_filename)
	else:	
		fit_analytical=np.ndarray(shape=(1000,4))
		for i in range(1000):
			print 'analytical fitting to map #',i,'bias=',bias
			ibmap=peak_mat[i]
			Y=np.mat(ibmap[:num]-fidu_mat_sim45[:num])
			del_p=((X*cov_inv*X.T).I)*(X*cov_inv*Y.T)		
			initguess=np.squeeze(np.array(del_p.T))+fidu_params
			del_N=Y-del_p.T*X
			fit_analytical[i,0]=(1000.0-2.0-bins)/(1000-1.0)*del_N*cov_inv*del_N.T#chisq
			fit_analytical[i,1:]=np.array(initguess)
		savetxt(fit_filename,fit_analytical)
	return fit_analytical

true_fit_analytical=analytical_fits(true_mat_noise, s=s, bins=bins, sigma=sigma, bias=0)
bias_fit_analytical=analytical_fits(bias_mat_noise, s=s, bins=bins, sigma=sigma, bias=1)

###########################################

print 'DONE-DONE-DONE-DONE-DONE-DONE'