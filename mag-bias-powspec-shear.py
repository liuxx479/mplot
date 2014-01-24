import numpy as np
from scipy import *
from scipy import interpolate
from pylab import *
import pyfits
import scipy.ndimage as snd
import os #math
from scipy import optimize
import sys

from scipy import fftpack

#inline
#python s ngal interpmethod 
########## constants & input ######

s = float(sys.argv[1])
#ngal_arcmin = float(sys.argv[2])

#ngal=12.0*60**2*ngal_arcmin/2048.0**2 # galaxies per pixel
#noise_width=sqrt((gamma*gamma)/ngal)
bins = 1000
z=1.0
gamma=0.15+0.035*z


s_arr = (0.0, 0.5, 1.0, 1.5, 2.0)
initguess = (0.27, -1.1, 0.76) # initial guesses for fitter
sigma = 6.976

kmapnum=1000 # number of maps per cosmology
fidu_params = np.array([0.26, -1.0, 0.798])

iseed=0
ps_dir = '/direct/astro+u/jia/jia/powspec_shear/'
ps_mat_dir = '/direct/astro+u/jia/jia/ps_mat_shear/'
		
kmap_dir = '/direct/astro+u/jia/kmap/' # where all the map directories locate
gmap_dir = '/direct/astro+u/jia/gmap/' 
#WL-shear2_m-512b240_Om0.230_Ol0.770_w-1.000_ns0.960_si0.798_4096xy_0716r_0029p_0100z_og.gre.fit
cosmo_all = os.listdir(kmap_dir) #a list of 7 cosmologies, suppose all the map folders are in kmap_dir

########### equations & functions that takes (bins, s, sigma) ########

biasgamma = lambda gmap, kmap, s: gmap*(1+(5*s-2)*kmap)

#biasnoise = lambda kmap, noise, s: kmap*(1+(5*s-2)*kmap)+noise
smooth = lambda kmap, sigma: snd.filters.gaussian_filter(kmap,sigma)
#biasnoisesmooth = lambda kmap, noise, sigma, s: snd.filters.gaussian_filter(kmap*(1+(5*s-2)*kmap)+noise,sigma)/snd.filters.gaussian_filter(1+(5*s-2)*kmap,sigma)

def azimuthalAverage(image, center=None):
	"""
	Calculate the azimuthally averaged radial profile.

	image - The 2D image
	center - The [x,y] pixel coordinates used as the center. The default is None, which then uses the center of the image (including fracitonal pixels).

	"""
	# Calculate the indices from the image
	y, x = np.indices(image.shape)

	if not center:
		center = np.array([(x.max()-x.min())/2.0, (x.max()-x.min())/2.0])

	r = np.hypot(x - center[0], y - center[1])#distance to center pixel, for each pixel

	# Get sorted radii
	ind = np.argsort(r.flat)
	r_sorted = r.flat[ind]
	i_sorted = image.flat[ind]
	ii_sorted = (image*r*(r+1)).flat[ind]#weighted by ell(ell+1)

	# find index that's corresponding to the lower edge of each bin
	kmin=1.0
	kmax=image.shape[0]/2.0
	#binsize=(kmax-kmin+0.001)/1000.0
	#r_arr=np.arange(kmin,kmax,binsize)
	hist_ind,ell_arr=np.histogram(r_sorted,bins=1000,range=(kmin,kmax+0.001)) # hist_ind: the number in each ell bins, sum them up is the index of lower edge of each bin
	
	ind2=hist_ind+4	# the first 4 is lower than kmin, so include in ell=1
	hist_ind[0]+=4 # per last line
	ind3=np.cumsum(hist_ind) # sum up to get index for lower edge of each bin
		
	# Cumulative sum to figure out sums for each radius bin
	csim = np.cumsum(i_sorted, dtype=float)	
	tbin=zeros(shape=1000)
	tbin[0]=csim[ind3[0]] # the first bin is the same as csim first bin
	tbin[1:]=csim[ind3[1:]]-csim[ind3[:-1]] # the other bins is csim[n] - csim[n-1], because they're cumulative
	tbin/=hist_ind # weighted by # of bins in that ell bin
	
	## normalize l(l+1) before binning, not used in the final code
	csim2 = np.cumsum(ii_sorted, dtype=float)	
	tbin2=zeros(shape=1000)
	tbin2[0]=csim2[ind3[0]]
	tbin2[1:]=csim2[ind3[1:]]-csim2[ind3[:-1]]
	tbin2/=hist_ind
	tbin2/=(ell_arr[:-1]*(ell_arr[:-1]+1))
	
	ell_arr*=360./sqrt(12) # normalized to our current map size
	
	return ell_arr[:-1], tbin, tbin2

def PowerSpectrum(img):
	F = fftpack.fftshift(fftpack.fft2(img))
	psd2D = np.abs(F)**2
	ell_arr,psd1D, psd1D2 = azimuthalAverage(psd2D)
	x=ell_arr
	norm = ((2*pi*sqrt(12.00)/360.0)**2)/(2048.0**2)**2
	powspec = x*(x+1)*psd1D/2/pi*norm
	return powspec

def powspec_gen(ikmap, igmap1, igmap2, s=s, bias=0):
	if bias:
		psfile = ps_dir+ikmap[-93:-4]+'_powspec_bias_s%s.ls'%(s)
		 
	else:
		psfile = ps_dir+ikmap[-93:-4]+'_powspec_true.ls'

	if os.path.isfile(psfile):
		powspec = np.genfromtxt(psfile)
	else:
		
		hdulist_g1 = pyfits.open(igmap1)
		gmap1 = np.array(hdulist_g1[0].data)
		hdulist_g2 = pyfits.open(igmap2)
		gmap2 = np.array(hdulist_g2[0].data)
		
		if bias:
			hdulist_k = pyfits.open(ikmap)
			kmap = np.array(hdulist_k[0].data)
			gmap1 = biasgamma(gmap1, kmap, s)
			gmap2 = biasgamma(gmap2, kmap, s)
		powspec = PowerSpectrum(smooth(gmap1,0))+PowerSpectrum(smooth(gmap2,0))		
		np.savetxt(psfile, powspec)
		print 'power spectrum for bias =',bool(bias),ikmap[-93:-4]
	return powspec


def batch_kmap_noise(kmap_dir,cosmo_num, s=s, bias=0):
	cosmo = cosmo_all[cosmo_num]
	kmaps = os.listdir(kmap_dir+cosmo)
	gmaps1 = os.listdir(gmap_dir+cosmo+'/shear1/')
	gmaps2 = os.listdir(gmap_dir+cosmo+'/shear2/')
	if bias:
		file_matrix = ps_mat_dir+'ps_matrix_bias_'+str(kmapnum)+cosmo+'_s%s.ls'%(s)
	else:
		file_matrix = ps_mat_dir+'ps_matrix_true_'+str(kmapnum)+cosmo+'.ls'
	if os.path.isfile(file_matrix):
		ps_matrix = np.genfromtxt(file_matrix)
	else:
		ps_matrix = np.ndarray (shape = (len(kmaps), bins), dtype = float32)
		for n in range(len(kmaps)):
			#seed(iseed+n) # not sure if I should seed here
			ikmap = kmap_dir+cosmo+'/'+kmaps[n]
			igmap1 = gmap_dir+cosmo+'/shear1/'+gmaps1[n]
			igmap2 = gmap_dir+cosmo+'/shear2/'+gmaps2[n]
			powspec = powspec_gen(ikmap, igmap1, igmap2, s=s, bias=bias)
			ps_matrix[n] = powspec
		print 'Now making powspec mat for cosmo ',cosmo_num, 'bias =',bool(bias)
		np.savetxt(file_matrix, ps_matrix)
	return ps_matrix
	

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

true_mat = batch_kmap_noise(kmap_dir,contr, bias=0)
bias_mat = batch_kmap_noise(kmap_dir,contr, bias=1)

#true_mat_noise_all = np.ndarray(shape=(len(cosmo_all), kmapnum, bins))
#for i in np.arange(len(cosmo_all)):
##	print 'cosmology ' + str(i) + ' of ' + str(len(cosmo_all))
	#true_mat_noise_all[i] = batch_kmap_noise(kmap_dir,i, bias=0)

savetxt('fidu_true_mat.ls',average(true_mat,axis=0))
savetxt('fidu_bias_mat_s%s.ls'%(s),average(bias_mat,axis=0))
print 'POWER SPECTRUM DONE-DONE-DONE-DONE-DONE-DONE'