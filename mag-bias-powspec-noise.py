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

#s, bins = np.array(sys.argv[1:3],dtype=float)
s = float(sys.argv[1])
ngal_arcmin = float(sys.argv[2])
interpmethod = str(sys.argv[3]) #spline,fw,or bw

ngal=12.0*60**2*ngal_arcmin/2048.0**2 # galaxies per pixel
bins = 1000
z=1.0
gamma=0.15+0.035*z
noise_width=sqrt((gamma*gamma)/ngal)

s_arr = (0.0, 0.5, 1.0, 1.5, 2.0)
initguess = (0.27, -1.1, 0.76) # initial guesses for fitter
sigma = 6.976
#sigma = 9.865

low, high = -0.02, 0.19 # -0.15 , 0.25 kappa range
kmapnum=1000 # number of maps per cosmology
fidu_params = np.array([0.26, -1.0, 0.798])

iseed=0
ps_dir = '/direct/astro+u/jia/jia/powspec_noise/'
#ps_mat_dir = '/direct/astro+u/jia/magbias/ps_mat_noise/'
ps_mat_dir = '/direct/astro+u/jia/jia/ps_mat_noise/'
#ps_fit_dir = '/direct/astro+u/jia/magbias/ps_fit_noise/%s/'%(interpmethod) # folder to save fit files
ps_fit_dir = '/direct/astro+u/jia/jia/ps_fit_noise/%s/'%(interpmethod)
		
kmap_dir = '/direct/astro+u/jia/kmap/' # where all the map directories locate
plot_dir = '/direct/astro+u/jia/magbias/plot/' # folder to save plots
#ps_analytical_dir = '/direct/astro+u/jia/magbias/ps_fit_noise_analytical/%s/'%(interpmethod)

cosmo_all = os.listdir(kmap_dir) #a list of 7 cosmologies, suppose all the map folders are in kmap_dir

########### equations & functions that takes (bins, s, sigma) ########
#bin_size = (high - low) / bins
#hist_bins = np.linspace (low+bin_size/2, high-bin_size/2, bins) # center of each bin

biasnoise = lambda kmap, noise, s: kmap*(1+(5*s-2)*kmap)+noise

smooth = lambda kmap, sigma: snd.filters.gaussian_filter(kmap,sigma)

biasnoisesmooth = lambda kmap, noise, sigma, s: snd.filters.gaussian_filter(kmap*(1+(5*s-2)*kmap)+noise,sigma)/snd.filters.gaussian_filter(1+(5*s-2)*kmap,sigma)

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

def PowerSpectrum(kmap):
	F = fftpack.fftshift(fftpack.fft2(kmap))
	psd2D = np.abs(F)**2
	ell_arr,psd1D, psd1D2 = azimuthalAverage(psd2D)
	x=ell_arr
	norm = ((2*pi*sqrt(12.00)/360.0)**2)/(2048.0**2)**2
	powspec = x*(x+1)*psd1D/2/pi*norm
	return powspec

def powspec_noise(fitsfile, s=s, bias=0):
	if bias:
		psfile = ps_dir+fitsfile[-93:-4]+'_powspec_bias_s%s_ngal%s.ls'%(s,int(ngal_arcmin))
		 
	else:
		psfile = ps_dir+fitsfile[-93:-4]+'_powspec_true_ngal%s.ls'%(int(ngal_arcmin))

	if os.path.isfile(psfile):
		powspec = np.genfromtxt(psfile)
	else:
		hdulist = pyfits.open(fitsfile)#, ignore_missing_end = True)
		kmap = np.array(hdulist[0].data)
		noise = np.array(np.random.normal(0.0, noise_width, kmap.shape),dtype=float32)
		if bias:
			kmap_noise = biasnoise(kmap,noise,s)
		else:
			kmap_noise = kmap + noise
		powspec = PowerSpectrum(kmap_noise)		
		np.savetxt(psfile, powspec)
		print 'power spectrum for bias =',bool(bias),fitsfile[-93:-4]
	return powspec


def batch_kmap_noise(kmap_dir,cosmo_num, s=s, bias=0):
	cosmo = cosmo_all[cosmo_num]
	kmaps = os.listdir(kmap_dir+cosmo)
	if bias:
		file_matrix = ps_mat_dir+'ps_matrix_bias_'+str(kmapnum)+cosmo+'_s%s_ngal%s.ls'%(s,int(ngal_arcmin))
	else:
		file_matrix = ps_mat_dir+'ps_matrix_true_'+str(kmapnum)+cosmo+'_ngal%s.ls'%(int(ngal_arcmin))
	if os.path.isfile(file_matrix):
		ps_matrix = np.genfromtxt(file_matrix)
	else:
		ps_matrix = np.ndarray (shape = (len(kmaps), bins), dtype = float32)
		for n in range(len(kmaps)):
			seed(iseed+n) # not sure if I should seed here
			fitsfile = kmap_dir+cosmo+'/'+kmaps[n]
			powspec = powspec_noise(fitsfile, s=s, bias=bias)
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

true_mat_noise = batch_kmap_noise(kmap_dir,contr, bias=0)
bias_mat_noise = batch_kmap_noise(kmap_dir,contr, bias=1)

true_mat_noise_all = np.ndarray(shape=(len(cosmo_all), kmapnum, bins))
for i in np.arange(len(cosmo_all)):
#	print 'cosmology ' + str(i) + ' of ' + str(len(cosmo_all))
	true_mat_noise_all[i] = batch_kmap_noise(kmap_dir,i, bias=0)

avg_noise_all = average(true_mat_noise_all,axis=1)
fidu_mat = avg_noise_all

########### extrapolate Ni_arr for certain bins ##############
if interpmethod=='spline':
	interpOm = lambda ibin: interpolate.UnivariateSpline(Om,fidu_mat[[lo_Om,contr,hi_Om], ibin], k=2, s=0)
	interpw  = lambda ibin: interpolate.UnivariateSpline(w, fidu_mat[[lo_w ,contr,hi_w ], ibin], k=2, s=0)
	interpsi = lambda ibin: interpolate.UnivariateSpline(si,fidu_mat[[lo_si,contr,hi_si], ibin], k=2, s=0)
else:
	def linearextrap(xs, ys):
		x0,x1 = xs
		y0,y1 = ys
		def pointwise(x):
			if y1 == y0:
				return y0
			else:
				return y0+(y1-y0)/(x1-x0)*(x-x0)
		return pointwise
	if interpmethod=='fw':
		interpOm = lambda ibin: linearextrap(Om[1:],fidu_mat[[contr,hi_Om], ibin])
		interpw  = lambda ibin: linearextrap(w[1:], fidu_mat[[contr,hi_w ], ibin])
		interpsi = lambda ibin: linearextrap(si[1:],fidu_mat[[contr,hi_si], ibin])
	elif interpmethod=='bw':
		interpOm = lambda ibin: linearextrap(Om[:-1],fidu_mat[[lo_Om,contr], ibin])
		interpw  = lambda ibin: linearextrap(w[:-1], fidu_mat[[lo_w ,contr], ibin])
		interpsi = lambda ibin: linearextrap(si[:-1],fidu_mat[[lo_si,contr], ibin])


Ni = lambda ibin, iOm, iw, isi: interpOm(ibin)(iOm)+interpw(ibin)(iw)+ interpsi(ibin)(isi)-2.0*fidu_mat[contr,ibin]

Ni_arr = lambda iOm, iw, isi: np.array(map(Ni, np.arange(bins),array([iOm]*bins),array([iw]*bins),array([isi]*bins)))

################ from here, using only limited ell for fit #######
################ chi-square fitting ##############################

#cov_mat = np.cov(true_mat[contr],rowvar=0)

kmin=1.0
kmax=2048./2
binsize=(kmax-kmin)/1000.0
ell_arr = (np.arange(kmin,kmax,binsize)*360./sqrt(12.0))[:1000]

num=203
cov_mat = np.cov(true_mat_noise[:,:num],rowvar=0)
cov_inv = np.mat(cov_mat).I

ps_analytical_dir = '/direct/astro+u/jia/magbias/ps_fit_noise_ell%s/%s/'%(int(ell_arr[num-1]),interpmethod)

def chisq (initguess, ibmap):
	iOm, iw, isi = initguess
	del_Ni = np.mat((Ni_arr(iOm, iw, isi) - ibmap)[:num])
	chisquare = float64(del_Ni * cov_inv * del_Ni.T)
	return chisquare

	
####### 4/15/2013 added analytical fit#####


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

def analytical_fits(ps_mat, s=s, bins=bins, bias=1):
	if bias:
		fit_filename = ps_analytical_dir+'powspec_bias_analytical_s%s_ngal%s.ls'%(s,int(ngal_arcmin))
	else:
		fit_filename = ps_analytical_dir+'powspec_true_analytical_ngal%s.ls'%(int(ngal_arcmin))
	if os.path.isfile(fit_filename):
		fit_analytical = np.genfromtxt(fit_filename)
	else:	
		fit_analytical=np.ndarray(shape=(1000,4))
		for i in range(1000):
			print 'analytical fitting to map #',i,'bias=',bool(bias)
			ibmap=ps_mat[i]
			Y=np.mat(ibmap[:num]-fidu_mat[contr,:num])
			del_p=((X*cov_inv*X.T).I)*(X*cov_inv*Y.T)		
			initguess=np.squeeze(np.array(del_p.T))+fidu_params
			fit_analytical[i,0]=chisq(initguess, ibmap)
			fit_analytical[i,1:]=np.array(initguess)
		savetxt(fit_filename,fit_analytical)
	return fit_analytical

true_fit_analytical=analytical_fits(true_mat_noise, bias=0)
bias_fit_analytical=analytical_fits(bias_mat_noise, s=s, bias=1)


################# commented out 7/29/2013 after found problem in matrix inversion ##########
################# chi-square fitting ##############################
##cov_mat = np.cov(true_mat[contr],rowvar=0)
#cov_mat = np.cov(true_mat_noise,rowvar=0)
#cov_inv = np.mat(cov_mat).I


#def chisq (initguess, ibmap):
	#iOm, iw, isi = initguess
	#del_Ni = np.mat(Ni_arr(iOm, iw, isi) - ibmap)
	#chisquare = float64(del_Ni * cov_inv * del_Ni.T)
	#return chisquare


######## 4/15/2013 added analytical fit#####


#if interpmethod=='fw':
	#dNOm=fidu_mat[hi_Om]-fidu_mat[contr]
	#dNw =fidu_mat[hi_w ]-fidu_mat[contr]
	#dNsi=fidu_mat[hi_si]-fidu_mat[contr]
	#dp = np.array([0.03,0.2,0.052])
#elif interpmethod=='bw':
	#dNOm=fidu_mat[contr]-fidu_mat[lo_Om]
	#dNw =fidu_mat[contr]-fidu_mat[lo_w ]
	#dNsi=fidu_mat[contr]-fidu_mat[lo_si]
	#dp = np.array([0.03,0.2,0.048])
#dNdOm=dNOm/dp[0]
#dNdw =dNw/dp[1] 
#dNdsi=dNsi/dp[2]
#X=np.mat([dNdOm,dNdw,dNdsi])

#def analytical_fits(ps_mat, s=s, bins=bins, bias=1):
	#if bias:
		#fit_filename = ps_analytical_dir+'powspec_bias_analytical_s%s_ngal%s.ls'%(s,int(ngal_arcmin))
	#else:
		#fit_filename = ps_analytical_dir+'powspec_true_analytical_ngal%s.ls'%(int(ngal_arcmin))
	#if os.path.isfile(fit_filename):
		#fit_analytical = np.genfromtxt(fit_filename)
	#else:	
		#fit_analytical=np.ndarray(shape=(1000,4))
		#for i in range(1000):
			#print 'analytical fitting to map #',i,'bias=',bool(bias)
			#ibmap=ps_mat[i]
			#Y=np.mat(ibmap-fidu_mat[contr])
			#del_p=((X*cov_inv*X.T).I)*(X*cov_inv*Y.T)		
			#initguess=np.squeeze(np.array(del_p.T))+fidu_params
			#fit_analytical[i,0]=chisq(initguess, ibmap)
			#fit_analytical[i,1:]=np.array(initguess)
		#savetxt(fit_filename,fit_analytical)
	#return fit_analytical

#true_fit_analytical=analytical_fits(true_mat_noise, bias=0)
#bias_fit_analytical=analytical_fits(bias_mat_noise, s=s, bias=1)


print 'POWER SPECTRUM DONE-DONE-DONE-DONE-DONE-DONE'