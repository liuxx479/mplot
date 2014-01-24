import numpy as np
from scipy import *
from scipy import interpolate, stats, optimize,linalg
from scipy.interpolate import UnivariateSpline as smoothz
from pylab import *
import os
from matplotlib.patches import Ellipse
import matplotlib.gridspec as gridspec
import pyfits
import scipy.ndimage as snd
from findpeak import findpeak
#import findpeak3 as fp2
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.transforms as transforms
from scipy import fftpack

s_arr = (0.0, 0.5, 1.0, 1.5, 2.0)
#s_arr = (0.1, 0.4, 1.5)
bins_arr =  (10, 20, 30, 40, 50, 60, 70, 80, 90, 100)
ngal_arr=[15,30,45]

sigma = 6.976 #9.865
low, high = -0.02, 0.19 # -0.15 , 0.2=5 kappa range
kmapnum=1000

s=0.0
ngal=30
ngal_arcmin = ngal
aaabins=30

smooth = lambda kmap, sigma: snd.filters.gaussian_filter(kmap,sigma)

interpmethod = 'fw'
kmap_dir = '/Users/jia/Documents/weaklensing/map2/'
peak_dir = '/Users/jia/Documents/weaklensing/magbias/peaks/'
plot_dir = '/Users/jia/Documents/weaklensing/magbias/plot/seed73000/%s/'%(interpmethod)
plot_dir2 = '/Users/jia/Documents/weaklensing/magbias/plot/'
fit_dir = '/Users/jia/Documents/weaklensing/magbias/fit/'
fit_noise_dir = '/Users/jia/Documents/weaklensing/magbias/fit_noise_seed73000/%s/'%(interpmethod)
fit_analytical_dir = '/Users/jia/Documents/weaklensing/magbias/fit_noise_sim45/fw_bins200/'
#fit_analytical_dir = '/Users/jia/Documents/weaklensing/magbias/fit_noise_analytical_seed0/%s/'%(interpmethod)

peaks_mat_noise_dir = '/Users/jia/Documents/weaklensing/magbias/peaks_mat_noise/'
peaks_mat_dir = peaks_mat_noise_dir

ps_mat_dir = '/Users/jia/Documents/weaklensing/magbias/ps_mat/'
#ps_fit_dir = '/Users/jia/Documents/weaklensing/magbias/ps_fit_covnoise/%s/'%(interpmethod)
#ps_mat_dir = '/Users/jia/Documents/weaklensing/magbias/ps_mat_noise/'

#ps_fit_dir = '/Users/jia/Documents/weaklensing/magbias/ps_fit_noise_ell21579/%s/'%(interpmethod)
ps_fit_dir = '/Users/jia/Documents/weaklensing/magbias/ps_fit_noise_sim45/fw/'

label_str=['omega_m','w','sigma_8']
label_latex=['\Omega_m','w','\sigma_8']
fidu_params = np.array([0.26, -1.0, 0.798])

cosmo_all = os.listdir(kmap_dir)

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

#labels=np.ndarray(shape=(7,),dtype=str)
labels=['']*7
labels[hi_Om] = 'High-$\Omega_m$'
labels[lo_Om] = 'Low-$\Omega_m$'
labels[hi_w ] = 'High-$w$'
labels[lo_w ] = 'Low-$w$'
labels[hi_si] = 'High-$\sigma_8$'
labels[lo_si] = 'Low-$\sigma_8$'
labels[contr] = 'Fiducial'
############## read & organize data tables ############

#Bells = np.genfromtxt(fit_dir+'bias_ells_brute.ls')
#Tells = np.genfromtxt(fit_dir+'true_ells_brute.ls')
# s, bins, x0, y0, z0, dx, dy, dz, dth, dph, d, e, f, g, h, x, y, z, u, v, Ns
#ChiCont = np.genfromtxt(fit_dir+'chisqcontour_true.ls')
# a.fun, dx, dy, dz, dth, dph, bins


########## uncomment begin#############
#Bfit = np.ndarray(shape=(5,10,1000,4)) # 4 values, chi^2, 3 params
#i=0
#for s in s_arr:
	#j=0
	#for bins in bins_arr:
		#f = np.genfromtxt(fit_dir+'bias_fit_s%s_g9.865_bins%s.ls'%(s,bins))
		#Bfit[i,j]=f
		#j += 1
	#i += 1

#Tfit = np.ndarray(shape=(10,1000,4))
#i = 0
#for bins in bins_arr:
	#f = np.genfromtxt(fit_noise_dir+'true_fit_g6.976_ngal30_bins%s.ls'%(bins))
	#Tfit[i]=f
	#i += 1

#Bfit_ngal15 = np.ndarray(shape=(5,10,1000,4))	
#Bfit_ngal30 = np.ndarray(shape=(5,10,1000,4))
#Bfit_ngal45 = np.ndarray(shape=(5,10,1000,4))
#i=0
#for s in s_arr:
	#j=0
	#for bins in bins_arr:
		#f15 = np.genfromtxt(fit_noise_dir+'bias_fit_s%s_g6.976_ngal15_bins%s.ls'%(s,bins))
		#Bfit_ngal15[i,j]=f15
		#f30 = np.genfromtxt(fit_noise_dir+'bias_fit_s%s_g6.976_ngal30_bins%s.ls'%(s,bins))
		#Bfit_ngal30[i,j]=f30
		#f45 = np.genfromtxt(fit_noise_dir+'bias_fit_s%s_g6.976_ngal45_bins%s.ls'%(s,bins))
		#Bfit_ngal45[i,j]=f45
		#j += 1
	#i += 1

###### turn on to get fit_noise map fits
Tfit_ngal15 = np.ndarray(shape=(10,1000,4))	
Tfit_ngal30 = np.ndarray(shape=(10,1000,4))
Tfit_ngal45 = np.ndarray(shape=(10,1000,4))
i = 0

for bins in bins_arr:
	f15 = np.genfromtxt(fit_noise_dir+'true_fit_g6.976_ngal15_bins%s.ls'%(bins))
	f30 = np.genfromtxt(fit_noise_dir+'true_fit_g6.976_ngal30_bins%s.ls'%(bins))
	f45 = np.genfromtxt(fit_noise_dir+'true_fit_g6.976_ngal45_bins%s.ls'%(bins))
	Tfit_ngal15[i]=f15
	Tfit_ngal30[i]=f30
	Tfit_ngal45[i]=f45
	i += 1
	
#####end: turn on to get fit_noise map fits
########## uncomment end#############	



########plots for paper########

bias_changes_cosmo = 0
bias_changes_cosmo_powspec = 0
bias_changes_cosmo_powspec_noise = 0
dN=0 # test abs value of dN, not dN/N
dN_mask=0 # test masked dN

delw_vs_s = 0
delw_vs_s_powspec = 0
delw_vs_s_powspec_noise = 0
delw_vs_s_powspec_noise_old = 0
del_w_vs_s_fw_bw = 0

ells_fidu_s01 = 0
ells_fidu_s0108 = 0 # plot s=0.2 and s=0.8
ells_fidu_s01b = 0 # different alpha for powspec and peaks
destroy_high_peaks_singlemap = 0


Schmidt09fig1 = 0
######## size bias plot #######

DEEP_test = 0

########plots for testing #####
peaks_err_vs_bins = 0
test_num_bins_cov = 0
peaks_err_vs_bins_old = 0
sigmak = 0
shear_fft = 0
shear_conv_fft = 0
bias_changes_cosmo_powspec_noise_4pt = 0
bias_changes_cosmo_powspec_shear = 0
test_actual_fit_powspec = 0
test_dN_noise_nonoise = 0
find_CFHT_err = 0
test_CFHT_err = 0
powspec_test = 0

find_alpha = 0 #alpha for si*Om^alpha
analytical_chisqmin=0#test zoltan's derivation
fw_bw_true_scatter=0

compare_04_true=0
destroy_high_peaks=0
#destroy_high_peaks_singlemap=0
cosmo_truemap=0
individual_components = 0

origin_vs_bins = 0
origin_vs_bins_noise = 0

brute_vs_chisqcontour_per_bin = 0
brute_vs_s_per_bin = 0
density_sum = 0
chisq_density_comparison = 0
chisq_density_comparison_noise = 0 # main plotting routine
chisq_density_comparison_noise_bias = 0 # test what's wrong with bias maps
margin2params_noise = 0

err_fisher_vs_fit_ngal30 = 0
chisq_density_comparison_xy = 0
test_kinks = 0
density_contour = 0
margin2params = 0
param_width = 0
fisher_err = 0
smoothcompxy = 0
all_errs = 0
## test errors
comp_xiuyuan_distribution = 0
Fisher_vs_Chisq_test = 0
noise_jia_vs_xy = 0
seeded100noises=0

seed(10)
mcolors = rand(100)

RotM = lambda theta, phi: np.mat([[cos(theta),-cos(phi)*sin(theta), sin(phi)*sin(theta)], [sin (theta), cos(phi)*cos(theta), -sin(phi)*cos(theta)],[0,sin(phi),cos(phi)]])

import scipy.ndimage as snd
smooth = lambda image, scale: snd.filters.gaussian_filter(image,scale)

biasnoisesmooth = lambda kmap, noise, sigma, s: snd.filters.gaussian_filter(kmap*(1+(5*s-2)*kmap)+noise,sigma)/snd.filters.gaussian_filter(1+(5*s-2)*kmap,sigma)

bias = lambda kmap, sigma, s: snd.filters.gaussian_filter(kmap*(1+(5*s-2)*kmap),sigma)/snd.filters.gaussian_filter(1+(5*s-2)*kmap,sigma)

bias_nosmooth = lambda kmap, s:kmap*(1+(5*s-2)*kmap)
########## get all the peak distributions########
bins=aaabins
def filename_to_map (fitsfile):
	hdulist = pyfits.open(fitsfile)
	kmap = np.array(hdulist[0].data)
	kmap = smooth(kmap,0)#otherwise will show warning "ValueError: type >f4 is not supported"
	return kmap

def batch_kmap_noise(kmap_dir,cosmo_num, s=s, bins=bins, sigma=sigma, bias=0,ngal=ngal):
	cosmo = cosmo_all[cosmo_num]
	kmaps = os.listdir(kmap_dir+cosmo)
	if bias:
		file_matrix = peaks_mat_dir+'peak_matrix_bias_'+str(kmapnum)+cosmo+'_s%s_g%s_ngal%s_bins%s.ls'%(s,sigma,ngal,bins)
	else:
		file_matrix = peaks_mat_dir+'peak_matrix_true_'+str(kmapnum)+cosmo+'_g%s_ngal%s_bins%s.ls'%(sigma,ngal,bins)
	peak_matrix = np.genfromtxt(file_matrix)
	return peak_matrix

def batch_kmap_powspec(kmap_dir,cosmo_num, ps_mat_dir=ps_mat_dir, s=s, bias=0, ngal=None):
	cosmo = cosmo_all[cosmo_num]
	kmaps = os.listdir(kmap_dir+cosmo)
	if ngal:
		if bias:
			file_matrix = ps_mat_dir+'ps_matrix_bias_'+str(kmapnum)+cosmo+'_s%s_ngal%s.ls'%(s,ngal)
		else:
			file_matrix = ps_mat_dir+'ps_matrix_true_'+str(kmapnum)+cosmo+'_ngal%s.ls'%(ngal)
	else:
		if bias:
			file_matrix = ps_mat_dir+'ps_matrix_bias_'+str(kmapnum)+cosmo+'_s%s.ls'%(s)
		else:
			file_matrix = ps_mat_dir+'ps_matrix_true_'+str(kmapnum)+cosmo+'.ls'
	ps_matrix = np.genfromtxt(file_matrix)
	return ps_matrix

def interp_distribution(interpmethod,fidu_mat,iOm, iw, isi,bins=bins):
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
	
	return Ni_arr(iOm, iw, isi)



##########################################################

def azimuthalAverage(image, center=None):
	"""
	Calculate the azimuthally averaged radial profile.

	image - The 2D image
	center - The [x,y] pixel coordinates used as the center. The default is 
		None, which then uses the center of the image (including 
		fracitonal pixels).

	"""
	# Calculate the indices from the image
	y, x = np.indices(image.shape)

	if not center:
		center = np.array([(x.max()-x.min())/2.0, (x.max()-x.min())/2.0])

	r = np.hypot(x - center[0], y - center[1])

	# Get sorted radii
	ind = np.argsort(r.flat)
	r_sorted = r.flat[ind]
	i_sorted = image.flat[ind]
	ii_sorted = (image*r*(r+1)).flat[ind]


	# jia change, bin differently
	kmin=1.0
	kmax=image.shape[0]/2.0
	binsize=(kmax-kmin+0.001)/1000.0
	#r_arr=np.arange(kmin,kmap,binsize)
	#r_arr=np.linspace(kmin,kmax,1000)
	#r_arr=np.linspace(kmin,kmax+binsize,1000)
	#ind2=np.histogram(r_sorted,bins=r_arr)[0]
	hist_ind,ell_arr=np.histogram(r_sorted,bins=1000,range=(kmin,kmax+0.001))
	
	ind2=hist_ind+4	
	hist_ind[0]+=4
	ind3=np.cumsum(hist_ind)
		
	# Cumulative sum to figure out sums for each radius bin
	csim = np.cumsum(i_sorted, dtype=float)	
	tbin=zeros(shape=1000)
	tbin[0]=csim[ind3[0]]
	tbin[1:]=csim[ind3[1:]]-csim[ind3[:-1]]
	tbin/=hist_ind
	
	## normalize l(l+1) before binning
	csim2 = np.cumsum(ii_sorted, dtype=float)	
	tbin2=zeros(shape=1000)
	tbin2[0]=csim2[ind3[0]]
	tbin2[1:]=csim2[ind3[1:]]-csim2[ind3[:-1]]
	tbin2/=hist_ind
	tbin2/=(ell_arr[:-1]*(ell_arr[:-1]+1))
	
	
	ell_arr*=360./sqrt(12)
	
	return ell_arr[:-1], tbin, tbin2

############### all the plotting routines ################

if dN_mask:
	
	#bins=200
	for bins in (10,20,50,100,150,200,500):
		f=figure(figsize=(12,4))
		SN_arr=np.linspace(low,high,bins)
		fn = lambda Om, w, si, imaks:average(genfromtxt('/Users/jia/weaklensing/magbias/peaks_mat_noise_mask/peak_matrix_true_1000m-512b240_Om%s0_Ol%s0_w%s00_ns0.960_si%.3f_g6.976_ngal30_bins%s_mask%s.ls'%(Om,1-float(Om), w, si, bins, imaks)),axis=0)
		
		params2=(hi_Om,hi_w,hi_si)
		for imask in (0.284, 0.531, 0.7):
			fn_fidu = fn(0.26,-1.0,0.798,imask)		
			for i in range(3):
				Om, w, si = params[params2[i]]
				X=fn(Om, w, si,imask)-fn_fidu
				ax = f.add_subplot(1,3,i+1)
				ax.plot(SN_arr,X/(1-imask),label=str(imask))
				ax.set_title(r'$%s$ bins%s'%(label_latex[i],bins))
				ax.set_xlabel('$\kappa$')
				leg=ax.legend(loc=0, ncol=1, labelspacing=.2, prop={'size':10})
				leg.get_frame().set_visible(False)
			
		savefig(plot_dir2+'dN_norm_mask_bins%s.jpg'%(bins))
			#savefig(plot_dir2+'dN_mask_bins%s.jpg'%(bins))
		close()
	
if Schmidt09fig1:

	kmin=1.0
	kmax=2048./2
	binsize=(kmax-kmin)/1000.0
	ell_arr = (np.arange(kmin,kmax,binsize)*360./sqrt(12.0))[:1000]
	bins=1000
	
	########## grab relevant peak distributions #####
	#powspec_mat=batch_kmap_powspec(kmap_dir,contr,ps_mat_dir=ps_mat_dir,bias=0, ngal=ngal)
	#yerr=np.std(powspec_mat,axis=0)
	fidu_truek_mat=np.genfromtxt('/Users/jia/Documents/weaklensing/magbias/ps_mat/fidu_true_mat.ls')[contr]
	fidu_true_mat=np.genfromtxt('/Users/jia/Documents/weaklensing/magbias/ps_mat_shear/fidu_true_mat.ls')
			
	fidu_bias_mat = np.genfromtxt('/Users/jia/Documents/weaklensing/magbias/ps_mat_shear/fidu_bias_mat_s0.8.ls')

	
	########## actual plots #######
	f = figure(figsize=(6,6))
	ax=f.add_subplot(111)
	ax.plot(ell_arr,fidu_bias_mat/fidu_true_mat-1.0,linewidth=2,label='Bias $s$ = 0.8')
	xlabel(r'$\ell$',fontsize=16)
	ax.set_ylabel(r'${\Delta}P/P$',fontsize=16)	
	ax.set_yscale('log')	
	ax.set_xscale('log')
	
	ax.set_xlim(60,2e4)
	ax.set_ylim(0.005,0.2)
	#show()	
	savefig(plot_dir2+'Schmidt09fig1.pdf')
	close()
	
if bias_changes_cosmo:
	bins=100
	bin_size = (high - low) / bins
	hist_bins = np.linspace (low+bin_size/2, high-bin_size/2, bins) # center of each bin
	for ngal in (30,):#ngal_arr:
		########## grab relevant peak distributions #####
		true_mat_noise_all = np.ndarray(shape=(len(cosmo_all), kmapnum, bins))

		for i in np.arange(len(cosmo_all)):
			true_mat_noise_all[i] = batch_kmap_noise(kmap_dir,i,bins=bins, bias=0,ngal=ngal)

		avg_noise_all = average(true_mat_noise_all,axis=1)
		fidu_mat = avg_noise_all

		s_arr2=(0.2, 0.4, 0.8)
		
		j=0
		bias_mat = np.ndarray(shape=(len(s_arr2),1000,bins))
		for s in s_arr2:	
			f = batch_kmap_noise(kmap_dir,contr,s=s,bias=1,bins=bins,ngal=ngal)
			bias_mat[j]=f
			j += 1	
		fidu_bias_mat=average(bias_mat,axis=1)
		
		
		########## actual plots #######
		markers=('k-','b--','r:','c*','g+','md','yd','2','3','1','g-.')
		markersizes=(5,5,5,8,5,5,5,5,5,5,5)
		fname=0
		for j in ((hi_si,lo_si),(hi_w,lo_w),(hi_Om,lo_Om)):
			f = figure(figsize=(8,6))
		
			gs = gridspec.GridSpec(2,1,height_ratios=[3,1])#width_ratios=[1,1]
			ax1=f.add_subplot(gs[0])
			ax2=f.add_subplot(gs[1],sharex=ax1)
			
			
			ax1.plot(hist_bins,fidu_mat[contr],markers[0],linewidth=1.5, label=labels[contr]+' (%i)'%(sum(fidu_mat[contr])))
			### add error bars here ###
			yerr = np.std(true_mat_noise_all[contr],axis=0)
			
			
			ax2.plot(hist_bins,zeros(len(hist_bins)),'k-',linewidth=1.5)
			
			n=1
			for i in j:
				ax1.plot(hist_bins,fidu_mat[i],markers[n],linewidth=2, label=labels[i]+' (%i)'%(sum(fidu_mat[i])))
				ax2.plot(hist_bins,fidu_mat[i]/fidu_mat[contr]-1.0,markers[n],linewidth=2, label=labels[i]+' (%i)'%(sum(fidu_mat[i])))
				n+=1
		
			for i in range(len(s_arr2)):#range(5):
				ax1.plot(hist_bins,fidu_bias_mat[i],markers[n],linewidth=2, label='Bias $s$ = %s (%i)'%(s_arr2[i],sum(fidu_bias_mat[i])),markersize=markersizes[n])
				ax2.plot(hist_bins,fidu_bias_mat[i]/fidu_mat[contr]-1.0,markers[n],linewidth=2,markersize=markersizes[n],label='Bias $s$ = %s (%i)'%(s_arr2[i],sum(fidu_bias_mat[i])))
				
				n+=1
			
			ebp = 16
			ax1.errorbar(hist_bins[ebp],fidu_mat[contr,ebp],yerr[ebp],fmt='k-',linewidth=1.5)	
			#ax.set_yscale('log',basey=5)
			plt.setp(ax1.get_xticklabels(), visible=False)
			plt.subplots_adjust(hspace=0.0)
			
			
			xlabel('convergence $\kappa$',fontsize=16)
			ax1.set_ylabel('peak counts $N$($\kappa$)',fontsize=16)
			ax2.set_ylabel(r'${\Delta}N/N$',fontsize=16)
			ax1.set_ylim(0,220)
			ax2.set_ylim(-0.15,0.15)
			#ax2.xaxis.set_label_coords(0.5, -0.35)

			xlim(-0.02,0.10)
			
			### get rid of error bars in legend
			handles, newlabels = ax2.get_legend_handles_labels()
			
			leg=ax1.legend(handles, newlabels, loc=0, ncol=1, labelspacing=.2, prop={'size':16})
			leg.get_frame().set_visible(False)
			
			#show()
			savefig(plot_dir+label_str[::-1][fname]+'ngal%s.pdf'%(ngal))
			
			close()
			fname+=1
if dN:
	bins=100
	bin_size = (high - low) / bins
	hist_bins = np.linspace (low+bin_size/2, high-bin_size/2, bins) # center of each bin
	for ngal in (30,):#ngal_arr:
		########## grab relevant peak distributions #####
		true_mat_noise_all = np.ndarray(shape=(len(cosmo_all), kmapnum, bins))

		for i in np.arange(len(cosmo_all)):
			true_mat_noise_all[i] = batch_kmap_noise(kmap_dir,i,bins=bins, bias=0,ngal=ngal)

		avg_noise_all = average(true_mat_noise_all,axis=1)
		fidu_mat = avg_noise_all

		s_arr2=(0.2, 0.4, 0.8)
		
		j=0
		bias_mat = np.ndarray(shape=(len(s_arr2),1000,bins))
		for s in s_arr2:	
			f = batch_kmap_noise(kmap_dir,contr,s=s,bias=1,bins=bins,ngal=ngal)
			bias_mat[j]=f
			j += 1	
		fidu_bias_mat=average(bias_mat,axis=1)
		
		
		########## actual plots #######
		markers=('k-','b--','r:','c*','g+','md','yd','2','3','1','g-.')
		markersizes=(5,5,5,8,5,5,5,5,5,5,5)
		fname=0
		
		for j in ((hi_si,lo_si),(hi_w,lo_w),(hi_Om,lo_Om)):
			f = figure(figsize=(8,6))
		
			gs = gridspec.GridSpec(2,1,height_ratios=[1,1])#width_ratios=[1,1]
			ax1=f.add_subplot(gs[0])
			ax2=f.add_subplot(gs[1],sharex=ax1)
			
			
			ax1.plot(hist_bins,fidu_mat[contr],markers[0],linewidth=1.5, label=labels[contr]+' (%i)'%(sum(fidu_mat[contr])))
			### add error bars here ###
			yerr = np.std(true_mat_noise_all[contr],axis=0)
			
			
			ax2.plot(hist_bins,zeros(len(hist_bins)),'k-',linewidth=1.5)
			
			n=1
			for i in j:
				ax1.plot(hist_bins,fidu_mat[i],markers[n],linewidth=2, label=labels[i]+' (%i)'%(sum(fidu_mat[i])))
				ax2.plot(hist_bins,fidu_mat[i]-fidu_mat[contr],markers[n],linewidth=2, label=labels[i]+' (%i)'%(sum(fidu_mat[i])))
				n+=1
		
			for i in range(len(s_arr2)):#range(5):
				ax1.plot(hist_bins,fidu_bias_mat[i],markers[n],linewidth=2, label='Bias $s$ = %s (%i)'%(s_arr2[i],sum(fidu_bias_mat[i])),markersize=markersizes[n])
				ax2.plot(hist_bins,fidu_bias_mat[i]-fidu_mat[contr],markers[n],linewidth=2,markersize=markersizes[n],label='Bias $s$ = %s (%i)'%(s_arr2[i],sum(fidu_bias_mat[i])))
				
				n+=1
				
			
			ebp = 16
			ax1.errorbar(hist_bins[ebp],fidu_mat[contr,ebp],yerr[ebp],fmt='k-',linewidth=1.5)	
			#ax.set_yscale('log',basey=5)
			plt.setp(ax1.get_xticklabels(), visible=False)
			plt.subplots_adjust(hspace=0.0)
			
			
			xlabel('convergence $\kappa$',fontsize=16)
			ax1.set_ylabel('peak counts $N$($\kappa$)',fontsize=16)
			ax2.set_ylabel(r'${\Delta}N$',fontsize=16)
			ax1.set_ylim(0,220)
			#ax2.set_ylim(-0.15,0.15)
			#ax2.xaxis.set_label_coords(0.5, -0.35)

			xlim(-0.02,0.10)
			
			### get rid of error bars in legend
			handles, newlabels = ax2.get_legend_handles_labels()
			
			leg=ax1.legend(handles, newlabels, loc=0, ncol=1, labelspacing=.2, prop={'size':16})
			leg.get_frame().set_visible(False)
			
			#show()
			#savefig(plot_dir+label_str[::-1][fname]+'ngal%s.pdf'%(ngal))
			
			savefig(plot_dir2+label_str[::-1][fname]+'ngal%s_dN.jpg'%(ngal))
			close()
			fname+=1
			
if bias_changes_cosmo_powspec:
	
	kmin=1.0
	kmax=2048./2
	binsize=(kmax-kmin)/1000.0
	ell_arr = (np.arange(kmin,kmax+0.001,binsize)*360./sqrt(12.0))[:1000]
	bins=1000
	########## grab relevant peak distributions #####
	powspec_mat=batch_kmap_powspec(kmap_dir,contr,bias=0)
	yerr=np.std(powspec_mat,axis=0)
	
	fidu_mat_fn=ps_mat_dir+'fidu_true_mat.ls'
	if os.path.isfile(fidu_mat_fn):
		fidu_mat=np.genfromtxt(fidu_mat_fn)
	else:
		true_mat = np.ndarray(shape=(len(cosmo_all), kmapnum, bins))

		for i in np.arange(len(cosmo_all)):
			true_mat[i] = batch_kmap_powspec(kmap_dir,i, s=s, bias=0)

		fidu_mat = average(true_mat,axis=1)
		savetxt(fidu_mat_fn,fidu_mat)
		

	s_arr2=(0.2, 0.4, 0.8)		
		
	fidu_bias_fn=ps_mat_dir+'fidu_bias_mat_08.ls'
	#fidu_bias_fn=ps_mat_dir+'fidu_bias_mat.ls' #for s_arr2=(0.1, 0.4, 0.99)
	if os.path.isfile(fidu_bias_fn):
		fidu_bias_mat=np.genfromtxt(fidu_bias_fn)
	else:		
		print 'making fidu_bias_fn'
		bias_mat = np.ndarray(shape=(len(s_arr2),1000, bins))
		j=0
		for s in s_arr2:	
			f = batch_kmap_powspec(kmap_dir,contr, s=s, bias=1)
			bias_mat[j]=f
			j += 1	
		fidu_bias_mat=average(bias_mat,axis=1)
		savetxt(fidu_bias_fn,fidu_bias_mat)
	
	
	########## actual plots #######
	markers=('k-','b--','r-.','c*:','g+','md:','yd','2','3','1','g-.')
#	markers=('k-','b-','r-','c--','g*','m--','yd','2','3','1','g-.')
	markersizes=(5,5,5,8,5,5,5,5,5,5,5)
	lw_arr=(2,1.5,1,2,1,1.5)
	fname=0
	noise_powspec=np.genfromtxt('/Users/jia/Documents/weaklensing/magbias/ps_mat_noise/fidu_true_mat_ngal30.ls')[contr]-fidu_mat[contr]
	for j in ((hi_si,lo_si),(hi_w,lo_w),(hi_Om,lo_Om)):
		f = figure(figsize=(8,6))
	
		gs = gridspec.GridSpec(2,1,height_ratios=[3,1])#width_ratios=[1,1]
		ax1=f.add_subplot(gs[0])
		ax2=f.add_subplot(gs[1],sharex=ax1)
		
		n=1
		
		ax1.plot(ell_arr,noise_powspec,'k:',linewidth=2)# add noise in the background
		ax1.plot(ell_arr,fidu_mat[contr],markers[0],linewidth=2, label=labels[contr])
		ax2.plot(ell_arr,zeros(len(ell_arr)),'k-',linewidth=2,label=labels[contr])
		logind=np.array(np.rint(np.logspace(0,3,20))-1,dtype=int)
		for i in j:
			ax1.plot(ell_arr,fidu_mat[i],markers[n],linewidth=1.5, label=labels[i])
			ax2.plot(ell_arr,fidu_mat[i]/fidu_mat[contr]-1.0,markers[n],linewidth=1.5, label=labels[i])
			n+=1
	
		for i in range(len(s_arr2)):#range(5):
		
			ax1.plot(ell_arr[logind],fidu_bias_mat[i][logind],markers[n],linewidth=2, label='Bias $s$ = %s'%(s_arr2[i]),markersize=markersizes[n])
			ax2.plot(ell_arr[logind],(fidu_bias_mat[i]/fidu_mat[contr]-1.0)[logind],markers[n],linewidth=2,markersize=markersizes[n],label='Bias $s$ = %s'%(s_arr2[i]))
			n+=1
		
		ebp=logind
		ax1.errorbar(ell_arr[ebp],fidu_mat[contr,ebp],yerr[ebp],fmt='k-',linewidth=1.5)
		
		#ax.set_yscale('log',basey=5)
		plt.setp(ax1.get_xticklabels(), visible=False)
		plt.subplots_adjust(hspace=0.0)
		#leg=ax1.legend(loc=0, ncol=1, labelspacing=.2, prop={'size':16})
		#leg.get_frame().set_visible(False)
		
		xlabel(r'$\ell$',fontsize=16)
		ax1.set_ylabel(r'$\ell(\ell+1)P(\ell)/2\pi$',fontsize=16)
		ax2.set_ylabel(r'${\Delta}P/P$',fontsize=16)
		
		ax1.set_ylim(1e-7,1e-3)
		ax1.set_yscale('log')
		ax1.set_xscale('log')
		ax2.set_xscale('log')
		ax2.set_xlim(100,1e5)
		ax2.set_ylim(-0.3,0.6)
		ax2.set_yticks(np.arange(-0.2,0.55,0.2))
		
		handles, newlabels = ax2.get_legend_handles_labels()
		leg=ax1.legend(handles, newlabels, loc=0, ncol=1, labelspacing=.2, prop={'size':16})
		leg.get_frame().set_visible(False)
		#show()
		
		# for noiseless, only covmat has noise
		savefig(plot_dir2+'powspec'+label_str[::-1][fname]+'.pdf')
		close()	
		
		fname+=1

if bias_changes_cosmo_powspec_noise:
	ngal=30 ###
	kmin=1.0
	kmax=2048./2
	binsize=(kmax-kmin)/1000.0
	ell_arr = (np.arange(kmin,kmax,binsize)*360./sqrt(12.0))[:1000]
	bins=1000
	########## grab relevant peak distributions #####
	powspec_mat=batch_kmap_powspec(kmap_dir,contr,bias=0, ngal=ngal)
	yerr=np.std(powspec_mat,axis=0)
	
	fidu_mat_fn=ps_mat_dir+'fidu_true_mat_ngal%s.ls'%(ngal)
	if os.path.isfile(fidu_mat_fn):
		fidu_mat=np.genfromtxt(fidu_mat_fn)
	else:
		true_mat = np.ndarray(shape=(len(cosmo_all), kmapnum, bins))

		for i in np.arange(len(cosmo_all)):
			true_mat[i] = batch_kmap_powspec(kmap_dir,i, s=s, bias=0, ngal=ngal)

		fidu_mat = average(true_mat,axis=1)
		savetxt(fidu_mat_fn,fidu_mat)
		

	s_arr2=(0.1, 0.4, 0.8)		
		
	fidu_bias_fn=ps_mat_dir+'fidu_bias_mat_08_ngal%s.ls'%(ngal)
	#fidu_bias_fn=ps_mat_dir+'fidu_bias_mat.ls' #for s_arr2=(0.1, 0.4, 0.99)
	if os.path.isfile(fidu_bias_fn):
		fidu_bias_mat=np.genfromtxt(fidu_bias_fn)
	else:		
		bias_mat = np.ndarray(shape=(len(s_arr2),1000, bins))
		j=0
		for s in s_arr2:	
			f = batch_kmap_powspec(kmap_dir,contr, s=s, bias=1,ngal=ngal)
			bias_mat[j]=f
			j += 1	
		fidu_bias_mat=average(bias_mat,axis=1)
		savetxt(fidu_bias_fn,fidu_bias_mat)
	
	
	########## actual plots #######
	markers=('k-','b--','r-.','c*:','g+','md:','yd','2','3','1','g-.')
#	markers=('k-','b-','r-','c--','g*','m--','yd','2','3','1','g-.')
	markersizes=(5,5,5,8,5,5,5,5,5,5,5)
	lw_arr=(2,1.5,1,2,1,1.5)
	fname=0
	for j in ((hi_si,lo_si),(hi_w,lo_w),(hi_Om,lo_Om)):
		f = figure(figsize=(8,6))
	
		gs = gridspec.GridSpec(2,1,height_ratios=[3,1])#width_ratios=[1,1]
		ax1=f.add_subplot(gs[0])
		ax2=f.add_subplot(gs[1],sharex=ax1)
		
		n=1
		ax1.plot(ell_arr,fidu_mat[contr],markers[0],linewidth=2, label=labels[contr])
		ax2.plot(ell_arr,zeros(len(ell_arr)),'k-',linewidth=2,label=labels[contr])
		logind=np.array(np.rint(np.logspace(0,3,20))-1,dtype=int)
		for i in j:
			ax1.plot(ell_arr,fidu_mat[i],markers[n],linewidth=1.5, label=labels[i])
			ax2.plot(ell_arr,fidu_mat[i]/fidu_mat[contr]-1.0,markers[n],linewidth=1.5, label=labels[i])
			n+=1
	
		for i in range(len(s_arr2)):#range(5):
		
			ax1.plot(ell_arr[logind],fidu_bias_mat[i][logind],markers[n],linewidth=2, label='Bias $s$ = %s'%(s_arr2[i]),markersize=markersizes[n])
			ax2.plot(ell_arr[logind],(fidu_bias_mat[i]/fidu_mat[contr]-1.0)[logind],markers[n],linewidth=2,markersize=markersizes[n],label='Bias $s$ = %s'%(s_arr2[i]))
			n+=1
		
		ebp=logind
		ax1.errorbar(ell_arr[ebp],fidu_mat[contr,ebp],yerr[ebp],fmt='k-',linewidth=1.5)
		
		#ax.set_yscale('log',basey=5)
		plt.setp(ax1.get_xticklabels(), visible=False)
		plt.subplots_adjust(hspace=0.0)
		#leg=ax1.legend(loc=0, ncol=1, labelspacing=.2, prop={'size':16})
		#leg.get_frame().set_visible(False)
		
		xlabel(r'$\ell$',fontsize=16)
		ax1.set_ylabel(r'$\ell(\ell+1)P(\ell)/2\pi$',fontsize=16)
		ax2.set_ylabel(r'${\Delta}P/P$',fontsize=16)
		
		ax1.set_yscale('log')
		ax1.set_xscale('log')
		ax2.set_xscale('log')
		ax2.set_xlim(100,1e5)
		ax2.set_ylim(-0.25,0.4)
		ax2.set_yticks(np.arange(-0.2,0.3,0.2))
		
		handles, newlabels = ax2.get_legend_handles_labels()
		leg=ax1.legend(handles, newlabels, loc=0, ncol=1, labelspacing=.2, prop={'size':16})
		leg.get_frame().set_visible(False)
		#show()
				
		savefig(plot_dir2+'powspec_noise_'+label_str[::-1][fname]+'.pdf')
		close()
		
		fname+=1

if bias_changes_cosmo_powspec_shear:

	kmin=1.0
	kmax=2048./2
	binsize=(kmax-kmin)/1000.0
	ell_arr = (np.arange(kmin,kmax,binsize)*360./sqrt(12.0))[:1000]
	bins=1000
	
	########## grab relevant peak distributions #####
	#powspec_mat=batch_kmap_powspec(kmap_dir,contr,ps_mat_dir=ps_mat_dir,bias=0, ngal=ngal)
	#yerr=np.std(powspec_mat,axis=0)
	fidu_truek_mat=np.genfromtxt('/Users/jia/Documents/weaklensing/magbias/ps_mat/fidu_true_mat.ls')[contr]
	fidu_true_mat=np.genfromtxt('/Users/jia/Documents/weaklensing/magbias/ps_mat_shear/fidu_true_mat.ls')
			
	fidu_bias_mat = np.genfromtxt('/Users/jia/Documents/weaklensing/magbias/ps_mat_shear/fidu_bias_mat_s0.8.ls')

	
	########## actual plots #######
	f = figure(figsize=(8,6))
	
	gs = gridspec.GridSpec(2,1,height_ratios=[3,1])#width_ratios=[1,1]
	ax1=f.add_subplot(gs[0])
	ax2=f.add_subplot(gs[1],sharex=ax1)
		
	ax1.plot(ell_arr,fidu_bias_mat,linewidth=2, label='Bias $s$ = 0.8')		
	ax1.plot(ell_arr,fidu_truek_mat,linewidth=2, label='True map (kappa)')
	ax1.plot(ell_arr,fidu_true_mat,linewidth=2, label='True map (shear)')
	ax1.set_title('Comparison to Schmidt+2009 (q=1, s.t. q+1=5s-2)')
	
	ax2.plot(ell_arr,fidu_bias_mat/fidu_true_mat-1.0,linewidth=1,label='Bias $s$ = 0.8')
	ax2.plot(ell_arr,fidu_truek_mat/fidu_true_mat-1.0,linewidth=1,label='True (kappa)')
	ax2.plot(ell_arr,zeros(len(ell_arr)),'k-')
	  
	plt.setp(ax1.get_xticklabels(), visible=False)
	plt.subplots_adjust(hspace=0.0)

	xlabel(r'$\ell$',fontsize=16)
	ax1.set_ylabel(r'$\ell(\ell+1)P(\ell)/2\pi$',fontsize=16)
	ax2.set_ylabel(r'${\Delta}P/P$',fontsize=16)
	
	ax1.set_yscale('log')
	ax1.set_xscale('log')
	#ax1.set_ylim(0.9e-6,5e-4)
	
	ax2.set_xscale('log')
#	ax2.set_yscale('log')
	ax2.set_xlim(100,1e5)
	ax2.set_ylim(-0.05,0.4)#(-0.25,0.25)
	
	leg1=ax1.legend(loc=0, ncol=1, labelspacing=.2, prop={'size':10})
	leg1.get_frame().set_visible(False)
	#show()	
	savefig(plot_dir2+'powspec_shear2.jpg')
	close()
		
	
	
if bias_changes_cosmo_powspec_noise_4pt:
	ngal=15 ###
	kmin=1.0
	kmax=2048./2
	binsize=(kmax-kmin)/1000.0
	ell_arr = (np.arange(kmin,kmax,binsize)*360./sqrt(12.0))[:1000]
	bins=1000
	
	########## grab relevant peak distributions #####
	#powspec_mat=batch_kmap_powspec(kmap_dir,contr,ps_mat_dir=ps_mat_dir,bias=0, ngal=ngal)
	#yerr=np.std(powspec_mat,axis=0)
	fidu_true_mat=np.genfromtxt('/Users/jia/Documents/weaklensing/magbias/ps_mat/fidu_true_mat.ls')[contr]
			
	fidu_bias_mat = np.genfromtxt('/Users/jia/Documents/weaklensing/magbias/ps_mat/fidu_bias_mat_08.ls')[-1]
			
	## 3pt and 4pt function
	fidu_3pt_mat=np.genfromtxt('/Users/jia/magbias/fidu_bispec_s0.8_ngal15.ls')
	fidu_3pt_mat=2*(fidu_3pt_mat-fidu_true_mat)
	#fidu_bias_mat_cos=np.genfromtxt('/Users/jia/magbias/fidu_bispec_cos_s0.8_ngal15.ls')
	fidu_4pt_mat=np.genfromtxt('/Users/jia/magbias/fidu_4pt_cos_s0.8_ngal15.ls')
	fidu_sum_mat=fidu_true_mat+fidu_3pt_mat+fidu_4pt_mat
	
	########## actual plots #######
	f = figure(figsize=(8,6))
	
	gs = gridspec.GridSpec(2,1,height_ratios=[3,1])#width_ratios=[1,1]
	ax1=f.add_subplot(gs[0])
	ax2=f.add_subplot(gs[1],sharex=ax1)
	
	
	
	ax1.plot(ell_arr,fidu_bias_mat,linewidth=2, label='Bias $s$ = 0.8')	
	ax1.plot(ell_arr,fidu_3pt_mat,linewidth=2,label='3pt')
	ax1.plot(ell_arr,fidu_4pt_mat,linewidth=2,label='4pt')	
	ax1.plot(ell_arr,fidu_sum_mat,'--',linewidth=2, label='True + 3pt + 4pt')
	ax1.plot(ell_arr,fidu_true_mat,linewidth=2, label='True map')
	
	ax2.plot(ell_arr,fidu_bias_mat/fidu_true_mat-1.0,linewidth=2,label='Bias $s$ = 0.8')
	ax2.plot(ell_arr,fidu_3pt_mat/fidu_true_mat,linewidth=2,label='3pt')
	ax2.plot(ell_arr,fidu_4pt_mat/fidu_true_mat,linewidth=2,label='4pt')
	ax2.plot(ell_arr,fidu_sum_mat/fidu_true_mat-1,'--',linewidth=2, label='True + 3pt + 4pt')
	
	#ax1.plot(ell_arr[logind],fidu_bias_mat[logind],markers[n],linewidth=2, label='Bias $s$ = 0.8',markersize=markersizes[n])
		
	#ax2.plot(ell_arr[logind],(fidu_bias_mat/fidu_mat[contr]-1.0)[logind],markers[n],linewidth=2,markersize=markersizes[n],label='Bias $s$ = 0.8')
	
	
	#ax2.plot(ell_arr[logind],(fourpoint/fidu_mat[contr]-1.0)[logind],markers[n],linewidth=2,markersize=markersizes[n],label='Bias 4pt fcn $s$ = 0.8')
		
	
	plt.setp(ax1.get_xticklabels(), visible=False)
	plt.subplots_adjust(hspace=0.0)

	xlabel(r'$\ell$',fontsize=16)
	ax1.set_ylabel(r'$\ell(\ell+1)P(\ell)/2\pi$',fontsize=16)
	ax2.set_ylabel(r'${\Delta}P/P$',fontsize=16)
	
	ax1.set_yscale('log')
	ax1.set_xscale('log')
	#ax1.set_ylim(0.9e-6,5e-4)
	
	ax2.set_xscale('log')
	ax2.set_xlim(100,1e5)
	ax2.set_ylim(-0.1,0.6)#(-0.25,0.25)
	#ax2.set_yticks(np.arange(-0.2,0.3,0.2))
	
	
	leg1=ax1.legend(loc=0, ncol=1, labelspacing=.2, prop={'size':10})
	leg1.get_frame().set_visible(False)
	leg2=ax2.legend(loc=2, ncol=1, labelspacing=.2, prop={'size':10})
	leg2.get_frame().set_visible(False)
	show()
				
	#savefig(plot_dir2+'powspec_noise_bispec_'+label_str[::-1][fname]+'.jpg')
	#close()
		
if compare_04_true:
	bias04fits=np.genfromtxt(fit_noise_dir+'bias_fit_s0.4_g6.976_ngal15_bins30.ls')
	truefits=np.genfromtxt(fit_noise_dir+'true_fit_g6.976_ngal15_bins30.ls')
	(bias04fits-truefits).nonzero()
	#ngal15 failed
	bias04peak_mat=np.genfromtxt(peaks_mat_noise_dir+'peak_matrix_bias_1000m-512b240_Om0.260_Ol0.740_w-1.000_ns0.960_si0.798_s0.4_g6.976_ngal15_bins30.ls')
	truepeak_mat=np.genfromtxt(peaks_mat_noise_dir+'peak_matrix_true_1000m-512b240_Om0.260_Ol0.740_w-1.000_ns0.960_si0.798_g6.976_ngal15_bins30.ls')
	(bias04peak_mat-truepeak_mat).nonzero()
	
	fitsfile="/Users/jia/Documents/weaklensing/magbias/test.fit"
	hdulist = pyfits.open(fitsfile)
	kmap = np.array(hdulist[0].data)
	true_map=smooth(kmap,sigma)
	ngal=15
	z=1.0
	gamma=0.15+0.035*z
	noise_width=sqrt((gamma*gamma)/ngal)
	sigma = 6.976
	seed(73001)
	noise = np.array(np.random.normal(0.0, noise_width, kmap.shape),dtype=float32)
	kmap_noise = kmap + noise
	s=0.4
	kmap_true_noise = smooth (kmap_noise, sigma)
	kmap_bias_noise = biasnoisesmooth(kmap,noise,sigma,s)
	
	(kmap_true_noise-kmap_bias_noise).nonzero()# result = 0 
	

def upper_lower(x,y,lower=-1,upper=1, num=0):
	# given an array, x, y, find the position of xmin, xmax, where y is in range of (lower,upper)
	# num is when there are > 1 intervals that matches the criteria, which one do we take, mostly designed for powspec plot
	#k=np.where(abs(x-0.4)<1e-3)[0] #find through x=0.4
	
	#kk=np.where(abs(y)<upper)[0] #find through y=0
	#return kk[[0,-1]]
	
	kk=np.where(abs(y)<1.0)[0] #find through y=0
	#print 'kk=',kk
	if len(kk) > 1:
		kk_ind=np.nonzero(kk[1:]-kk[:-1]-1)[0]+1
		kk_ind=np.insert(kk_ind,0,0)
		kk=kk[kk_ind]
	else:
		num=0
	print 'kk,num=',kk,num
	k=kk[num]
	
	i=0
	j=0
	while y[k+i]>lower and y[k+i]<upper:
		if k+i==0 or k+i == len(x)-1:
			break
		i-=1
	while y[k+j]<upper and y[k+j]>lower:		
		if k+j==0 or k+j == len(x)-1:
			break
		j+=1
	return int(k+i+1),int(k+j-1)
	
if delw_vs_s:
	bins=200#30
	ngal=30
	s_arr80=np.arange(-0.5,1.0,0.01)
	analytical=1
	alpha=0.48
	
	def upper_lower(x,y,lower=-1,upper=1, num=0):
		# given an array, x, y, find the position of xmin, xmax, where y is in range of (lower,upper)
		# num is when there are > 1 intervals that matches the criteria, which one do we take, mostly designed for powspec plot
		#k=np.where(abs(x-0.4)<1e-3)[0] #find through x=0.4
		kk=np.where(abs(y)<upper)[0] #find through y=0
		return kk[[0,-1]]

	def gen_bias_fit_mat(ngal,bins):
		bias_fit_mat_file =fit_analytical_dir+'bias_analytical_mat_bins%s_ngal%s.ls'%(bins,ngal)
		if os.path.isfile(bias_fit_mat_file):
			bias_fit_mat=np.genfromtxt(bias_fit_mat_file)
		else:
			bias_fit_mat = np.ndarray(shape=(len(s_arr80),9))
			i=0
			for s in s_arr80:
				s = float('%.2f'%(s))	
				fn=np.genfromtxt(fit_analytical_dir+'bias_analytical_s%s_g6.976_ngal%s_bins%s.ls'%(s,ngal,bins)) # cols 0, 1, 2, 3: chsq, Om, w, si8
				
				chisq_f, iOm, iw, isi = average(fn, axis=0)	
				iOm_wid ,iw_wid, isi_wid = np.std(fn[:,1:],axis=0)
				
				iOmsi = iOm**alpha*isi
				iOmsi_wid = np.std(fn[:,1]**alpha*fn[:,-1], axis=0) 
				
				print iOm, iOm_wid, iw, iw_wid, isi, isi_wid, iOmsi, iOmsi_wid, chisq_f
				bias_fit_mat[i]=np.array([iOm,iw,isi,iOm_wid,iw_wid,isi_wid,iOmsi,iOmsi_wid, chisq_f])
				i+=1
				
			savetxt(bias_fit_mat_file,bias_fit_mat)
		return bias_fit_mat

	bias_fit_mat = np.ndarray(shape=(3,len(s_arr80),9)) # 3 ngals, 150 s, 9 params(iOm,iw,isi,iOm_wid,iw_wid,isi_wid,iOmsi,iOmsi_wid, chisq_f)	
	i=0
	
	for ngal in (15,30,45):
		bias_fit_mat[i]=gen_bias_fit_mat(ngal,bins)
		i+=1
	icenter=0.798*(0.26**alpha)
	
	mmarkers=('r--','g-','m:')
	mcolors=('r','g','m')
	hatchs=('-','x','.')
	xranges=((35,55),(35,55),(35,50))
	gs = gridspec.GridSpec(2,1,height_ratios=[3,1])
	gs2 = gridspec.GridSpec(3,1,height_ratios=[3,3,1])
	
	########### plot the ultimate plot for this work #####
	for i in (1,):#range(3): #params count
		f = figure(figsize=(8,10))
		
		ax1=f.add_subplot(gs2[0])
		ax2=f.add_subplot(gs2[1],sharex=ax1)
		ax3=f.add_subplot(gs2[2],sharex=ax1)
		#x0,x1=(0,-1)
		x=s_arr80#[x0:x1]

		for j in range(3): #ngal count
			isigma = bias_fit_mat[j,90,i+3]/sqrt(20000./12)
			y=(bias_fit_mat[j,:,i]-fidu_params[i])/isigma
			z=bias_fit_mat[j,:,-1]/bins#chisq	
			z=smoothz(x,z)(x)## z get rid of the spikes
	
			y2=(bias_fit_mat[j,:,6]-icenter)/bias_fit_mat[j,90,7]*sqrt(20000./12)#Om^alpha*si8
			ax1.plot(x,y,mmarkers[j],label='$n_{gal}$ = %s arcmin$^{-2}$'%(ngal_arr[j]),linewidth=2)
			ax2.plot(x,y2,mmarkers[j],linewidth=2)
			if j==1:
				for ilim in (3.0,2.0,1.0):
					a,b=upper_lower(x,y,lower=-ilim,upper=ilim)
					print '#std at 0.2, s at -/+ 1 sigma:',y[70],x[a],x[b],label_str[i]
					xwidth=x[b]-x[a]
					trans = transforms.blended_transform_factory(ax1.transData, ax1.transAxes)
					rect = patches.Rectangle((x[a],0), width=xwidth, height=1, transform=trans, color=mcolors[j], alpha=0.8/ilim)#1.0, fill=False,hatch=hatchs[j])
					ax1.add_patch(rect)
					
					a2,b2=upper_lower(x,y2,lower=-ilim,upper=ilim)
					print '#std at 0.2, s at -/+ 1 sigma:',y2[70],x[a2],x[b2],'omsi'
					xwidth2=x[b2]-x[a2]
					trans2 = transforms.blended_transform_factory(ax2.transData, ax2.transAxes)
					rect2 = patches.Rectangle((x[a2],0), width=xwidth2, height=1, transform=trans2, color=mcolors[j], alpha=0.8/ilim)#1.0, fill=False,hatch=hatchs[j])
					ax2.add_patch(rect2)
			ax3.plot(x,z,mmarkers[j],linewidth=2)

		leg=ax1.legend(loc=0, ncol=1, labelspacing=.2, prop={'size':14})
		leg.get_frame().set_visible(False)
		plt.setp(ax1.get_xticklabels(), visible=False)
		plt.setp(ax2.get_xticklabels(), visible=False)
		plt.subplots_adjust(hspace=0.0)
		
		ax1.set_ylim(-25,25)
		ax2.set_ylim(-25,30)
		ax3.set_ylim(0.65,1.2)
		ax1.set_yticks(np.arange(-20,25,5))
		ax2.set_yticks(np.arange(-20,30,5))
		ax3.set_yticks(np.arange(0.75,1.2,0.15))
		ax3.set_xlabel('$s$',fontsize=16)
		ax3.set_xlim(-0.5,1.0)
		ax2.set_ylabel(r'${\Delta}\left(\sigma_8\Omega_m^{%s}\right) / \sigma_{\left(\sigma_8\Omega_m^{%s}\right)}$'%(alpha,alpha),fontsize=16)
		ax1.set_ylabel('${\Delta}%s$ / $\sigma_{%s}$'%(label_latex[i],label_latex[i]),fontsize=16)
		ax3.set_ylabel('$\chi^2$ / dof',fontsize=16)
		
		#show()
		savefig(plot_dir2+'ana_del_%s_bins%s.pdf'%(label_str[i],bins))
		close()
			
	########## junk plot the ultimate plot for this work #####
	#for i in (1,):#range(3): #params count
		#f = figure(figsize=(8,6))
		
		#ax1=f.add_subplot(gs[0])
		#ax2=f.add_subplot(gs[1],sharex=ax1)
		##ax1.errorbar(s_arr80,(bias_fit_mat[:,i]-centers[i]),yerr=bias_fit_mat[:,i+3],fmt='--o')
		##x0,x1=(35,55)
		#x0,x1=(0,-1)
		#x=s_arr80[x0:x1]

		#for j in range(3): #ngal count
			#y=(bias_fit_mat[j,x0:x1,i]-centers[j,i])/(sigmas[j,i])
			#z=bias_fit_mat[j,x0:x1,-1]/bins
			#ax1.plot(x,y,mmarkers[j],label='$n_{gal}$ = %s arcmin$^{-2}$'%(ngal_arr[j]),linewidth=2)
			#if j==1:
				#a,b=upper_lower(x,y)
				#print a, b
				#print y[a], y[b]
				#xwidth=x[b]-x[a]
				#trans = transforms.blended_transform_factory(ax1.transData, ax1.transAxes)
				#rect = patches.Rectangle((x[a],0), width=xwidth, height=1, transform=trans, color=mcolors[j], alpha=0.5)#1.0, fill=False,hatch=hatchs[j])
				#ax1.add_patch(rect)
			#ax2.plot(x,z,mmarkers[j],linewidth=2)

		#leg=ax1.legend(loc=0, ncol=1, labelspacing=.2, prop={'size':16})
		#leg.get_frame().set_visible(False)
		#plt.setp(ax1.get_xticklabels(), visible=False)
		#plt.subplots_adjust(hspace=0.0)
			
		#ax2.set_xlabel('$s$',fontsize=16)
		#ax1.set_ylabel('${\Delta}%s$ / $\sigma_{%s}$'%(label_latex[i],label_latex[i]),fontsize=16)
		#ax2.set_ylabel('$\chi^2$ / dof',fontsize=16)
		
		#show()
		
		#if analytical:
			#savefig(plot_dir+'ana_del_%s_bins%s.pdf'%(label_str[i],bins))
			#savefig(plot_dir+'ana_del_%s_bins%s.jpg'%(label_str[i],bins))
			#close()
		#else:
			#savefig(plot_dir+'fit_del_%s_bins%s.pdf'%(label_str[i],bins))		
			#savefig(plot_dir+'fit_del_%s_bins%s.jpg'%(label_str[i],bins))
			#close()
		
	###### JUNK - plot combining sigma8 & omega mass to beat denegeracy ##########
	
	#icenter=0.798*(0.26**alpha)
	#f = figure(figsize=(8,6))
	#ax1=f.add_subplot(gs[0])
	#ax2=f.add_subplot(gs[1],sharex=ax1)
	#x=s_arr80
	#Tfit3=(Tfit_ngal15,Tfit_ngal30,Tfit_ngal45)
	#for j in range(3): #ngal count	
		#a,b=Tfit3[j][2,:,[1,-1]]
		#c=a**alpha*b
		#c=c[~np.isnan(c)]
		#isigma=np.std(c)/sqrt(20000./12)
		#y=(bias_fit_mat[j,:,2]*(bias_fit_mat[j,:,0]**alpha)-icenter)/isigma
		
		#z=bias_fit_mat[j,:,-1]/bins
		#ax1.plot(x,y,color=mcolors[j],label='$n_{gal}$=%s/ arcmin'%(ngal_arr[j]),linewidth=2)
		#ax2.plot(x,z,color=mcolors[j],linewidth=2)

	#leg=ax1.legend(loc=0, ncol=1, labelspacing=.2, prop={'size':16})
	#leg.get_frame().set_visible(False)
	
	#ax2.set_ylim(0.6,1.3)
	#ax2.set_xlabel('$s$')
	##ax1.set_ylabel('${\Delta}\sigma_8/\sigma_8+0.6\Delta\Omega_m/\Omega_m$')
	#ax1.set_ylabel('$\sigma_8\Omega_m^{%s}$'%(alpha))
	#ax2.set_ylabel('$\chi^2$ / dof')
	
	#show()
	#savefig(plot_dir2+'chisq_zoom_si8om06_%s.jpg'%(interpmethod))	


if del_w_vs_s_fw_bw:
	ngal = 30
	bins = 200
	#fit_analytical_dir = '/Users/jia/Documents/weaklensing/magbias/fit_noise_sim45/'
	fit_analytical_dir = '/Users/jia/Documents/weaklensing/magbias/ps_fit_noise_sim45/'
	alpha = 0.62
	
	fn = lambda interpmethod, s: np.genfromtxt(fit_analytical_dir+'%s/powspec_bias_analytical_s%s_ngal30.ls'%(interpmethod, s))
	
	fn_true = np.genfromtxt(fit_analytical_dir+'fw/powspec_true_analytical_ngal30.ls')
	
	#fn = lambda interpmethod, s: np.genfromtxt(fit_analytical_dir+'%s/bias_analytical_s%s_g6.976_ngal%s_bins%s.ls'%(interpmethod, s,ngal,bins))
	
	#fn_true = np.genfromtxt(fit_analytical_dir+'fw/true_analytical_g6.976_ngal30_bins200.ls')
	iOm_wid, iw_wid, isi_wid = np.std(fn_true[:,1:],axis=0)/sqrt(20000./12)
#	iw_wid /= sqrt(20000./12)
	iOmsi_wid = np.std(fn_true[:,1]**alpha*fn_true[:,-1], axis=0)/sqrt(20000./12)
	
	omsicentr = 0.798*0.26**alpha
	for s in (0.2, 0.8):
		for im in ('fw','bw'):
			chisq_f, iOm, iw, isi = average(fn(im,s), axis=0)	
			iOmsi = iOm**alpha*isi
			print '%s\t%.1f\t%.1f\t%.1f\t%.1f\t%.1f'%(im, s, (iOm-0.26)/iOm_wid, (iw+1.0)/iw_wid, (isi-0.798)/isi_wid, (iOmsi-omsicentr)/iOmsi_wid)
		
		

if delw_vs_s_powspec_noise:
	bins=202
	ngal=30
	s_arr80=np.arange(-0.5,1.0,0.01)
	analytical=1
	alpha=0.62
	 
	def gen_bias_fit_mat(ngal,bins):
		bias_fit_mat_file =ps_fit_dir+'mat_powspec_bias_analytical_ngal%s.ls'%(ngal)
		if os.path.isfile(bias_fit_mat_file):
			bias_fit_mat=np.genfromtxt(bias_fit_mat_file)
		else:
			bias_fit_mat = np.ndarray(shape=(len(s_arr80),9))
			i=0
			for s in s_arr80:
				s = float('%.2f'%(s))	
				fn=np.genfromtxt(ps_fit_dir+'powspec_bias_analytical_s%s_ngal%s.ls'%(s,ngal)) # cols 0, 1, 2, 3: chsq, Om, w, si8
				
				chisq_f, iOm, iw, isi = average(fn, axis=0)	
				iOm_wid ,iw_wid, isi_wid = np.std(fn[:,1:],axis=0)
				
				iOmsi = iOm**alpha*isi
				iOmsi_wid = np.std(fn[:,1]**alpha*fn[:,-1], axis=0) 
				
				print iOm, iOm_wid, iw, iw_wid, isi, isi_wid, iOmsi, iOmsi_wid, chisq_f
				bias_fit_mat[i]=np.array([iOm,iw,isi,iOm_wid,iw_wid,isi_wid,iOmsi,iOmsi_wid, chisq_f])
				i+=1
				
			savetxt(bias_fit_mat_file,bias_fit_mat)
		return bias_fit_mat

	bias_fit_mat = np.ndarray(shape=(3,len(s_arr80),9)) # 3 ngals, 150 s, 9 params(iOm,iw,isi,iOm_wid,iw_wid,isi_wid,iOmsi,iOmsi_wid, chisq_f)	
	i=0
	
	for ngal in (15,30,45):#ngal_arr:
		bias_fit_mat[i]=gen_bias_fit_mat(ngal,bins)
		i+=1
	icenter=0.798*(0.26**alpha)
	
	mmarkers=('r--','g-','m:')
	mcolors=('r','g','m')
	hatchs=('-','x','.')
	xranges=((35,55),(35,55),(35,50))
	gs = gridspec.GridSpec(2,1,height_ratios=[3,1])
	gs2 = gridspec.GridSpec(3,1,height_ratios=[3,3,1])
	
	########### plot the ultimate plot for this work #####
	for i in (1,):#range(3): #params count
		f = figure(figsize=(8,10))
		
		ax1=f.add_subplot(gs2[0])
		ax2=f.add_subplot(gs2[1],sharex=ax1)
		ax3=f.add_subplot(gs2[2],sharex=ax1)
		#x0,x1=(0,-1)
		x=s_arr80#[x0:x1]

		for j in range(3): #ngal count
			isigma = bias_fit_mat[j,90,i+3]/sqrt(20000./12)
			y=(bias_fit_mat[j,:,i]-fidu_params[i])/isigma
			z=bias_fit_mat[j,:,-1]/bins#chisq	
			z=smoothz(x,z)(x)## z get rid of the spikes
	
			y2=(bias_fit_mat[j,:,6]-icenter)/bias_fit_mat[j,90,7]*sqrt(20000./12)#Om^alpha*si8
			ax1.plot(x,y,mmarkers[j],label='$n_{gal}$ = %s arcmin$^{-2}$'%(ngal_arr[j]),linewidth=2)
			ax2.plot(x,y2,mmarkers[j],label='$n_{gal}$ = %s arcmin$^{-2}$'%(ngal_arr[j]),linewidth=2)
			if j==1:
				for ilim in (3.0,2.0,1.0):
					a,b=upper_lower(x,y,lower=-ilim,upper=ilim)
					print '#std at 0.1, s at -/+ 1 sigma:',y[60],x[a],x[b],label_str[i]
					xwidth=x[b]-x[a]
					trans = transforms.blended_transform_factory(ax1.transData, ax1.transAxes)
					rect = patches.Rectangle((x[a],0), width=xwidth, height=1, transform=trans, color=mcolors[j], alpha=0.8/ilim)#1.0, fill=False,hatch=hatchs[j])
					ax1.add_patch(rect)
					
					a2,b2=upper_lower(x,y2,lower=-ilim,upper=ilim)
					print '#std at 0.1, s at -/+ 1 sigma:',y2[60],x[a2],x[b2]
					xwidth2=x[b2]-x[a2]
					trans2 = transforms.blended_transform_factory(ax2.transData, ax2.transAxes)
					rect2 = patches.Rectangle((x[a2],0), width=xwidth2, height=1, transform=trans2, color=mcolors[j], alpha=0.8/ilim)#1.0, fill=False,hatch=hatchs[j])
					ax2.add_patch(rect2)
			ax3.plot(x,z,mmarkers[j],linewidth=2)

		leg=ax2.legend(loc=0, ncol=1, labelspacing=.2, prop={'size':14})
		leg.get_frame().set_visible(False)
		plt.setp(ax1.get_xticklabels(), visible=False)
		plt.setp(ax2.get_xticklabels(), visible=False)
		plt.subplots_adjust(hspace=0.0)
		
		ax1.set_ylim(-22,6)
		ax1.set_yticks(np.arange(-20,6,5))
		ax2.set_ylim(-105,60)
		ax2.set_yticks(np.arange(-100,60,25))
		ax3.set_ylim(0.7,1.2)
		ax3.set_yticks(np.arange(0.75,1.2,0.15))
		ax3.set_xlim(-0.5,1.0)
		
		ax3.set_xlabel('$s$',fontsize=16)
		ax3.set_xlim(-0.5,1.0)
		ax2.set_ylabel(r'${\Delta}\left(\sigma_8\Omega_m^{%s}\right) / \sigma_{\left(\sigma_8\Omega_m^{%s}\right)}$'%(alpha,alpha),fontsize=16)
		ax1.set_ylabel('${\Delta}%s$ / $\sigma_{%s}$'%(label_latex[i],label_latex[i]),fontsize=16)
		ax3.set_ylabel('$\chi^2$ / dof',fontsize=16)
		
		#show()
		savefig(plot_dir2+'powspec_noise_del_%s_%s.pdf'%(label_str[i],interpmethod))
		close()

if delw_vs_s_powspec:
	bins=203
	ngal=30
	s_arr80=np.arange(-0.5,1.0,0.01)
	alpha=0.57#0.66
	
	def gen_bias_fit_mat(ngal):
		bias_fit_mat_file =ps_fit_dir+'mat_powspec_bias_analytical_ngal%s.ls'%(ngal)
		if os.path.isfile(bias_fit_mat_file):
			bias_fit_mat=np.genfromtxt(bias_fit_mat_file)
		else:
			bias_fit_mat = np.ndarray(shape=(len(s_arr80),7))
			##7 columns: iOm,iw,isi,iOm_wid,iw_wid,isi_wid,chisq_f
			i=0
			for s in s_arr80:
				s = float('%.2f'%(s))
				f=np.genfromtxt(ps_fit_dir+'powspec_bias_analytical_s%s_ngal%s.ls'%(s,ngal))
				
				bias_fits=f[:,1:]
				chisq_f=average(f[:,0])
				Oms	=	stats.bayes_mvs(bias_fits[:,0])
				ws	=	stats.bayes_mvs(bias_fits[:,1])
				sis	=	stats.bayes_mvs(bias_fits[:,2])
				
				iOm	=	Oms[0][0]
				iOm_wid	=	0.5*(Oms[2][1][0]+Oms[2][1][1])
				iw	=	ws[0][0]
				iw_wid	=	0.5*(ws[2][1][0]+ws[2][1][1])
				isi	=	sis[0][0]
				isi_wid	=	0.5*(sis[2][1][0]+sis[2][1][1])	
				print iOm,iOm_wid,iw,iw_wid,isi,isi_wid,chisq_f
				bias_fit_mat[i]=np.array([iOm,iw,isi,iOm_wid,iw_wid,isi_wid,chisq_f])
				i+=1
				
			savetxt(bias_fit_mat_file,bias_fit_mat)
		return bias_fit_mat
	
	bias_fit_mat = np.ndarray(shape=(3,len(s_arr80),7))	
	sigmas=np.ndarray(shape=(3,4))
	
	############ get Tfit to compute error #########
	f15 = np.genfromtxt(ps_fit_dir+'powspec_true_analytical_ngal15.ls')
	f30 = np.genfromtxt(ps_fit_dir+'powspec_true_analytical_ngal30.ls')
	f45 = np.genfromtxt(ps_fit_dir+'powspec_true_analytical_ngal45.ls')
	Tfit3=(f15,f30,f45)
	###################################################
	
	i=0	
	for ngal in ngal_arr:
		bias_fit_mat[i]=gen_bias_fit_mat(ngal)	
		sigmas[i,:-1]=bias_fit_mat[i,90,3:6]/sqrt(20000./12) # 3 sigmas for 3 params, 4th is for Omsi8 combination, 90 is the 0.4 one		
		a,b=(Tfit3[i][:,[1,-1]]).T
		c=a**alpha*b
		c=c[~np.isnan(c)]
		sigmas[i,-1]=np.std(c)/sqrt(20000./12)		
		i+=1
		
	icenter=0.798*(0.26**alpha)	
	
	mmarkers=('r--','g-','m:')
	mcolors=('r','g','m')
	hatchs=('-','x','.')
	xranges=((35,55),(35,55),(35,50))
	gs = gridspec.GridSpec(2,1,height_ratios=[3,1])
	gs2 = gridspec.GridSpec(3,1,height_ratios=[3,3,1])
	

	########### plot2 the ultimate plot for this work #####
	for i in (1,):#range(3): #params count
		f = figure(figsize=(8,10))
		
		ax1=f.add_subplot(gs2[0])
		ax2=f.add_subplot(gs2[1],sharex=ax1)
		ax3=f.add_subplot(gs2[2],sharex=ax1)
		#x0,x1=(0,-1)
		x=s_arr80#[x0:x1]

		for j in range(3): #ngal count
			y=(bias_fit_mat[j,:,i]-fidu_params[i])/(sigmas[j,i])
			z=bias_fit_mat[j,:,-1]/bins			
			z=smoothz(x,z)(x)## z get rid of the spikes
			y2=(bias_fit_mat[j,:,2]*(bias_fit_mat[j,:,0]**alpha)-icenter)/sigmas[j,3]
			ax1.plot(x,y,mmarkers[j],label='$n_{gal}$ = %s arcmin$^{-2}$'%(ngal_arr[j]),linewidth=2)
			ax2.plot(x,y2,mmarkers[j],label='$n_{gal}$ = %s arcmin$^{-2}$'%(ngal_arr[j]),linewidth=2)
			if j==1:
				yyy=y
				for ilim in (3.0,2.0,1.0):
					for num in (0,):#1):
						a,b=upper_lower(x,y,lower=-ilim,upper=ilim,num=num)
						print '#std at 0.1, s at -/+ 1 sigma:',y[60],x[a],x[b],label_str[i]
						xwidth=x[b]-x[a]
						trans = transforms.blended_transform_factory(ax1.transData, ax1.transAxes)
						rect = patches.Rectangle((x[a],0), width=xwidth, height=1, transform=trans, color=mcolors[j], alpha=1.0/ilim)#1.0, fill=False,hatch=hatchs[j])
						if ilim==3.0 and num==1:
							# this is to avoid layer too many shades together
							continue
						else:
							ax1.add_patch(rect)
					
					a2,b2=upper_lower(x,y2,lower=-ilim,upper=ilim)
					print '#std at 0.1, s at -/+ 1 sigma:',y2[60],x[a2],x[b2]
					xwidth2=x[b2]-x[a2]
					trans2 = transforms.blended_transform_factory(ax2.transData, ax2.transAxes)
					rect2 = patches.Rectangle((x[a2],0), width=xwidth2, height=1, transform=trans2, color=mcolors[j], alpha=1.0/ilim)#1.0, fill=False,hatch=hatchs[j])
					ax2.add_patch(rect2)
			ax3.plot(x,z,mmarkers[j],linewidth=2)

		leg=ax2.legend(loc=0, ncol=1, labelspacing=.2, prop={'size':14})
		leg.get_frame().set_visible(False)
		plt.setp(ax1.get_xticklabels(), visible=False)
		plt.setp(ax2.get_xticklabels(), visible=False)
		plt.subplots_adjust(hspace=0.0)
		
		#ax1.set_yticks(np.arange(-40,20,10))
		#ax2.set_yticks(np.arange(-80,20,20))
		
		#ax1.set_ylim(-20,15)
		#ax2.set_ylim(-20,15)
		#ax1.set_yticks(np.arange(-15,15,5))
		#ax2.set_yticks(np.arange(-15,15,5))
		
		ax3.set_xlabel('$s$',fontsize=16)
		ax2.set_ylabel(r'${\Delta}\left(\sigma_8\Omega_m^{%s}\right) / \sigma_{\left(\sigma_8\Omega_m^{%s}\right)}$'%(alpha,alpha),fontsize=16)
		ax1.set_ylabel('${\Delta}%s$ / $\sigma_{%s}$'%(label_latex[i],label_latex[i]),fontsize=16)
		ax3.set_ylabel('$\chi^2$ / dof',fontsize=16)
		#ax3.set_ylim(0,18)
		#ax3.set_yticks(np.arange(0,18,5))
		
		show()
		
		#savefig(plot_dir2+'powspec_del_%s_%s.pdf'%(label_str[i],interpmethod))
		#close()

if delw_vs_s_powspec_noise_old:
	bins=203
	ngal=30
	s_arr80=np.arange(-0.5,1.0,0.01)
	alpha=0.62 #noiseless: 0.66
	
	def gen_bias_fit_mat(ngal):
		bias_fit_mat_file =ps_fit_dir+'mat_powspec_bias_analytical_ngal%s.ls'%(ngal)
		if os.path.isfile(bias_fit_mat_file):
			bias_fit_mat=np.genfromtxt(bias_fit_mat_file)
		else:
			bias_fit_mat = np.ndarray(shape=(len(s_arr80),7))
			##7 columns: iOm,iw,isi,iOm_wid,iw_wid,isi_wid,chisq_f
			i=0
			for s in s_arr80:
				s = float('%.2f'%(s))
				f=np.genfromtxt(ps_fit_dir+'powspec_bias_analytical_s%s_ngal%s.ls'%(s,ngal))
				
				bias_fits=f[:,1:]
				chisq_f=average(f[:,0])
				Oms	=	stats.bayes_mvs(bias_fits[:,0])
				ws	=	stats.bayes_mvs(bias_fits[:,1])
				sis	=	stats.bayes_mvs(bias_fits[:,2])
				
				iOm	=	Oms[0][0]
				iOm_wid	=	0.5*(Oms[2][1][0]+Oms[2][1][1])
				iw	=	ws[0][0]
				iw_wid	=	0.5*(ws[2][1][0]+ws[2][1][1])
				isi	=	sis[0][0]
				isi_wid	=	0.5*(sis[2][1][0]+sis[2][1][1])	
				print iOm,iOm_wid,iw,iw_wid,isi,isi_wid,chisq_f
				bias_fit_mat[i]=np.array([iOm,iw,isi,iOm_wid,iw_wid,isi_wid,chisq_f])
				i+=1
				
			savetxt(bias_fit_mat_file,bias_fit_mat)
		return bias_fit_mat
	
	bias_fit_mat = np.ndarray(shape=(3,len(s_arr80),7))	
	sigmas=np.ndarray(shape=(3,4))
	
	############ get Tfit to compute error #########
	f15 = np.genfromtxt(ps_fit_dir+'powspec_true_analytical_ngal15.ls')
	f30 = np.genfromtxt(ps_fit_dir+'powspec_true_analytical_ngal30.ls')
	f45 = np.genfromtxt(ps_fit_dir+'powspec_true_analytical_ngal45.ls')
	Tfit3=(f15,f30,f45)
	###################################################
	
	i=0	
	for ngal in ngal_arr:
		bias_fit_mat[i]=gen_bias_fit_mat(ngal)	
		sigmas[i,:-1]=bias_fit_mat[i,90,3:6]/sqrt(20000./12) # 3 sigmas for 3 params, 4th is for Omsi8 combination, 90 is the 0.4 one		
		a,b=(Tfit3[i][:,[1,-1]]).T
		c=a**alpha*b
		c=c[~np.isnan(c)]
		sigmas[i,-1]=np.std(c)/sqrt(20000./12)		
		i+=1
		
	icenter=0.798*(0.26**alpha)	
	
	mmarkers=('r--','g-','m:')
	mcolors=('r','g','m')
	hatchs=('-','x','.')
	xranges=((35,55),(35,55),(35,50))
	gs = gridspec.GridSpec(2,1,height_ratios=[3,1])
	gs2 = gridspec.GridSpec(3,1,height_ratios=[3,3,1])
	

	########### plot2 the ultimate plot for this work #####
	for i in (1,):#range(3): #params count
		f = figure(figsize=(8,10))
		
		ax1=f.add_subplot(gs2[0])
		ax2=f.add_subplot(gs2[1],sharex=ax1)
		ax3=f.add_subplot(gs2[2],sharex=ax1)
		#x0,x1=(0,-1)
		x=s_arr80#[x0:x1]

		for j in range(3): #ngal count
			y=(bias_fit_mat[j,:,i]-fidu_params[i])/(sigmas[j,i])
			z=bias_fit_mat[j,:,-1]/bins			
			z=smoothz(x,z)(x)## z get rid of the spikes
			y2=(bias_fit_mat[j,:,2]*(bias_fit_mat[j,:,0]**alpha)-icenter)/sigmas[j,3]
			ax1.plot(x,y,mmarkers[j],label='$n_{gal}$ = %s arcmin$^{-2}$'%(ngal_arr[j]),linewidth=2)
			ax2.plot(x,y2,mmarkers[j],label='$n_{gal}$ = %s arcmin$^{-2}$'%(ngal_arr[j]),linewidth=2)
			if j==1:
				yyy=y
				for ilim in (3.0,2.0,1.0):
					for num in (0,):#1):
						a,b=upper_lower(x,y,lower=-ilim,upper=ilim,num=num)
						if ilim==1.0 and num ==0:
							ind01=np.where(x==0.1+np.min(abs(x-0.1)))[0][0]
							ind08=np.where(x==0.8+np.min(abs(x-0.8)))[0][0]
							print '#std at 0.1, 0.8, s at -/+ 1 sigma:',y[ind01], y[ind08], x[a],x[b],label_str[i]
							print '#std at 0.1, 0.8, s at -/+ 1 sigma:',y2[ind01], y2[ind08], x[a2],x[b2],'omega^alpha*sigma8'
						xwidth=x[b]-x[a]
						trans = transforms.blended_transform_factory(ax1.transData, ax1.transAxes)
						rect = patches.Rectangle((x[a],0), width=xwidth, height=1, transform=trans, color=mcolors[j], alpha=1.0/ilim)#1.0, fill=False,hatch=hatchs[j])
						
						
						if ilim!=1.0 and num==1:
							# this is to avoid layer too many shades together
							continue
						else:
							ax1.add_patch(rect)
					
					a2,b2=upper_lower(x,y2,lower=-ilim,upper=ilim)
					xwidth2=x[b2]-x[a2]
					trans2 = transforms.blended_transform_factory(ax2.transData, ax2.transAxes)
					rect2 = patches.Rectangle((x[a2],0), width=xwidth2, height=1, transform=trans2, color=mcolors[j], alpha=1.0/ilim)#1.0, fill=False,hatch=hatchs[j])
					ax2.add_patch(rect2)
			ax3.plot(x,z,mmarkers[j],linewidth=2)

		leg=ax2.legend(loc=0, ncol=1, labelspacing=.2, prop={'size':14})
		leg.get_frame().set_visible(False)
		plt.setp(ax1.get_xticklabels(), visible=False)
		plt.setp(ax2.get_xticklabels(), visible=False)
		plt.subplots_adjust(hspace=0.0)
		
		
		ax1.set_ylim(-22,6)
		ax1.set_yticks(np.arange(-20,6,5))
		ax2.set_ylim(-105,60)
		ax2.set_yticks(np.arange(-100,60,25))
		
		ax3.set_xlabel('$s$',fontsize=16)
		ax2.set_ylabel(r'${\Delta}\left(\sigma_8\Omega_m^{%s}\right) / \sigma_{\left(\sigma_8\Omega_m^{%s}\right)}$'%(alpha,alpha),fontsize=16)
		ax1.set_ylabel('${\Delta}%s$ / $\sigma_{%s}$'%(label_latex[i],label_latex[i]),fontsize=16)
		ax3.set_ylabel('$\chi^2$ / dof',fontsize=16)
		
		ax3.set_ylim(0.7,1.2)
		ax3.set_yticks(np.arange(0.7,1.2,0.15))
		ax3.set_xlim(-0.5,1.0)
		show()
		#savefig(plot_dir2+'powspec_noise_del_%s_%s.pdf'%(label_str[i],interpmethod))
		#close()
		
if powspec_test:


#####################################
	#import radialProfile_jia
	from scipy import fftpack
	import healpy as hp

	fitsfile="/Users/jia/Dropbox/test_z1.fit"
	#"/Users/jia/Documents/weaklensing/magbias/test.fit"
	hdulist = pyfits.open(fitsfile)
	kmap = np.array(hdulist[0].data)
	kmap=smooth(kmap,0.0)
	true_map=smooth(kmap,sigma)
	
	#F = fftpack.fftshift(fftpack.fft2(true_map))
	F = fftpack.fftshift(fftpack.fft2(kmap))
	psd2D = np.abs(F)**2
	ell_arr,psd1D, psd1D2 = azimuthalAverage(psd2D)
	x=ell_arr
	
	noise_width=0.33
	noise = np.array(np.random.normal(0.0, noise_width, kmap.shape),dtype=float32)
	F_noise=fftpack.fftshift(fftpack.fft2(kmap+noise))
	psd2D_noise = np.abs(F_noise)**2
	ell_arr_noise,psd1D_noise, psd1D2_noise = azimuthalAverage(psd2D_noise)
	psjan=genfromtxt('/Users/jia/magbias/powspec/power_spectrum_test_z1.txt').T
	
	f=figure(figsize=(8,6))
	ax=f.add_subplot(111)
	
	######normalization
	norm = ((2*pi*sqrt(12.00)/360.0)**2)/(2048.0**2)**2
	
	ax.plot(x,x*(x+1)*psd1D/2/pi*norm,label="jia (no noise)",linewidth=1)
	ax.plot(x,x*(x+1)*psd1D2/2/pi*norm,label="jia (l(l+1) before binning, no noise)",linewidth=1)
	
	ax.plot(x,x*(x+1)*psd1D_noise/2/pi*norm,label="jia (noise)",linewidth=1)
	
	x2=psjan[0]
	ax.plot(psjan[0],x2*(x2+1)*psjan[1]/2/pi,label="jan",linewidth=1)
	legend()
	
	xlim(100,1.5e5)
	ax.set_yscale('log')
	ax.set_xscale('log')
	ax.set_xlabel('l')
	ax.set_ylabel(r'l(l+1)P(l)/2$\pi$')
	show()
	

if sigmak:
	sigma = 6.976 
	fitsfile=lambda n: "/Users/jia/weaklensing/map2/m-512b240_Om0.260_Ol0.740_w-1.000_ns0.960_si0.798/WL-conv_m-512b240_Om0.260_Ol0.740_w-1.000_ns0.960_si0.798_4096xy_00%sr_0029p_0100z_og.gre.fit"%(n)
	def sigmak(n):
		
		hdulist = pyfits.open(fitsfile(n))
		kmap = np.array(hdulist[0].data)
		return np.std(kmap)
	x = range(10,99)
	m = 0
	for n in x:
		sk = sigmak(n)
		print n, sk
		m += sk
	print m/len(x)
	
if ells_fidu_s01:
	
	alpha=0.44
	os0=0.26**alpha*0.798
	ngal=30
	bins=30
	fits_fidu=np.genfromtxt(fit_analytical_dir+'bias_analytical_s%s_g6.976_ngal%s_bins%s.ls'%(0.4,ngal,bins))
	fits_s01=np.genfromtxt(fit_analytical_dir+'bias_analytical_s%s_g6.976_ngal%s_bins%s.ls'%(0.1,ngal,bins))
	fits_s08=np.genfromtxt(fit_analytical_dir+'bias_analytical_s%s_g6.976_ngal%s_bins%s.ls'%(0.8,ngal,bins))
	fits_pow_fidu=np.genfromtxt(ps_fit_dir+'powspec_bias_analytical_s0.4_ngal30.ls')
	fits_pow_s01=np.genfromtxt(ps_fit_dir+'powspec_bias_analytical_s0.1_ngal30.ls')
	fits_pow_s08=np.genfromtxt(ps_fit_dir+'powspec_bias_analytical_s0.8_ngal30.ls')
	
	
	w1=fits_fidu.T[2]
	w2=fits_s01.T[2]
	
	w1_pow=fits_pow_fidu.T[2]
	w2_pow=fits_pow_s01.T[2]
	
	os20=average(fits_s01.T[1])**alpha*average(fits_s01.T[3])
	os20_pow=average(fits_pow_s01.T[1])**alpha*average(fits_pow_s01.T[3])
	
	os1=fits_fidu.T[1]**alpha*fits_fidu.T[3]
	os2=fits_s01.T[1]**alpha*fits_s01.T[3]	
	os1_pow=fits_pow_fidu.T[1]**alpha*fits_pow_fidu.T[3]
	os2_pow=fits_pow_s01.T[1]**alpha*fits_pow_s01.T[3]	
	#os1=os1[~np.isnan(os1)]
	#os2=os2[~np.isnan(os2)]
	
	w1b=fidu_params[1]+(w1-fidu_params[1])/sqrt(20000./12)
	w2b=fidu_params[1]+(w2-fidu_params[1])/sqrt(20000./12)
	w1b_pow=fidu_params[1]+(w1_pow-fidu_params[1])/sqrt(20000./12)
	w2b_pow=fidu_params[1]+(w2_pow-fidu_params[1])/sqrt(20000./12)
	
	os1b=os0+(os1-average(os1[~np.isnan(os1)]))/sqrt(20000./12)
	os2b=os0+(os2-average(os2[~np.isnan(os2)]))/sqrt(20000./12)
	os1b_pow=os0+(os1_pow-average(os1_pow[~np.isnan(os1)]))/sqrt(20000./12)
	os2b_pow=os0+(os2_pow-average(os2_pow[~np.isnan(os2)]))/sqrt(20000./12)
	
	os1c=os1[~np.isnan(os1)]
	w1c=w1[~np.isnan(os1)]
	os1c_pow=os1_pow[~np.isnan(os1_pow)]
	w1c_pow=w1_pow[~np.isnan(os1_pow)]
	
	os2c=os2[~np.isnan(os2)]
	w2c=w2[~np.isnan(os2)]
	os2c_pow=os2_pow[~np.isnan(os2_pow)]
	w2c_pow=w2_pow[~np.isnan(os2_pow)]
	
	#stdw1=np.std(w1c)/sqrt(20000./12)
	#stdw2=np.std(w2c)/sqrt(20000./12)
	#stdos1=np.std(os1c)/sqrt(20000./12)
	#stdos2=np.std(os2c)/sqrt(20000./12)
	
	
	def plotEllipse(pos,P,edge,ls,ilabel):
		LSSTr=sqrt(20000./12)
		U, s, Vh = svd(P)
		orient = math.atan2(U[1,0],U[0,0])*180/pi
		ellipsePlot = Ellipse(xy=pos, width=2.0*math.sqrt(s[0])/LSSTr, height=2.0*math.sqrt(s[1])/LSSTr, angle=orient,edgecolor=edge, fill = False,label=ilabel,ls=ls,linewidth=1.5)#facecolor=edge)
		ax = gca()
		ax.add_patch(ellipsePlot)
		return ellipsePlot
	
	f=figure(figsize=(8,6))
	ax = f.add_subplot(111) 
	
	pos1=(-1,os0)
	pos2=(average(w2),os20)
	P1=np.cov(w1c,os1c)
	P2=np.cov(w2c,os2c)
	
	pos1_pow=(-1,os0)
	pos2_pow=(average(w2_pow),os20_pow)
	P1_pow=np.cov(w1c_pow,os1c_pow)
	P2_pow=np.cov(w2c_pow,os2c_pow)
	
	plotEllipse(pos1,P1,'r','solid','Fiducial (peak)')
	plotEllipse(pos2,P2,'b','dashdot','$s=0.1$ (peak)')
	
	plotEllipse(pos1_pow,P1_pow,'m','dotted','Fiducial (pow spec)')
	plotEllipse(pos2_pow,P2_pow,'g','dashed','$s=0.1$ (pow spec)')
	
	ax.plot(pos1[0],pos1[1],'ko')
	ax.plot(pos2[0],pos2[1],'ko')	
	ax.plot(pos2_pow[0],pos2_pow[1],'ko')
	
	#legend()
	ax.set_xlim(-1.1,-0.95)
	#ax.set_ylim(0.453,0.47)
	ax.set_ylim(0.423,0.445)
	
	ax.set_xlabel('w',fontsize=16)
	ax.set_ylabel(r'$\sigma_8\Omega_m^{%s}$'%(alpha),fontsize=16)
	
	ax.text(-1.0,0.443,'fiducial(peak)',color='r',fontsize=14)
	ax.text(-1.02,0.439,'fiducial(pow spec)',color='m',fontsize=14)
	
	#fw
	ax.text(-1.08,0.4353,'$s=0.1$ (peak)',color='b',fontsize=14)	
	ax.text(-1.0,0.428,'$s=0.1$ (pow spec)',color='g',fontsize=14)
	
	## for alpha=0.4
	#ax.text(-1.0,0.467,'fiducial(peak)',color='r',fontsize=14)
	#ax.text(-1.02,0.4635,'fiducial(pow spec)',color='m',fontsize=14)
	
	##fw
	#ax.text(-1.08,0.46,'$s=0.1$ (peak)',color='b',fontsize=14)	
	#ax.text(-1.0,0.456,'$s=0.1$ (pow spec)',color='g',fontsize=14)
	
	#bw
	#ax.text(-1.16,0.46,'$s=0.1$ (peak)',color='b',fontsize=14)	
	#ax.text(-1.04,0.453,'$s=0.1$ (pow spec)',color='g',fontsize=14)
	#show()
	savefig(plot_dir2+'w_os_ells_peak_powspec_%s.pdf'%(interpmethod))
	close()


if ells_fidu_s0108:
	
	alpha=0.48
	os0=0.26**alpha*0.798
	ngal=30
	bins=200
	fits_fidu=np.genfromtxt(fit_analytical_dir+'bias_analytical_s%s_g6.976_ngal%s_bins%s.ls'%(0.4,ngal,bins))
	fits_s01=np.genfromtxt(fit_analytical_dir+'bias_analytical_s%s_g6.976_ngal%s_bins%s.ls'%(0.2,ngal,bins))
	fits_s08=np.genfromtxt(fit_analytical_dir+'bias_analytical_s%s_g6.976_ngal%s_bins%s.ls'%(0.8,ngal,bins))
	fits_pow_fidu=np.genfromtxt(ps_fit_dir+'powspec_bias_analytical_s0.4_ngal30.ls')
	fits_pow_s01=np.genfromtxt(ps_fit_dir+'powspec_bias_analytical_s0.2_ngal30.ls')
	fits_pow_s08=np.genfromtxt(ps_fit_dir+'powspec_bias_analytical_s0.8_ngal30.ls')
	
	
	def prepareplot(fits_fidu):

		w_arr=fits_fidu.T[2]
		w_origin=average(w_arr)
		os_origin=average(fits_fidu.T[1])**alpha*average(fits_fidu.T[3])
		os_arr=fits_fidu.T[1]**alpha*fits_fidu.T[3]
		#w_os_mat=np.array([w_arr[~np.isnan(os_arr)],os_arr[~np.isnan(os_arr)]])
		P=np.cov(w_arr[~np.isnan(os_arr)],os_arr[~np.isnan(os_arr)])
		
		return w_origin, os_origin, P
		
		
	def plotEllipse(pos,P,edge='r',ls='solid'):#,ilabel=None):
		LSSTr=sqrt(20000./12)
		U, s, Vh = svd(P)
		orient = math.atan2(U[1,0],U[0,0])*180/pi
		ellipsePlot = Ellipse(xy=pos, width=2.0*math.sqrt(s[0])/LSSTr, height=2.0*math.sqrt(s[1])/LSSTr, angle=orient, fill = False,linewidth=1.5,edgecolor=edge,ls=ls)#label=ilabel,ls=ls,facecolor=edge)
		ax = gca()
		ax.add_patch(ellipsePlot)
		return ellipsePlot
	
	f=figure(figsize=(8,6))
	ax = f.add_subplot(111) 

	edgecolors=('r','b','m','g','c','y')
	lss=['solid']*3+['dotted']*3
	#locs=np.array([[-1.04,0.443],[-1.1,0.435],[-0.9755,0.447],[-1.02,0.438],[-1.0,0.426],[-1.13,0.447]])
	locs=np.array([[-1.0,0.416],[-1.08,0.416],[-0.9755,0.421],[-1.05,0.42],[-1.02,0.407],[-1.11,0.425]])
	
	ell_labels=('fiducial (peak)','$s=0.2$ (peak)','$s=0.8$ (peak)','fiducial (pow spec)','$s=0.2$ (pow spec)','$s=0.8$ (pow spec)')
	
	i=0
	for fits_fidu in (fits_fidu,fits_s01,fits_s08,fits_pow_fidu,fits_pow_s01,fits_pow_s08):
		w_origin, os_origin, P = prepareplot(fits_fidu)
		ax.plot(w_origin,os_origin,'ko', ms=3)		
		plotEllipse((w_origin,os_origin),P,edgecolors[i],lss[i])
		ax.text(locs[i,0],locs[i,1],ell_labels[i],color=edgecolors[i],fontsize=14)
		
		i+=1
	
	ax.set_xlabel('w',fontsize=16)
	ax.set_ylabel(r'$\sigma_8\Omega_m^{%s}$'%(alpha),fontsize=16)
	
	ax.set_xlim(-1.15,-0.9)
	ax.set_ylim(0.4,0.435)
	
	#show()
	savefig(plot_dir2+'w_os_ells_peak_powspec_0108.pdf')
	close()
	
if ells_fidu_s01b:
	
	alpha1=0.44
	alpha2=0.66
	os0_1=0.26**alpha1*0.798
	os0_2=0.26**alpha2*0.798
	ngal=30
	bins=30
	fits_fidu=np.genfromtxt(fit_analytical_dir+'bias_analytical_s%s_g6.976_ngal%s_bins%s.ls'%(0.4,ngal,bins))
	fits_s01=np.genfromtxt(fit_analytical_dir+'bias_analytical_s%s_g6.976_ngal%s_bins%s.ls'%(0.1,ngal,bins))
	
	fits_pow_fidu=np.genfromtxt(ps_fit_dir+'powspec_bias_analytical_s0.4_ngal30.ls')
	fits_pow_s01=np.genfromtxt(ps_fit_dir+'powspec_bias_analytical_s0.1_ngal30.ls')
	
	
	w1=fits_fidu.T[2]
	w2=fits_s01.T[2]
	w1_pow=fits_pow_fidu.T[2]
	w2_pow=fits_pow_s01.T[2]
	
	os20=average(fits_s01.T[1])**alpha1*average(fits_s01.T[3])
	os20_pow=average(fits_pow_s01.T[1])**alpha2*average(fits_pow_s01.T[3])
	
	os1=fits_fidu.T[1]**alpha1*fits_fidu.T[3]
	os2=fits_s01.T[1]**alpha1*fits_s01.T[3]	
	os1_pow=fits_pow_fidu.T[1]**alpha2*fits_pow_fidu.T[3]
	os2_pow=fits_pow_s01.T[1]**alpha2*fits_pow_s01.T[3]	
	#os1=os1[~np.isnan(os1)]
	#os2=os2[~np.isnan(os2)]
	
	w1b=fidu_params[1]+(w1-fidu_params[1])/sqrt(20000./12)
	w2b=fidu_params[1]+(w2-fidu_params[1])/sqrt(20000./12)
	w1b_pow=fidu_params[1]+(w1_pow-fidu_params[1])/sqrt(20000./12)
	w2b_pow=fidu_params[1]+(w2_pow-fidu_params[1])/sqrt(20000./12)
	
	os1b=os0_1+(os1-average(os1[~np.isnan(os1)]))/sqrt(20000./12)
	os2b=os0_1+(os2-average(os2[~np.isnan(os2)]))/sqrt(20000./12)
	os1b_pow=os0_2+(os1_pow-average(os1_pow[~np.isnan(os1)]))/sqrt(20000./12)
	os2b_pow=os0_2+(os2_pow-average(os2_pow[~np.isnan(os2)]))/sqrt(20000./12)
	
	os1c=os1[~np.isnan(os1)]
	w1c=w1[~np.isnan(os1)]
	os1c_pow=os1_pow[~np.isnan(os1_pow)]
	w1c_pow=w1_pow[~np.isnan(os1_pow)]
	
	os2c=os2[~np.isnan(os2)]
	w2c=w2[~np.isnan(os2)]
	os2c_pow=os2_pow[~np.isnan(os2_pow)]
	w2c_pow=w2_pow[~np.isnan(os2_pow)]
	
	#stdw1=np.std(w1c)/sqrt(20000./12)
	#stdw2=np.std(w2c)/sqrt(20000./12)
	#stdos1=np.std(os1c)/sqrt(20000./12)
	#stdos2=np.std(os2c)/sqrt(20000./12)
	
	
	def plotEllipse(pos,P,edge,ls,ilabel):
		LSSTr=sqrt(20000./12)
		U, s, Vh = svd(P)
		orient = math.atan2(U[1,0],U[0,0])*180/pi
		ellipsePlot = Ellipse(xy=pos, width=2.0*math.sqrt(s[0])/LSSTr, height=2.0*math.sqrt(s[1])/LSSTr, angle=orient,edgecolor=edge, fill = False,label=ilabel,ls=ls,linewidth=1.5)#facecolor=edge)
		ax = gca()
		ax.add_patch(ellipsePlot)
		return ellipsePlot
	
	f=figure(figsize=(8,6))
	ax = f.add_subplot(111) 
	
	pos1=(-1,os0_1)
	pos2=(average(w2),os20)
	P1=np.cov(w1c,os1c)
	P2=np.cov(w2c,os2c)
	
	pos1_pow=(-1,os0_2)
	pos2_pow=(average(w2_pow),os20_pow)
	P1_pow=np.cov(w1c_pow,os1c_pow)
	P2_pow=np.cov(w2c_pow,os2c_pow)
	
	plotEllipse(pos1,P1,'r','solid','Fiducial (peak)')
	plotEllipse(pos2,P2,'b','dashdot','$s=0.1$ (peak)')
	
	plotEllipse(pos1_pow,P1_pow,'m','dotted','Fiducial (pow spec)')
	plotEllipse(pos2_pow,P2_pow,'g','dashed','$s=0.1$ (pow spec)')
	
	ax.plot(pos1[0],pos1[1],'ko')
	ax.plot(pos2[0],pos2[1],'ko')	
	ax.plot(pos2_pow[0],pos2_pow[1],'ko')
	
	#legend()
	ax.set_xlim(-1.1,-0.95)
	ax.set_ylim(0.453,0.47)
	ax.set_xlabel('w',fontsize=16)
	ax.set_ylabel(r'$\sigma_8\Omega_m^{%s}$'%(r'$\alpha$'),fontsize=16)
	ax.text(-1.0,0.467,'fiducial(peak)',color='r',fontsize=14)
	ax.text(-1.02,0.4635,'fiducial(pow spec)',color='m',fontsize=14)
	
	#forwards text location
	ax.text(-1.08,0.46,'$s=0.1$ (peak)',color='b',fontsize=14)	
	ax.text(-1.0,0.456,'$s=0.1$ (pow spec)',color='g',fontsize=14)
	
	#backwards text location
	#ax.text(-1.16,0.46,'$s=0.1$ (peak)',color='b',fontsize=14)	
	#ax.text(-1.04,0.453,'$s=0.1$ (pow spec)',color='g',fontsize=14)
	show()
	
	#savefig(plot_dir2+'w_os_ells_peak_powspec_%s.pdf'%(interpmethod))
	#close()
	
if find_alpha:
	bins,ngal=200,30
	
	#find alpha for power spectrum
	Oms,sis=np.genfromtxt(ps_fit_dir+'powspec_true_analytical_ngal30.ls')[:,[1,3]].T
	#find alpha for peak counts
	#Oms,sis=np.genfromtxt(fit_analytical_dir+'true_analytical_g6.976_ngal30_bins%s.ls'%(bins))[:,[1,3]].T
	
	### threw out 2nd term
	alpha_residue=lambda alpha: sum(abs(alpha*(Oms-0.26)/0.26+(sis-0.798)/0.798))
	
	### keep 2nd term
	#alpha_residue=lambda alpha: sum(abs(alpha*(Oms-0.26)/0.26+(sis-0.798)/0.798+alpha*(Oms-0.26)/0.26*(sis-0.798)/0.798))
	
	# just do brute one
	#def alpha_residue(alpha):
		#a=abs(Oms**alpha*sis-0.26**alpha*0.798)
		#return sum(a[~np.isnan(a)])
	
	alpha_arr=linspace(0.2,1.0,100)
	residue_arr=np.ndarray(shape=alpha_arr.shape)
	for i in range(len(alpha_arr)):
		residue_arr[i]=alpha_residue(alpha_arr[i])
	alpha_minimize=optimize.minimize(alpha_residue,0.4,method='L-BFGS-B')
	alpha_min,alpha_fun=alpha_minimize.x,alpha_minimize.fun
	
	f=figure(figsize=(6,8))
	subplot(211)
	plot(alpha_arr,residue_arr)
	plot(alpha_min,alpha_fun,'go',label=r'$\alpha=%.2f$'%(alpha_min))
	xlabel(r'$\alpha$')
	ylabel(r'$\alpha \frac{\Delta\Omega_m}{\Omega_m}+\frac{\Delta\sigma}{\sigma}$')
	title ( r'best fit $\alpha=%.2f$ %s'%(alpha_min,interpmethod))
	
	subplot(212)
	const_comb=0.798*0.26**alpha_min
	scatter(Oms,sis)
	xlabel('$\Omega_m$')
	ylabel('$\sigma_8$')
	Oms_arr=linspace(np.amin(Oms),np.amax(Oms),100)
	sis_arr=const_comb/Oms_arr**alpha_min
	plot(Oms_arr,sis_arr,'g-',linewidth=2)
	show()
	#savefig(plot_dir2+'find_alpha_%s.jpg'%(interpmethod))
	#close()

if peaks_err_vs_bins:
	alpha=0.48
	#fit_noise_dir = '/Users/jia/Documents/weaklensing/magbias/peaks_fit_bins/'
	#fit_noise_dir = '/Users/jia/Documents/weaklensing/magbias/fit_noise_lowk/'
	#fit_noise_dir = '/Users/jia/Documents/weaklensing/magbias/fit_noise_sim45/'
	fit_noise_dir = '/Users/jia/Documents/weaklensing/magbias/fit_noise_sim45/fw/'
	
	#s_arr=np.linspace(0,1,11)
	#bins_arr1=np.arange(10,310,10)
	#bins_arr2=np.arange(350,1000,50)
	#bins_arr=np.concatenate((bins_arr1,bins_arr2))
	s_arr=(0.0, 0.2, 0.6, 0.8)
	#bins_arr1=np.arange(10,210,10)
	#bins_arr2=np.arange(200,900,50)
	#bins_arr=np.concatenate((bins_arr1,bins_arr2))
	bins_arr = np.arange(10,500,20)
	
	Bfit_ngal30 = np.ndarray(shape=(len(s_arr),len(bins_arr),1000,4)) #s,bins,#maps,4cols
	
	i=0
	for s in s_arr:
		j=0
		for bins in bins_arr:
			
			Bfit_ngal30[i,j]= np.genfromtxt(fit_noise_dir+'bias_analytical_s%.1f_g6.976_ngal30_bins%s.ls'%(s,bins))

			j += 1
		i += 1

	Tfit_ngal30 = np.ndarray(shape=(len(bins_arr),1000,4))
	i = 0

	for bins in bins_arr:

		f30 = np.genfromtxt(fit_noise_dir+'true_analytical_g6.976_ngal30_bins%s.ls'%(bins))

		Tfit_ngal30[i]=f30
		i += 1
	
	x = bins_arr
	y = average(Bfit_ngal30[:,:,:,1:],axis=2).T
	markers=('g--','y-','r:','c-.','m.')
	#markers=('go--','md-','r*:','c*-','k-','y2-','yd','2','3','1','g-.')

	######## begin - plot derived MB params vs number of bins ###########
	label_latex = ('\sigma_8\Omega_m^{%s}'%(alpha),'w')
	fidu_params = (0.798*0.26**alpha,-1.0)
	yerr = np.std(Tfit_ngal30[:,:,1:],axis=1).T
	yerr /= sqrt(20000./12)
	#ylims=(0.006, 0.2)
	
	for i in (1,0):#range(3):
		f = figure(figsize=(8,6))
		gs = gridspec.GridSpec(2,1,height_ratios=[3,1])
		ax=f.add_subplot(gs[0])
		ax2=f.add_subplot(gs[1],sharex=ax)
		iyerr=yerr[i]
		
		for j in range(len(s_arr)):
			if i == 1:
				iy=y[i,:,j]
			else:
				iy=y[0,:,j]**alpha*y[2,:,j]
			ax.plot(x, iy, markers[j], label='s = %s'%(s_arr[j]),linewidth=2)
			
		ax2.plot(x, iyerr, 'm-',linewidth=2)
		ax2.set_xlabel('Number of bins',fontsize=16)
		
		ax2.set_ylabel(r'$\sigma_{%s}$'%(label_latex[i]),fontsize=16)
		ax.set_ylabel(r'$%s$'%(label_latex[i]),fontsize=16)
		
		
		leg=ax.legend(loc=0, ncol=1, labelspacing=.2, prop={'size':16})
		leg.get_frame().set_visible(False)
		ax.errorbar(x,zeros(len(bins_arr))+fidu_params[i], iyerr, fmt='m-', linewidth=2)
		plt.setp(ax.get_xticklabels(), visible=False)
		plt.subplots_adjust(hspace=0.0)
		
		
		
		if i == 1:
			ax2.set_yticks(np.linspace(0.005, 0.020,4))
			ax2.set_ylim(0.002, 0.024)
			
			savefig(plot_dir2+'MB_vs_bins_w_sim45_fw.pdf')
		else:
			ax2.set_yticks(np.linspace(0.001, 0.004, 3))
			ax2.set_ylim(0.0, 0.005)
			savefig(plot_dir2+'MB_vs_bins_Omsi_sim45_fw.pdf')
		close()

######## begin - plot err size vs number of bins for OmSi8 comb ###########
	#f = figure(figsize=(8,6))
	#ax=f.add_subplot(111)
	#for j in range(len(s_arr)):
		##ax.plot(x, y[i,:,j], markers[j/2], label='s = %s'%(s_arr[j]),linewidth=2, markersize=8)
		#ax.plot(x,y[0,:,j]**0.48*y[2,:,j],markers[j],label='s = %s'%(s_arr[j]),linewidth=1.5)
	#ax.set_xlabel('Number of bins',fontsize=16)
	#ax.set_ylabel(r'$\Omega_m^{0.48}\sigma_8$',fontsize=16)
	#leg=ax.legend(loc=0, ncol=1, labelspacing=.2, prop={'size':16})
	#leg.get_frame().set_visible(False)

	#show()

	##savefig(plot_dir2+'MB_vs_bins1000_OmSi8_sim45.jpg')
	##close()
	
######## begin - plot err size vs number of bins ###########
	#x = bins_arr
	#y = np.std(Tfit_ngal30[:,:,1:],axis=1).T

	#y /= sqrt(20000./12)
	#f = figure(figsize=(8,6))
	#ax=f.add_subplot(111)
	#markers=('go--','md-','r*:','c*','k-','md','yd','2','3','1','g-.')
	#for i in range(3):
		#ax.plot(x, y[i], markers[i], label=r'$%s$'%(label_latex[i]),linewidth=2, markersize=8)
	#ax.set_xlabel('Number of bins',fontsize=16)
	#ax.set_ylabel('$\sigma$',fontsize=16)
	#leg=ax.legend(loc=0, ncol=1, labelspacing=.2, prop={'size':16})
	#leg.get_frame().set_visible(False)
	##show()
	
	#savefig(plot_dir2+'err_vs_bins1000_fw.jpg')
	#close()
	######## end - plot err size vs number of bins ###########
	
if test_num_bins_cov:
	bins=30#700
	hist_bins = np.linspace (low, high, bins)
	peaks_mat_dir = '/Users/jia/magbias/peaks_mat_noise/'
	GetPeaksMat = lambda bins: genfromtxt(peaks_mat_dir+ 'peak_matrix_true_1000m-512b240_Om0.260_Ol0.740_w-1.000_ns0.960_si0.798_g6.976_ngal30_bins%s.ls'%(bins))
	test = GetPeaksMat(bins)
	f=figure()
	plot(hist_bins,average(test,axis=0),linewidth=2,drawstyle = 'steps-mid')
	show()

if peaks_err_vs_bins_old: 
	x = np.arange(10,110,10)
	y = np.std(Tfit_ngal30[:,:,1:],axis=1).T
	y /= sqrt(20000./12)
	f = figure(figsize=(8,6))
	ax=f.add_subplot(111)
	markers=('go--','md-','r*:','c*','k-','md','yd','2','3','1','g-.')
	for i in range(3):
		ax.plot(x, y[i], markers[i], label=r'$%s$'%(label_latex[i]),linewidth=2, markersize=8)
	ax.set_xlabel('Number of bins',fontsize=16)
	ax.set_ylabel('$\sigma$',fontsize=16)
	leg=ax.legend(loc=0, ncol=1, labelspacing=.2, prop={'size':16})
	leg.get_frame().set_visible(False)
	#show()
	savefig(plot_dir2+'err_vs_bins.pdf')
	close()
	
if shear_fft:
	mapdir = '/Users/jia/weaklensing/map_conv_shear_sample/'
	conv = filename_to_map  (mapdir+'WL-conv_m-512b240_Om0.260_Ol0.740_w-1.000_ns0.960_si0.798_4096xy_0001r_0029p_0100z_og.gre.fit')
	shear1 = filename_to_map(mapdir+'WL-shear1_m-512b240_Om0.260_Ol0.740_w-1.000_ns0.960_si0.798_4096xy_0001r_0029p_0100z_og.gre.fit')
	shear2 = filename_to_map(mapdir+'WL-shear2_m-512b240_Om0.260_Ol0.740_w-1.000_ns0.960_si0.798_4096xy_0001r_0029p_0100z_og.gre.fit')
	
	F_conv = fftpack.fftshift(fftpack.fft2(conv))
	F_shear1 = fftpack.fftshift(fftpack.fft2(shear1))
	F_shear2 = fftpack.fftshift(fftpack.fft2(shear2))
	
	## if we want to change by 45 degree
	#F_shear1 = fftpack.fftshift(fftpack.fft2(-shear2))
	#F_shear2 = fftpack.fftshift(fftpack.fft2(shear1))
	
	
	y, x = np.indices(F_shear1.shape)
	center = np.array([(x.max()-x.min())/2.0, (x.max()-x.min())/2.0])
	r = np.hypot(x - center[0], y - center[1])
	x, y = x - center[0], y - center[1]
	cosphi = x/r
	sinphi = y/r
	cos2phi = (cosphi**2-sinphi**2)
	sin2phi = (2*cosphi*sinphi)
	
	#lhs = F_shear1 / (cosphi**2-sinphi**2)
	#rhs = F_shear2 / (2*cosphi*sinphi)
	Emode = F_shear1*cos2phi+F_shear2*sin2phi
	Bmode = F_shear1*cos2phi-F_shear1*sin2phi
	#wrongEmode = F_shear2*(cosphi**2-sinphi**2)+F_shear1*(2*cosphi*sinphi)
	#wrongBmode = F_shear2*(cosphi**2-sinphi**2)-F_shear1*(2*cosphi*sinphi)
	def plotimshow(img,titlename, filename, vmin=-10,vmax=10):
		im=imshow(img, interpolation='nearest', origin='lower', vmin=vmin, vmax=vmax)
		im.set_extent([x[0,0],x[-1,-1], y[0,0],y[-1,-1]])
		title(titlename,fontsize=22)
		xlabel(r'$\ell_1$',fontsize=22)
		ylabel(r'$\ell_2$',fontsize=22)
		colorbar()
		savefig('/Users/jia/Dropbox/shear/'+filename)
		close()
		#savetxt('/Users/jia/Dropbox/shear/'+filename[:-4]+'.ls',img)
	
	plotimshow(real(Emode),'E mode','Emode.jpg')
	plotimshow(real(Bmode),'B mode','Bmode.jpg')
	plotimshow(real(Bmode)/real(Emode),'B/E','BEratio.jpg')
	
	#plotimshow(real(Bmode)/real(Emode),'B/E','BEratio2.jpg',vmin=-1,vmax=1)
	
	#plot power spectrum
	#Epsd2D = np.abs(Emode)**2
	#Bpsd2D = np.abs(Bmode)**2
	#x1, psd1DE, psd1DE2 = azimuthalAverage(Epsd2D)#
	#x1, psd1DB, psd1DB2 = azimuthalAverage(Bpsd2D)#x is ell_arr, psd1D is 
	#norm = ((2*pi*sqrt(12.00)/360.0)**2)/(2048.0**2)**2
	#f=figure()
	#ax=f.add_subplot(211)
	#ax.plot(x1,x1*(x+1)*psd1DE/2/pi*norm,label='E mode',linewidth=1)
	#ax.plot(x1,x1*(x+1)*psd1DB/2/pi*norm,label='B mode',linewidth=1)
	#ax.set_yscale('log')
	#ax.set_xscale('log')
	#ax.set_xlabel('l')
	#ax.set_ylabel(r'l(l+1)P(l)/2$\pi$')
	##show()
	#legend()
	#ax2=f.add_subplot(212)
	#ax2.plot(x,psd1DB/psd1DE,label='B/E',linewidth=1)
	#legend()
	##ax2.set_yscale('log')
	#ax2.set_xlabel('l')
	#ax2.set_xscale('log')
	#savefig('/Users/jia/Dropbox/shear/BEpowspec.jpg')
	#close()
	


	#plotimshow(real(F_shear1),r'$Re[\hat{\gamma_1}]$','real_shear1.jpg')
	#plotimshow(imag(F_shear1),r'$Im[\hat{\gamma_1}]$','imag_shear1.jpg')
	#plotimshow(real(F_shear2),r'$Re[\hat{\gamma_2}]$','real_shear2.jpg')
	#plotimshow(imag(F_shear2),r'$Im[\hat{\gamma_2}]$','imag_shear2.jpg')
	#plotimshow(real(lhs),r'$Re[\hat{\gamma_1}] / cos(2\phi)$','real_lhs.jpg')
	#plotimshow(imag(lhs),r'$Im[\hat{\gamma_1}] / cos(2\phi)$','imag_lhs.jpg')
	#plotimshow(real(rhs),r'$Re[\hat{\gamma_2}] / sin(2\phi)$','real_rhs.jpg')
	#plotimshow(imag(rhs),r'$Im[\hat{\gamma_2}] / sin(2\phi)$','imag_rhs.jpg')
	#plotimshow(real(lhs)/real(rhs),r'$Re [\hat{\gamma_1} / cos(2\phi)] / Re [\hat{\gamma_2} / sin(2\phi)$]','real_lr_ratio.jpg', -3, 3)
	#plotimshow(imag(lhs)/imag(rhs),r'$Im [\hat{\gamma_1} / cos(2\phi)]/ Im [\hat{\gamma_2} / sin(2\phi)$]', 'imag_lr_ratio.jpg', -3, 3)


if shear_conv_fft:
	mapdir = '/Users/jia/weaklensing/map_conv_shear_sample/'
	conv = filename_to_map  (mapdir+'WL-conv_m-512b240_Om0.260_Ol0.740_w-1.000_ns0.960_si0.798_4096xy_0001r_0029p_0100z_og.gre.fit')
	shear1 = filename_to_map(mapdir+'WL-shear1_m-512b240_Om0.260_Ol0.740_w-1.000_ns0.960_si0.798_4096xy_0001r_0029p_0100z_og.gre.fit')
	shear2 = filename_to_map(mapdir+'WL-shear2_m-512b240_Om0.260_Ol0.740_w-1.000_ns0.960_si0.798_4096xy_0001r_0029p_0100z_og.gre.fit')
	
	F_conv = fftpack.fftshift(fftpack.fft2(conv))
	F_shear1 = fftpack.fftshift(fftpack.fft2(shear1))
	F_shear2 = fftpack.fftshift(fftpack.fft2(shear2))
	
	#F_conv = fftpack.fft2(conv)
	#F_shear1 = fftpack.fft2(shear1)
	#F_shear2 = fftpack.fft2(shear2)
	
	psd2D_conv = np.abs(F_conv)**2
	psd2D_shear = np.abs(F_shear1)**2+np.abs(F_shear2)**2
	
	hist((psd2D_shear/psd2D_conv).flatten(),range=(-0.1,3),bins=100)
	xlabel(r'[$\hat{\gamma_1}^2 + \hat{\gamma_2}^2] / \hat{\kappa}^2$')
	show()
	#savefig(plot_dir2+"shear_conv_fft.jpg")



if test_dN_noise_nonoise:
	contr_mat=np.genfromtxt("/Users/jia/Documents/weaklensing/magbias/ps_mat/ps_matrix_true_1000m-512b240_Om0.260_Ol0.740_w-1.000_ns0.960_si0.798.ls")

	contr_mat_noise=np.genfromtxt("/Users/jia/Documents/weaklensing/magbias/ps_mat_noise/ps_matrix_true_1000m-512b240_Om0.260_Ol0.740_w-1.000_ns0.960_si0.798_ngal15.ls")

	hi_si_mat=np.genfromtxt("/Users/jia/Documents/weaklensing/magbias/ps_mat/ps_matrix_true_1000m-512b240_Om0.260_Ol0.740_w-1.000_ns0.960_si0.850.ls")

	hi_si_mat_noise=np.genfromtxt("/Users/jia/Documents/weaklensing/magbias/ps_mat_noise/ps_matrix_true_1000m-512b240_Om0.260_Ol0.740_w-1.000_ns0.960_si0.850_ngal15.ls")

	dN_nonoise=contr_mat-hi_si_mat
	dN_noise=contr_mat_noise-hi_si_mat_noise

	fidu_dN_nonoise=average(contr_mat,axis=0)-average(hi_si_mat,axis=0)
	fidu_dN_noise=average(contr_mat_noise,axis=0)-average(hi_si_mat_noise,axis=0)

	difference=np.min(abs((fidu_dN_nonoise/fidu_dN_noise-1)[:203]))
	print difference

if find_CFHT_err:
	alpha=0.59
	def prepareplot(fits_fidu,alpha):

		w_arr=fits_fidu.T[2]
		w_origin=average(w_arr)
		os_origin=average(fits_fidu.T[1])**alpha*average(fits_fidu.T[3])
		os_arr=fits_fidu.T[1]**alpha*fits_fidu.T[3]
		P=np.cov(w_arr[~np.isnan(os_arr)],os_arr[~np.isnan(os_arr)])
	
		return w_origin, os_origin, P

	fits_pow_fidu=np.genfromtxt(ps_fit_dir+'powspec_bias_analytical_s0.4_ngal15.ls')
	fits_pow_s02 =np.genfromtxt(ps_fit_dir+'powspec_bias_analytical_s0.2_ngal15.ls')
			
	fits_peaks_s02_15=genfromtxt(fit_analytical_dir+'bias_analytical_s0.2_g6.976_ngal15_bins200.ls')
	fits_peaks_fidu_15=genfromtxt(fit_analytical_dir+'true_analytical_g6.976_ngal15_bins200.ls')
	
	fits_peaks_s02_45=genfromtxt(fit_analytical_dir+'bias_analytical_s0.2_g6.976_ngal45_bins200.ls')
	fits_peaks_fidu_45=genfromtxt(fit_analytical_dir+'true_analytical_g6.976_ngal45_bins200.ls')

	
	def find_sigmas(fidu_fits,bias_fits,area,alpha):
		om_arr,w_arr,si_arr=fidu_fits.T[1:]
		w_std=np.std(w_arr)/sqrt(area/12.)
		omsi_std=np.std((om_arr/.26)**alpha*si_arr)/sqrt(area/12.)
		
		#omsi_std=np.std((om_arr)**alpha*si_arr)/sqrt(area/12.)
		
		om,w,si=average(bias_fits,axis=0)[1:]
		omsi0=0.26**alpha*0.798
		print 'w_std, omsi_std', w_std, omsi_std
		return 'dw, domsi = %.2f %.2f'%( (w+1)/w_std, (om**alpha*si-omsi0)/omsi_std)
		
	print 'CFHT powspec', find_sigmas(fits_pow_fidu,fits_pow_s02, 154.*.75,0.59)
	
	print 'CFHT peaks', find_sigmas(fits_peaks_fidu_15,fits_peaks_s02_15, 154.*.75,0.48)
	
	print 'COSMOS powspec', find_sigmas(fits_pow_fidu,fits_pow_s02, 1.0,0.62)
	
	print 'COSMOS peaks', find_sigmas(fits_peaks_fidu_45,fits_peaks_s02_45, 1.0,0.48)
	
	
	
	##fits_pow_s08 =np.genfromtxt(ps_fit_dir+'powspec_bias_analytical_s0.8_ngal15.ls')
	##os_arr=fits_pow_fidu.T[1]**alpha*fits_pow_fidu.T[3]
	#os_arr=(fits_pow_fidu.T[1]/0.26)**alpha*fits_pow_fidu.T[3]
	#os_arr=os_arr[~np.isnan(os_arr)]
	#os_std= np.std(os_arr)/sqrt(154.*.75/12)

	#os_fidu = prepareplot(fits_pow_fidu)[1]
	#os_s02  = prepareplot(fits_pow_s02 )[1]
	#os_s08  = prepareplot(fits_pow_s08 )[1]
	
	##print 'our std (CFHT=0.3) = ', os_std
	##print 'CFHTLS sigmas', (os_s02-os_fidu)/0.03	

if test_CFHT_err:
	#test if we limite our ell to up to 4000, error size would increase by 2
	ngal=15 ## for noise maps
	num=203 # the index at which ell ends
	
	bins = 1000
	kmin=1.0
	kmax=2048./2
	binsize=(kmax-kmin)/1000.0
	ell_arr = (np.arange(kmin,kmax,binsize)*360./sqrt(12.0))[:1000]
	bins=1000
	########## grab relevant peak distributions #####
	ps_mat=batch_kmap_powspec(kmap_dir,contr,bias=0,ngal=ngal)[:,:num]
	#ps_mat=batch_kmap_powspec(kmap_dir,contr,bias=0)[:,:num] #noiseless
	
	fidu_mat_fn=ps_mat_dir+'fidu_true_mat_ngal%s.ls'%(ngal)
	fidu_mat=genfromtxt(fidu_mat_fn)[:,:num]
	
	cov_mat=np.cov(ps_mat,rowvar=0) # for noise maps
	cov_inv=np.mat(cov_mat).I
	#cov_inv = np.genfromtxt('cov_inv_CFHT.ls') # for noiseless
	
	dNOm=fidu_mat[hi_Om]-fidu_mat[contr]
	dNw =fidu_mat[hi_w ]-fidu_mat[contr]
	dNsi=fidu_mat[hi_si]-fidu_mat[contr]
	dp = np.array([0.03,0.2,0.052])
	dNdOm=dNOm/dp[0]
	dNdw =dNw/dp[1] 
	dNdsi=dNsi/dp[2]
	X=np.mat([dNdOm,dNdw,dNdsi])
	
	fit_filename='/Users/jia/Documents/weaklensing/magbias/CFHT_test_fit_ell%s.ls'%(int(ell_arr[num]))
	fit_analytical=np.ndarray(shape=(1000,4))
	for i in range(1000):
		#print 'analytical fitting to map #',i
		ibmap=ps_mat[i]
		Y=np.mat(ibmap-fidu_mat[contr])
		del_p=((X*cov_inv*X.T).I)*(X*cov_inv*Y.T)		
		initguess=np.squeeze(np.array(del_p.T))+fidu_params
		fit_analytical[i,0]=0#chisq(initguess, ibmap)
		fit_analytical[i,1:]=np.array(initguess)
	savetxt(fit_filename,fit_analytical)
	Om_arr,si_arr=fit_analytical.T[[1,3]]
	os_arr=(Om_arr/0.26)**0.59*si_arr#(Om_arr/0.27)**0.59*si_arr
	os_std= np.std(os_arr[~np.isnan(os_arr)])/sqrt(154.*.75/12)
	print os_std


if test_actual_fit_powspec:
	ngal=30
	bins=203
	fidu_mat_fn=ps_mat_dir+'fidu_true_mat_ngal%s.ls'%(ngal)
	if os.path.isfile(fidu_mat_fn):
		fidu_mat=np.genfromtxt(fidu_mat_fn)
	else:
		true_mat = np.ndarray(shape=(len(cosmo_all), kmapnum, bins))

		for i in np.arange(len(cosmo_all)):
			true_mat[i] = batch_kmap_powspec(kmap_dir,i, s=s, bias=0, ngal=ngal)

		fidu_mat = average(true_mat,axis=1)
		savetxt(fidu_mat_fn,fidu_mat)
	ibmap01=genfromtxt(ps_mat_dir+"ps_matrix_bias_1000m-512b240_Om0.260_Ol0.740_w-1.000_ns0.960_si0.798_s0.1_ngal%s.ls"%(ngal))[0,:bins]
	ibmap08=genfromtxt(ps_mat_dir+"ps_matrix_bias_1000m-512b240_Om0.260_Ol0.740_w-1.000_ns0.960_si0.798_s0.8_ngal%s.ls"%(ngal))[0,:bins]
	ps_params_01=genfromtxt(ps_fit_dir+"powspec_bias_analytical_s0.1_ngal%s.ls"%(ngal))[0]
	ps_params_08=genfromtxt(ps_fit_dir+"powspec_bias_analytical_s0.8_ngal%s.ls"%(ngal))[0]
	
	kmin=1.0
	kmax=2048./2
	binsize=(kmax-kmin)/1000.0
	ell_arr = (np.arange(kmin,kmax,binsize)*360./sqrt(12.0))[:bins]
	
	iOm01, iw01, isi01=ps_params_01[1:]
	ps_fit01=interp_distribution('fw', fidu_mat, iOm01, iw01, isi01,bins=bins)
	iOm08, iw08, isi08=ps_params_08[1:]
	ps_fit08=interp_distribution('fw', fidu_mat, iOm08, iw08, isi08,bins=bins)
	
	f=figure(figsize=(6,8))
	ax1=f.add_subplot(211)
	ax2=f.add_subplot(212,sharex=ax1)
	ax1.plot(ell_arr,ps_fit01,'r-',label='s=0.1(fit)')
	ax1.plot(ell_arr,ps_fit08,'b-',label='s=0.8(fit)')
	ax1.plot(ell_arr,ibmap01,'r--',label='s=0.1(original)')
	ax1.plot(ell_arr,ibmap08,'b--',label='s=0.8(original)')
	
	ax2.plot(ell_arr,ps_fit01/ibmap01-1,'r-',label='s=0.1')
	ax2.plot(ell_arr,ps_fit08/ibmap08-1,'b-',label='s=0.8')
	leg=ax1.legend(loc=0, ncol=1, labelspacing=.2, prop={'size':12})
	leg.get_frame().set_visible(False)
	ax1.set_yscale('log')
	ax1.set_xscale('log')
	ax2.set_xscale('log')
	ax2.set_xlabel(r'$\ell$')
	ax1.set_ylabel(r'$\ell(\ell+1)P(\ell)/2\pi$')
	ax2.set_ylabel(r'${\Delta}P/P$',fontsize=16)
	ax1.set_xlim(1e2,1e4)
	ax1.set_ylim(1e-5,1e-3)
	show()
	


if analytical_chisqmin:

	bins=30
	ngal=30
	
	peak_mat_true=np.ndarray(shape=(7,1000,bins))
	fidu_mat=np.ndarray(shape=(7,30))
	for i in range(7):
		peak_mat_temp=batch_kmap_noise(kmap_dir,i, bins=bins, bias=0,ngal=30)
		peak_mat_true[i]=peak_mat_temp
		f=average(peak_mat_temp,axis=0)
		fidu_mat[i]=f
	cov_mat=np.cov(peak_mat_true[contr],rowvar=0)
	cov_inv=np.mat(cov_mat).I
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
	X=np.mat([dNdOm,dNdw,dNdsi])
	
	peak_mat_bias=batch_kmap_noise(kmap_dir,contr,s=1.5,bias=1,bins=bins,ngal=ngal)
	bias_fit_py=np.genfromtxt(fit_noise_dir+'bias_fit_s1.5_g6.976_ngal30_bins30.ls')
	true_fit_py=np.genfromtxt(fit_noise_dir+'true_fit_g6.976_ngal30_bins30.ls')
		
	def chisq_analytical (initguess, ibmap):
		iOm, iw, isi = initguess
		Ni_arr=interp_distribution(interpmethod,fidu_mat,iOm, iw, isi)
		del_Ni = np.mat(Ni_arr - ibmap)
		chisquare = float64(del_Ni * cov_inv * del_Ni.T)
	#	print chisquare, initguess
		return chisquare
	
	def analytical_fits(peak_mat,cov_inv,X):
		fit_analytical=np.ndarray(shape=(1000,4))
		for i in range(1000):
			Y=np.mat(peak_mat[i]-fidu_mat[contr])
			del_p=((X*cov_inv*X.T).I)*(X*cov_inv*Y.T)
			
			initguess=np.squeeze(np.array(del_p.T))+fidu_params
			fit_analytical[i,0]=chisq_analytical(initguess,peak_mat[i])
			fit_analytical[i,1:]=np.array(initguess)
		return fit_analytical
	bias_fit_analytical=analytical_fits(peak_mat_bias,cov_inv,X)
	true_fit_analytical=analytical_fits(peak_mat_true[contr],cov_inv,X)
	
	##############find alpha###########
	#Oms,sis=true_fit_analytical[:,[1,3]].T
	#alpha_residue=lambda alpha: sum(abs(alpha*(Oms-0.26)/0.26+(sis-0.798)/0.798+alpha*(Oms-0.26)/0.26*(sis-0.798)/0.798))
	
	
	#alpha_minimize=optimize.minimize(alpha_residue,0.4,method='L-BFGS-B')
	#alpha_x=alpha_minimize.x
	#####################
	
	
	########### plot fit and analytical scatter plot ############
	cov_true=np.cov(true_fit_analytical[:,1:].T)	
	eig_true=linalg.eig(cov_true)
	P=np.mat(eig_true[1])
	A=np.mat(cov_true)
	D=P.I*A*P
	R=D*A.I#RA=D
	
	#f=figure(figsize=(8,6))
	#ax=subplot(111)
	#ax.scatter(true_fit_analytical[:,1],true_fit_analytical[:,3],c='b',marker='o',label='analytical $<\chi^2>$=%.1f'%(average(true_fit_analytical[:,0])))
	#ax.scatter(true_fit_py[:,1],true_fit_py[:,3],c='r',marker='o',label='fit $<\chi^2>$=%.1f'%(average(true_fit_py[:,0])))
		
	#ax.legend(loc=0, scatterpoints = 1)
	#ax.set_xlabel('$\Omega_m$')
	
	#ax.set_ylabel('$\sigma_8$')
	#ax.set_title(interpmethod)
	
	#savefig(plot_dir2+'analytical_vs_fit_%s_correct.jpg'%(interpmethod))
	#close()
	
	######################## plot w vs Om ##########
	##ax.scatter(true_fit_analytical[:,1],true_fit_analytical[:,2],c='b',marker='o',label='analytical $<\chi^2>$=%.1f'%(average(true_fit_analytical[:,0])))
	##ax.scatter(true_fit_py[:,1],true_fit_py[:,2],c='r',marker='o',label='fit $<\chi^2>$=%.1f'%(average(true_fit_py[:,0])))
	##ax.set_ylabel('$w$')
	##savefig(plot_dir2+'analytical_vs_fit_%s_w.jpg'%(interpmethod))
	#close()
	
	##find_alpha2=lambda alpha:alpha*sqrt(cov_true[0,0])/0.26+sqrt(cov_true[2,2])/0.798-sqrt(D[2,2])
	##mata=np.mat(cov_true[::2,::2])
	##matb=np.mat(D[::2,::2])
	##c=np.mat([0,0,1])
	##P=np.mat(eig_true[1])
	
if fw_bw_true_scatter:
	bins=30
	ngal=30
	s=0.4
	fw_fit=np.genfromtxt('/Users/jia/Documents/weaklensing/magbias/fit_noise_analytical_seed0/fw/bias_analytical_s%s_g6.976_ngal%s_bins%s.ls'%(s,ngal,bins))
	bw_fit=np.genfromtxt('/Users/jia/Documents/weaklensing/magbias/fit_noise_analytical_seed0/bw/bias_analytical_s%s_g6.976_ngal%s_bins%s.ls'%(s,ngal,bins))
	
	############ plot the scatter plot individually for fw/bw #####
	#f = figure(figsize=(8,6))
	#ax=f.add_subplot(111)
	#ax.scatter(fw_fit[:,1],fw_fit[:,3],c='b',marker='o',label='fw')
	#ax.scatter(bw_fit[:,1],bw_fit[:,3],c='r',marker='o',label='bw')
	#ax.legend(loc=0, scatterpoints = 1)
	#ax.set_xlabel('$\Omega_m$')
	##ax.set_ylabel('$w$')
	#ax.set_ylabel('$\sigma_8$')
	#savefig(plot_dir2+'fw_bw_true_scatter_w.jpg')
	#close()
	
	############# plot the offset between fw&bw as function of Om ######
	#f = figure(figsize=(8,6))
	#ax=f.add_subplot(111)
	#ax.scatter(fw_fit[:,1]-bw_fit[:,1],fw_fit[:,3]-bw_fit[:,3],c='b',marker='o',label='fw-bw')
	#ax.legend(loc=0, scatterpoints = 1)
	#ax.set_xlabel('$\Delta\Omega_m$')
	##ax.set_ylabel('$w$')
	#ax.set_ylabel('$\Delta\sigma_8$')
	##show()
	#savefig(plot_dir2+'fw_bw_del_true_scatter.jpg')
	#close()
	
	############# take into account of degeneracy ######
	#f = figure(figsize=(8,6))
	#ax=f.add_subplot(111)
	#ax.scatter(fw_fit[:,1]**0.4*fw_fit[:,3],bw_fit[:,1]**0.4*bw_fit[:,3],c='b',marker='o',label='fw-bw')
	##ax.legend(loc=0, scatterpoints = 1)
	#ax.set_xlabel('$\Omega_m^{0.6}\sigma_8$ fw')
	##ax.set_ylabel('$w$')
	#ax.set_ylabel('$\Omega_m^{0.6}\sigma_8$ bw')
	##show()
	#savefig(plot_dir2+'fw_bw_degeneracy_true_scatter.jpg')
	#close()
	
	
	############# plot the offset between fw&bw as function of Om ######
	
	for i in range(3):
		f = figure(figsize=(8,6))
		ax=f.add_subplot(111)
		ax.plot(fw_fit[:,i+1],bw_fit[:,i+1],'bo')
		ax.set_xlabel('fw')
		ax.set_ylabel('bw')
		ax.set_title(label_str[i])
		savefig(plot_dir2+'fw_bw_%s.jpg'%(label_str[i]))
		close()

if cosmo_truemap:
	def getfitsmap(fitsfile):
		hdulist = pyfits.open(fitsfile)
		kmap = np.array(hdulist[0].data)
		return kmap
	
	vmin,vmax=-0.034,0.54
	
	contr_kmap=getfitsmap(kmap_dir+cosmo_all[contr]+'/WL-conv_m-512b240_Om0.260_Ol0.740_w-1.000_ns0.960_si0.798_4096xy_0001r_0029p_0100z_og.gre.fit')
	hi_w_kmap=getfitsmap(kmap_dir+cosmo_all[hi_w]+'/WL-conv_m-512b240_Om0.260_Ol0.740_w-0.800_ns0.960_si0.798_4096xy_0001r_0029p_0100z_og.gre.fit')
	lo_si_kmap=getfitsmap(kmap_dir+cosmo_all[lo_si]+'/WL-conv_m-512b240_Om0.260_Ol0.740_w-1.000_ns0.960_si0.750_4096xy_0001r_0029p_0100z_og.gre.fit')
	
	img_labels=('contr','hi_w','lo_si')
	i=0
	#for img in (contr_kmap,hi_w_kmap,lo_si_kmap):
		#f=figure(figsize=(8,8))
		#ax=f.add_subplot(111)
		#imshow(img,norm=None,interpolation='nearest',vmin=vmin,vmax=vmax)
		#plt.setp(ax.get_xticklabels(), visible=False)
		#plt.setp(ax.get_yticklabels(), visible=False)
		#plt.subplots_adjust(top=1.0,bottom=0.0,left=0.0,right=1.0)
		#savefig(plot_dir2+'kmap_%s.jpg'%(img_labels[i]))
		#close()
		
		#f=figure(figsize=(8,8))
		#ax=f.add_subplot(111)
		#imshow(smooth(img,sigma),norm=None,interpolation='nearest',vmin=vmin,vmax=vmax)
		#plt.setp(ax.get_xticklabels(), visible=False)
		#plt.setp(ax.get_yticklabels(), visible=False)
		#plt.subplots_adjust(top=1.0,bottom=0.0,left=0.0,right=1.0)
		#savefig(plot_dir2+'kmap_smooth_%s.jpg'%(img_labels[i]))
		#close()		
		
		#i+=1
	
	def imshow_kmap(kmap,filename):
		f=figure(figsize=(8,8))
		ax=f.add_subplot(111)
		imshow(kmap,norm=None,interpolation='nearest')#,vmin=vmin,vmax=vmax)
		plt.setp(ax.get_xticklabels(), visible=False)
		plt.setp(ax.get_yticklabels(), visible=False)
		plt.subplots_adjust(top=1.0,bottom=0.0,left=0.0,right=1.0)
		savefig(plot_dir2+filename)
		close()
	
	ngal=15
	z=1.0
	gamma=0.15+0.035*z
	noise_width=sqrt((gamma*gamma)/ngal)
	sigma = 6.976
	noise = np.array(np.random.normal(0.0, noise_width,contr_kmap.shape),dtype=float32)
	
	imshow_kmap(contr_kmap,'sample_raw_map.jpg')
	imshow_kmap(noise,'sample_noise_map.jpg')
	imshow_kmap(biasnoisesmooth(contr_kmap,noise,sigma,1.5),'sample_biasnoisesmooth_map.jpg')
	
	imshow_kmap(bias_nosmooth(contr_kmap,1.5),'sample_bias_map.jpg')
		
	
if individual_components:
	
	bins=30
	ngal=30
	bin_size = (high - low) / bins
	hist_bins = np.linspace (low+bin_size/2, high-bin_size/2, bins)
	markers=('k-','b--','r-','b-','g-','m-d','yd','2','3','1','g-.')
	
	fidu_mat=np.ndarray(shape=(7,30))
	for i in range(7):
		peak_mat_temp=batch_kmap_noise(kmap_dir,i, bins=bins, bias=0,ngal=30)
		f=average(peak_mat_temp,axis=0)
		fidu_mat[i]=f
	
	for interpmethod in ('spline','fw','bw'):
		for s in (0.1,-0.5,1.5,2.0):
			bias_fit_temp=np.genfromtxt('/Users/jia/magbias/fit_noise_seed73000/%s/bias_fit_s%s_g6.976_ngal30_bins30.ls'%(interpmethod,s))
			
			bias_hist=batch_kmap_noise(kmap_dir,contr,s=s,bins=bins, bias=1,ngal=30)
			bias_hist=bias_hist[0]
			
			bias_fit_temp=bias_fit_temp[0]
			
			#average(bias_hist,axis=0)
			#average(bias_fit_temp,axis=0)
			
			iOm, iw, isi = bias_fit_temp[1:]
			
			bias_all = interp_distribution(interpmethod,fidu_mat,iOm,iw,isi)
			bias_iOm = interp_distribution(interpmethod,fidu_mat,iOm,-1.0,0.798)
			bias_iw = interp_distribution(interpmethod,fidu_mat,0.26,iw,0.798)
			bias_isi = interp_distribution(interpmethod,fidu_mat,0.26,-1.0,isi)
			
			f=figure(figsize=(8,6))
			gs = gridspec.GridSpec(2,1,height_ratios=[3,1])
			ax1=f.add_subplot(gs[0])
			ax2=f.add_subplot(gs[1],sharex=ax1)
			
			ax1.plot(hist_bins,fidu_mat[contr],'k-',linewidth=1.5, label=labels[contr])
			ax1.plot(hist_bins,bias_hist,'k--',linewidth=1.5, label='Bias')
			ax2.plot(hist_bins,zeros(len(hist_bins)),'k-',linewidth=1.5,label='Bias')
			ax2.plot(hist_bins,bias_hist/fidu_mat[contr]-1,'k--',linewidth=1.5)
			
			temp_labels=('Fit all (%.3f %.2f %.3f)'%(iOm,iw,isi),'Om only','w only','si only')
			i=0
			for bias_arr in (bias_all,bias_iOm,bias_iw,bias_isi):
				ax1.plot(hist_bins,bias_arr,markers[i+1],linewidth=1.5, label=temp_labels[i])
				ax2.plot(hist_bins,bias_arr/fidu_mat[contr]-1,markers[i+1],linewidth=1.5)
				i+=1
			
			ax1.set_title('s=%s bins%s, ngal%s,%s'%(s,bins,ngal,interpmethod))
			plt.setp(ax1.get_xticklabels(), visible=False)
			plt.subplots_adjust(hspace=0.0)
			leg=ax1.legend(loc=0, ncol=1, labelspacing=.2, prop={'size':14})
			leg.get_frame().set_visible(False)
			ax1.set_ylabel('peak counts $N$($\kappa$)')
			ax2.set_ylabel(r'${\Delta}N/N$')
			ax2.set_ylim(-0.5,0.5)
			ax2.set_xlim(-0.02,0.10)
			plt.setp(ax1.get_xticklabels(), visible=False)
			savefig(plot_dir2+'individual_components_1map_s%s_%s.jpg'%(s,interpmethod))
		

if destroy_high_peaks:
	
	ngal=45
	bins=100
	s_arr4=(0.0,1.0,1.5,2.0)
	bin_size = (high - low) / bins
	hist_bins = np.linspace (low+bin_size/2, high-bin_size/2, bins)
	plot_dir2 = '/Users/jia/Documents/weaklensing/magbias/plot/'
	destroy_dir='/Users/jia/magbias/destroy_hipeaks/'
	for ngal in ngal_arr:
		c=np.genfromtxt(destroy_dir+'peak_matrix_true_1000m-512b240_Om0.260_Ol0.740_w-1.000_ns0.960_si0.798_g6.976_bins100.ls')
		d=np.genfromtxt(destroy_dir+'peak_matrix_true_1000m-512b240_Om0.260_Ol0.740_w-1.000_ns0.960_si0.798_g6.976_ngal%s_bins100.ls'%(ngal))
		peaks_true=average(c,axis=0)
		peaks_true_noise=average(d,axis=0)
		
		peaks_bias=np.ndarray(shape=(len(s_arr4),bins))
		peaks_bias_noise=np.ndarray(shape=(len(s_arr4),bins))
		
		i=0
		for s in s_arr4:		
			a=np.genfromtxt(destroy_dir+'peak_matrix_bias_1000m-512b240_Om0.260_Ol0.740_w-1.000_ns0.960_si0.798_s%s_g6.976_bins100.ls'%(s))
			b=np.genfromtxt(destroy_dir+'peak_matrix_bias_1000m-512b240_Om0.260_Ol0.740_w-1.000_ns0.960_si0.798_s%s_g6.976_ngal%s_bins100.ls'%(s,ngal))
			
			peaks_bias[i]=average(a,axis=0)
			peaks_bias_noise[i]=average(b,axis=0)
			i+=1
			
		####### plot bias - true, useless ###
		#f = figure(figsize=(8,6))
		#gs = gridspec.GridSpec(2,1,height_ratios=[3,1])
		#ax1=f.add_subplot(gs[0])
		#ax2=f.add_subplot(gs[1],sharex=ax1)
		#ax1.plot(hist_bins,peaks_true_noise-peaks_true,linewidth=1.5, label='Fiducial (Noiseless)')
		#ax2.plot(hist_bins,zeros(bins),linewidth=1.5)
		#for n in range(len(s_arr4)):
			#ax1.plot(hist_bins,peaks_bias_noise[n]-peaks_bias[n],linewidth=1.5, label='Bias s=%s'%(s_arr4[n]))
			#ax2.plot(hist_bins,(peaks_bias_noise[n]-peaks_bias[n])/(peaks_true_noise-peaks_true)-1,linewidth=1.5, label='Bias s=%s'%(s_arr4[n]))
		#ax2.set_xlabel('convergence $\kappa$')
		#ax1.set_ylabel('$N_{noise}-N_{noiseless}$')
		#leg1=ax1.legend(loc=0, ncol=1, labelspacing=.2, prop={'size':12})
		#leg1.get_frame().set_visible(False)
		#ax1.set_xlim(-0.02,0.1)
		#ax2.set_ylim(-1,1)
		#plt.setp(ax1.get_xticklabels(), visible=False)
		#plt.subplots_adjust(hspace=0.0)
		#savefig(plot_dir2+'destroypeaks_def_ngal%s.jpg'%(ngal))
		#close()
		
		f = figure(figsize=(12,6))
		gs = gridspec.GridSpec(2,2,height_ratios=[3,1])
		ax1=f.add_subplot(gs[0])
		ax2=f.add_subplot(gs[1])
		ax3=f.add_subplot(gs[2],sharex=ax1)
		ax4=f.add_subplot(gs[3],sharex=ax2)
		
		ax1.plot(hist_bins,peaks_true,linewidth=1.5, label='Fiducial (Noiseless) %s' %(int(sum(peaks_true))))
		ax2.plot(hist_bins,peaks_true_noise,linewidth=1.5, label='Fiducial (Noise)%s'%(int(sum(peaks_true_noise))))
		ax3.plot(hist_bins,zeros(bins),linewidth=1.5)
		ax4.plot(hist_bins,zeros(bins),linewidth=1.5)
		for n in range(len(s_arr4)):
			ax1.plot(hist_bins,peaks_bias[n],linewidth=1.5, label='Bias s=%s (%s)'%(s_arr4[n], int(sum(peaks_bias[n]))))
			ax2.plot(hist_bins,peaks_bias_noise[n],linewidth=1.5, label='Bias s=%s (%s)'%(s_arr4[n], int(sum(peaks_bias_noise[n]))))	
			ax3.plot(hist_bins,peaks_bias[n]/peaks_true-1,linewidth=1.5, label='Bias s=%s'%(s_arr4[n]))
			ax4.plot(hist_bins,peaks_bias_noise[n]/peaks_true_noise-1,linewidth=1.5, label='Bias s=%s'%(s_arr4[n]))
		leg1=ax1.legend(loc=0, ncol=1, labelspacing=.2, prop={'size':12})
		leg1.get_frame().set_visible(False)
		leg2=ax2.legend(loc=0, ncol=1, labelspacing=.2, prop={'size':12})
		leg2.get_frame().set_visible(False)
		
		
		ax3.set_xlabel('convergence $\kappa$')
		ax1.set_xlabel('convergence $\kappa$')
		
		ax1.set_ylabel('peak counts $N$($\kappa$)')
		ax3.set_ylabel(r'${\Delta}N/N$')
		
		ax1.set_xlim(-0.02,0.1)
		ax2.set_xlim(-0.02,0.1)
		ax3.set_ylim(-1,1)
		ax4.set_ylim(-0.5,0.5)
		
		plt.setp(ax1.get_xticklabels(), visible=False)	
		plt.setp(ax2.get_xticklabels(), visible=False)
		plt.subplots_adjust(hspace=0.0)
		
		savefig(plot_dir2+'destroypeaks_ngal%s.jpg'%(ngal))
		close()


	
if destroy_high_peaks_singlemap:
	from findpeak import findpeak
	import findpeak3 as fp2
	fitsfile="/Users/jia/Documents/weaklensing/magbias/test.fit"
	
	ngal_arcmin=30
	ngal=12.0*60**2*ngal_arcmin/2048.0**2
	z=1.0
	gamma=0.15+0.035*z
	noise_width=sqrt((gamma*gamma)/ngal)
	sigma = 6.976
	
	low, high = -0.02, 0.19
	bins=100
	bin_size = (high - low) / bins
	hist_bins = np.linspace (low+bin_size/2, high-bin_size/2, bins)
	
	s_arr3=(1.0,2.0)#np.arange(-0.5,2.5,0.5)
	biasmooth = lambda kmap, sigma, s: snd.filters.gaussian_filter(kmap*(1+(5*s-2)*kmap),sigma)/snd.filters.gaussian_filter(1+(5*s-2)*kmap,sigma)
	
	biasnoisesmooth = lambda kmap, noise, sigma, s: snd.filters.gaussian_filter(kmap*(1+(5*s-2)*kmap)+noise,sigma)/snd.filters.gaussian_filter(1+(5*s-2)*kmap,sigma)
	
	smooth = lambda kmap, sigma: snd.filters.gaussian_filter(kmap,sigma)
		
	hdulist = pyfits.open(fitsfile)
	kmap = np.array(hdulist[0].data)
	seed(0)#seed(73001)
	noise = np.array(np.random.normal(0.0, noise_width, kmap.shape),dtype=float32)
	
	kmap_GRF=smooth(noise,sigma)	
	kmap_noise = kmap + noise
	kmap_true = smooth (kmap, sigma)	
	kmap_true_noise = smooth (kmap_noise, sigma)
	
	s=1.5
	kmap_bias = biasmooth(kmap, sigma, s)
	kmap_bias_noise = biasnoisesmooth(kmap,noise,sigma,s)	
	
	 
	def peakloc(kmap):
		positions = fp2.findpeak(kmap)[1]
		return positions
		
	peak_loc_true = peakloc(kmap_true)
	peak_loc_true_noise = peakloc(kmap_true_noise)
	peak_loc_bias = peakloc(kmap_bias)
	peak_loc_bias_noise = peakloc(kmap_bias_noise)
	peak_loc_GRF = peakloc(kmap_GRF)
	
	kappa_diff=kmap_bias_noise*peak_loc_true_noise-kmap_true_noise*peak_loc_true_noise
	
	def returnpeakhist(kmap,bins=bins):
		peaks = findpeak(kmap)
		peakhist=np.histogram(peaks, range=(low,high), bins=bins)[0]
		return peakhist
		
	#peaks_due_to_GRF = peak_loc_true_noise*peak_loc_GRF
	all_true_noise_map = kmap_true_noise*peak_loc_true_noise
	remain_loc_noise = peak_loc_true_noise*peak_loc_bias_noise
	remain_true_noise_map = kmap_true_noise*remain_loc_noise
	remain_bias_noise_map = kmap_bias_noise*remain_loc_noise
	
	remain_true_noise = remain_true_noise_map[remain_true_noise_map.nonzero()]
	remain_bias_noise = remain_bias_noise_map[remain_bias_noise_map.nonzero()]
	true_noise_peaks = kmap_true_noise[peak_loc_true_noise.nonzero()]
	true_noise_peaks_in_bias = kmap_bias_noise[peak_loc_true_noise.nonzero()]
	
	####### paper plot1: scatter plot of all pixels before & after MB for paper ###
	
	#all_true = kmap_true_noise.reshape(1,-1).squeeze()
	#all_bias = kmap_bias_noise.reshape(1,-1).squeeze()
	
	#f = figure(figsize=(8,6))
	#ax=f.add_subplot(111)

	####### use heat map to plot density profile, not so good looking, so give up####
	##heatmap, xedges, yedges = np.histogram2d(all_bias, all_true, bins=100)#, range=[[-0.1,0.15],[-0.1,0.2]])
	###extent = [xedges[0], xedges[-1], yedges[0], yedges[-1]]
	##extent = [yedges[0], yedges[-1], xedges[0], xedges[-1]]
	##imshow(heatmap, interpolation='nearest', origin='lower',aspect='auto',cmap = 'binary',extent=extent)#[yedges[0], yedges[-1], xedges[0], xedges[-1]])
	####################
	
	#pixsample=randint(0, all_true.shape[0], 10000)
	
	#ax.scatter(all_true[pixsample],all_bias[pixsample],s=1,edgecolor='r',facecolor='r',label='All pixels (selected sample)')
	#ax.scatter(true_noise_peaks,true_noise_peaks_in_bias,s=1,edgecolor='g',facecolor='g',label='Peaks in true map')	
	#ax.scatter(remain_true_noise,remain_bias_noise,s=1,edgecolor='y',facecolor='y',label='Remaining peaks')
	#a=np.linspace(-0.1,0.2,100)
	#b=a*(1.0+(5*s-2.0)*a)
	#c=1.1*a
	
	#plot(a,a,'k-',linewidth=1.5,label='$\kappa_{true}=\kappa_{bias}$')
		
	#xlim(-0.065,0.15)
	#ylim(-0.065,0.15)
	##xlim(-0.1,0.2)
	##ylim(-0.105,0.25)
	#xlabel(r'$\kappa_{true}$',fontsize=16)
	#ylabel(r'$\kappa_{bias}$',fontsize=16)
	#leg=ax.legend(loc=2, ncol=1, labelspacing=.2, prop={'size':16})
	#leg.get_frame().set_visible(False)
	#savefig(plot_dir2+'all_pix_kappa_scatter.pdf')
	#close()
	##show()

	
	##### paper plot2: sample high/low peaks#####
	hipeak_loc=np.where((remain_true_noise_map>0.07)&(remain_true_noise_map<0.11))
	lopeak_disapear_loc=np.where((peak_loc_true_noise==1)&(remain_loc_noise==0))#(all_true_noise_map<-0.001)&
	peak_appear_loc=np.where((peak_loc_true_noise==0)&(peak_loc_bias_noise==1))
	psize=15

	lostpeakcount=0
	gainpeakcount=0
	lostx_arr=np.array([])
	losty_arr=np.array([])
	gainx_arr=np.array([])
	gainy_arr=np.array([])
	ss=5#shiftsize,within which we count as same peak
	for i in range(len(lopeak_disapear_loc[0])):#range(500):
		for hilo in ('appear','lo'):#'lo','hi',
			#print i, hilo
			#vmin, vmax=-0.02, 0.08
			
			if hilo=='hi':
				if i > len(hipeak_loc[0])-1:
					continue
				#else:
					#x,y=hipeak_loc[0][i],hipeak_loc[1][i]
					#filename=plot_dir2+'imshowpeak/sample_hipeak%s.pdf'%(i)
			elif hilo=='lo':
				if i > len(lopeak_disapear_loc[0])-1:
					continue
				x,y=lopeak_disapear_loc[0][i],lopeak_disapear_loc[1][i]
				if sum(peak_loc_bias_noise[x-ss+1:x+ss,y-ss+1:y+ss])==0:
					#print 'disapear',i
					lostpeakcount+=1
					lostx_arr=np.append(lostx_arr,x)
					losty_arr=np.append(losty_arr,y)
					filename=plot_dir2+'imshowpeak/sample_dispear%s_%s.jpg'%(i,ss)
					#filename=plot_dir2+'imshowpeak/sample_dispear%s.pdf'%(i)
				else:
					continue
			elif hilo=='appear':
				if i > len(peak_appear_loc[0])-1:
					continue
				x,y=peak_appear_loc[0][i],peak_appear_loc[1][i]
				if sum(peak_loc_true_noise[x-ss+1:x+ss,y-ss+1:y+ss])==0:
					#print 'appear',i
					gainpeakcount+=1
					gainx_arr=np.append(gainx_arr,x)
					gainy_arr=np.append(gainy_arr,y)
					filename=plot_dir2+'imshowpeak/sample_apear%s_%s.jpg'%(i,ss)
					#filename=plot_dir2+'imshowpeak/sample_dispear%s.pdf'%(i)
				else:
					continue

			
			if x+psize>2048 or x-psize<0 or y+psize>2048 or y-psize<0:
				print 'at map edge',hilo
				continue
			#else:
				##print 'plotting'
				#f = figure(figsize=(10,5))
			
				#ax1=f.add_subplot(121)
				#ax2=f.add_subplot(122)
				
				#truen_im=kmap_true_noise[x-psize:x+psize+1,y-psize:y+psize+1]
				#biasn_im=kmap_bias_noise[x-psize:x+psize+1,y-psize:y+psize+1]
				
				#vmin, vmax=np.amin(truen_im),np.amax(biasn_im)
				
				#im=ax1.imshow(truen_im,norm=None,interpolation='nearest',vmin=vmin,vmax=vmax)
				
				#ax2.imshow(biasn_im,norm=None,interpolation='nearest',vmin=vmin,vmax=vmax)
				
				#cax = f.add_axes([0.83, 0.125, 0.03, 0.75])#from left, bottom, width
				#f.colorbar(im, cax=cax)
				#plt.subplots_adjust(wspace=0.0,bottom=0.05,top=0.95,left=0.05,right=0.8)#right=0.95)
				
				#plt.setp(ax1.get_xticklabels(), visible=False)	
				#plt.setp(ax2.get_xticklabels(), visible=False)
				
				#plt.setp(ax1.get_yticklabels(), visible=False)	
				#plt.setp(ax2.get_yticklabels(), visible=False)
	
				#ax1.set_title('True map  $\kappa_{center}=%.3f$'%(truen_im[psize,psize]))
				#ax2.set_title('Bias map  $\kappa_{center}=%.3f$'%(biasn_im[psize,psize]))
				
				##show()
				#savefig(filename)
				#close()
	print 'lostpeakcount',lostpeakcount
	print 'gainpeakcount',gainpeakcount
	
	######### plot stats for lost & new peaks ####
	#lostx_arr=lostx_arr.astype('int')
	#losty_arr=losty_arr.astype('int')
	#gainx_arr=gainx_arr.astype('int')
	#gainy_arr=gainy_arr.astype('int')
	#kappa_lost_bias=kmap_bias_noise[lostx_arr,losty_arr]
	#kappa_gain_bias=kmap_bias_noise[gainx_arr,gainy_arr]
	#kappa_lost_true=kmap_true_noise[lostx_arr,losty_arr]
	#kappa_gain_true=kmap_true_noise[gainx_arr,gainy_arr]
	
	#fig = plt.figure()
	#ax1 = fig.add_subplot(111)
	#kappa_arr=np.linspace(-0.04,0.14,10)
	#ax1.plot(kappa_arr,kappa_arr,'k-')	
	#ax1.scatter(kappa_lost_true,kappa_lost_bias,edgecolor='r',facecolor='r',s=3,label='lost peaks (%s)'%(len(kappa_lost_true)))
	#ax1.scatter(kappa_gain_true,kappa_gain_bias,edgecolor='b',facecolor='b',s=3,label='new peaks (%s)'%(len(kappa_gain_true)))	
	#ax1.set_ylim(-0.1,)
	#ax1.set_xlabel(r'$\kappa_{true}$')
	#ax1.set_ylabel(r'$\kappa_{bias}$')
	
	#ax2 = ax1.twinx()
	#ax2.hist(kappa_lost_true,range=(-0.04,0.14),align='mid',edgecolor='r',fill=False)	
	#ax2.hist(kappa_gain_true,range=(-0.04,0.14),align='mid',edgecolor='b',fill=False)
	#ax2.set_ylim(0,150)
	#ax2.set_ylabel(r'$N(\kappa_t)$')
	#leg=ax1.legend(loc=2, prop={'size':14})
	#leg.get_frame().set_visible(False)
	
	##show()
	#savefig(plot_dir2+'gainlostkappa.jpg')
	#close()
	
	
	
	######### plot individual high/low peaks #####
	#hipeak_loc=np.where((remain_true_noise_map>0.075)&(remain_true_noise_map<0.11))
	#lopeak_loc=np.where((remain_true_noise_map>-0.02)&(remain_true_noise_map<-0.001))
	
	#for i in range(40):#range(len(hipeak_loc[0])):
		#f = figure(figsize=(16,5))
		#ax1=f.add_subplot(151)
		#ax2=f.add_subplot(152)
		#ax3=f.add_subplot(153)
		#ax4=f.add_subplot(154)
		#ax5=f.add_subplot(155)
		
		#vmin,vmax=0.02,0.11
		#x,y=hipeak_loc[0][i],hipeak_loc[1][i]
		
		##vmin,vmax=-0.02,0.0
		##x,y=lopeak_loc[0][i],lopeak_loc[1][i]
		
		#true_im =kmap_true[x-10:x+11,y-10:y+11]
		#bias_im =kmap_bias[x-10:x+11,y-10:y+11]
		#truen_im=kmap_true_noise[x-10:x+11,y-10:y+11]
		#biasn_im=kmap_bias_noise[x-10:x+11,y-10:y+11]
		#GRF_im  =kmap_GRF[x-10:x+11,y-10:y+11]
		
		#ax1.imshow(true_im,norm=None,interpolation='nearest',vmin=vmin,vmax=vmax)		
		#ax2.imshow(bias_im,norm=None,interpolation='nearest',vmin=vmin,vmax=vmax)
		#ax3.imshow(truen_im,norm=None,interpolation='nearest',vmin=vmin,vmax=vmax)
		#ax4.imshow(biasn_im,norm=None,interpolation='nearest',vmin=vmin,vmax=vmax)
		#ax5.imshow(GRF_im,norm=None,interpolation='nearest',vmin=vmin,vmax=vmax)
		
		
		#ax1.set_title('Peak# %s Noisless Fiducial \n $\kappa=%s$'%(i,true_im[10,10]))
		#ax2.set_title('Noiseless Bias s%s \n $\kappa=%s$'%(s,bias_im[10,10]))
		#ax3.set_title('Noise Fiducial \n $\kappa=%s$'%(truen_im[10,10]))
		#ax4.set_title('Noise Bias s=%s \n $\kappa=%s$'%(s,biasn_im[10,10]))
		#ax5.set_title('GRF \n $\kappa=%s$'%(GRF_im[10,10]))
		
		#plt.subplots_adjust(wspace=0.0,left=0.05,right=0.95)
		#plt.setp(ax1.get_xticklabels(), visible=False)	
		#plt.setp(ax2.get_xticklabels(), visible=False)
		#plt.setp(ax3.get_xticklabels(), visible=False)	
		#plt.setp(ax4.get_xticklabels(), visible=False)
		#plt.setp(ax5.get_xticklabels(), visible=False)
		
		#plt.setp(ax1.get_yticklabels(), visible=False)	
		#plt.setp(ax2.get_yticklabels(), visible=False)
		#plt.setp(ax3.get_yticklabels(), visible=False)	
		#plt.setp(ax4.get_yticklabels(), visible=False)
		#plt.setp(ax5.get_yticklabels(), visible=False)
		
		#savefig(plot_dir2+'imshowpeak/hipeak%s.jpg'%(i))
		##show()
		##savefig(plot_dir2+'imshowpeak/lopeak%s.jpg'%(i))
		#close()
	
	######## plot ratio of bias/true map #######
	
	#remain_ratio=remain_bias_noise/remain_true_noise
	#scatter(remain_true_noise,remain_bias_noise)
	#ylabel('remain_bias_noise')
	#xlabel('remain_true_noise')
	#x=linspace(-0.03,0.15,100)
	##y=x*(1+(5*s-2)*x)
	#y=1.09*x+0.005
	#plot(x,y,'r-',linewidth=2,label='y=1.09*x+0.005')
	#legend()
	#title('s=%s'%(s))
	#savefig('singlemap_MB_remain_ratio_s%s.jpg'%(s))
	#close()
	
	#true_hist=returnpeakhist(kmap_true_noise,bins=10)
	#remain_hist=np.histogram(remain_true_noise, range=(low,high), bins=10)[0]
	#x=linspace(low,high,10)
	#plot(x,true_hist/float(sum(true_hist)),label='true hist',linewidth=2,drawstyle = 'steps-mid')
	#plot(x,remain_hist/float(sum(remain_hist)),label='remain hist',linewidth=2,drawstyle = 'steps-mid')
	#legend()
	#savefig('singlemap2_MB_remain_ratio_s%s.jpg'%(s))
	#close()
	
	######### plot individual high/low peaks #####
	#hipeak_loc=np.where((remain_true_noise_map>0.075)&(remain_true_noise_map<0.11))
	#lopeak_loc=np.where((remain_true_noise_map>-0.02)&(remain_true_noise_map<-0.001))
	
	#for i in range(40):#range(len(hipeak_loc[0])):
		#f = figure(figsize=(16,5))
		#ax1=f.add_subplot(151)
		#ax2=f.add_subplot(152)
		#ax3=f.add_subplot(153)
		#ax4=f.add_subplot(154)
		#ax5=f.add_subplot(155)
		
		#vmin,vmax=0.02,0.11
		#x,y=hipeak_loc[0][i],hipeak_loc[1][i]
		
		##vmin,vmax=-0.02,0.0
		##x,y=lopeak_loc[0][i],lopeak_loc[1][i]
		
		#true_im =kmap_true[x-10:x+11,y-10:y+11]
		#bias_im =kmap_bias[x-10:x+11,y-10:y+11]
		#truen_im=kmap_true_noise[x-10:x+11,y-10:y+11]
		#biasn_im=kmap_bias_noise[x-10:x+11,y-10:y+11]
		#GRF_im  =kmap_GRF[x-10:x+11,y-10:y+11]
		
		#ax1.imshow(true_im,norm=None,interpolation='nearest',vmin=vmin,vmax=vmax)		
		#ax2.imshow(bias_im,norm=None,interpolation='nearest',vmin=vmin,vmax=vmax)
		#ax3.imshow(truen_im,norm=None,interpolation='nearest',vmin=vmin,vmax=vmax)
		#ax4.imshow(biasn_im,norm=None,interpolation='nearest',vmin=vmin,vmax=vmax)
		#ax5.imshow(GRF_im,norm=None,interpolation='nearest',vmin=vmin,vmax=vmax)
		
		
		#ax1.set_title('Peak# %s Noisless Fiducial \n $\kappa=%s$'%(i,true_im[10,10]))
		#ax2.set_title('Noiseless Bias s%s \n $\kappa=%s$'%(s,bias_im[10,10]))
		#ax3.set_title('Noise Fiducial \n $\kappa=%s$'%(truen_im[10,10]))
		#ax4.set_title('Noise Bias s=%s \n $\kappa=%s$'%(s,biasn_im[10,10]))
		#ax5.set_title('GRF \n $\kappa=%s$'%(GRF_im[10,10]))
		
		#plt.subplots_adjust(wspace=0.0,left=0.05,right=0.95)
		#plt.setp(ax1.get_xticklabels(), visible=False)	
		#plt.setp(ax2.get_xticklabels(), visible=False)
		#plt.setp(ax3.get_xticklabels(), visible=False)	
		#plt.setp(ax4.get_xticklabels(), visible=False)
		#plt.setp(ax5.get_xticklabels(), visible=False)
		
		#plt.setp(ax1.get_yticklabels(), visible=False)	
		#plt.setp(ax2.get_yticklabels(), visible=False)
		#plt.setp(ax3.get_yticklabels(), visible=False)	
		#plt.setp(ax4.get_yticklabels(), visible=False)
		#plt.setp(ax5.get_yticklabels(), visible=False)
		
		#savefig(plot_dir2+'imshowpeak/hipeak%s.jpg'%(i))
		##show()
		##savefig(plot_dir2+'imshowpeak/lopeak%s.jpg'%(i))
		#close()
	
	
	######### plot remained peak distribution #######
	#f = figure(figsize=(8,6))
	#hist_remain_true=np.histogram(remain_true_noise, range=(low,high), bins=100)[0]
	#hist_remain_bias=np.histogram(remain_bias_noise, range=(low,high), bins=100)[0]
	#hist_remain_true=np.array(hist_remain_true,dtype=float)
	#hist_remain_bias=np.array(hist_remain_bias,dtype=float)
	
	#ax1=f.add_subplot(211)
	#ax2=f.add_subplot(212,sharex=ax1)
	#ax1.plot(hist_bins,hist_remain_true,'*',linewidth=1.5,label='Fiducial')
	#ax1.plot(hist_bins,hist_remain_bias,'d',linewidth=1.5,label='Bias s=%s'%(s))
	
	#ax2.plot(hist_bins,zeros(bins))
	#ax2.plot(hist_bins,hist_remain_bias/hist_remain_true-1,'o',linewidth=1.5)
	
	#leg1=ax1.legend(loc=0, ncol=1, labelspacing=.2, prop={'size':12})
	#leg1.get_frame().set_visible(False)
	#ax1.set_xlabel('convergence $\kappa$')
	
	#ax1.set_xlim(-0.02,0.15)
	#ax2.set_ylim(-2,5)
	##ax1.set_yscale('log',basey=5)
	#plt.subplots_adjust(hspace=0.0)
	#ax1.set_title('Remained peak distribution Ngal=%s'%(ngal))
	#show()	
	
	
	############### plot GRF peak histogram ######
	#peak_noise=returnpeakhist(kmap_GRF,1000)
	#f = figure(figsize=(8,6))
	#ax1=f.add_subplot(111)
	#ax1.plot(hist_bins,peak_noise,'*',linewidth=1.5)
	#ax1.set_xlabel('convergence $\kappa$')
	#ax1.set_xlim(-0.005,0.01)
	#ax1.set_yscale('log',basey=5)
	#title('Noise map peak distribution Ngal=%s'%(ngal))
	#show()
	
	
	
	############# plot peak distribution on single map ###
	#peaks_true=returnpeakhist(kmap_true,bins)
	#peaks_true_noise=returnpeakhist(kmap_true_noise,bins)
	
	#peaks_bias=np.ndarray(shape=(len(s_arr3),bins))
	#peaks_bias_noise=np.ndarray(shape=(len(s_arr3),bins))
	#i=0
	#for s in s_arr3:		
		#kmap_bias = biasmooth(kmap, sigma, s)
		#kmap_bias_noise = biasnoisesmooth(kmap,noise,sigma,s)
		#peaks_bias[i]=returnpeakhist(kmap_bias,bins)
		#peaks_bias_noise[i]=returnpeakhist(kmap_bias_noise,bins)
		#i+=1
	
	
	#f = figure(figsize=(12,6))
	#gs = gridspec.GridSpec(2,2,height_ratios=[3,1])
	#ax1=f.add_subplot(gs[0])
	#ax2=f.add_subplot(gs[1])
	#ax3=f.add_subplot(gs[2],sharex=ax1)
	#ax4=f.add_subplot(gs[3],sharex=ax2)
	
	#ax1.plot(hist_bins,peaks_true,linewidth=1.5, label='Fiducial (Noiseless)')
	#ax2.plot(hist_bins,peaks_true_noise,linewidth=1.5, label='Fiducial (Noise)')
	#ax3.plot(hist_bins,zeros(bins),linewidth=1.5)
	#ax4.plot(hist_bins,zeros(bins),linewidth=1.5)
	#for n in range(len(s_arr3)):
		#ax1.plot(hist_bins,peaks_bias[n],linewidth=1.5, label='Bias s=%s'%(s_arr3[n]))
		#ax2.plot(hist_bins,peaks_bias_noise[n],linewidth=1.5, label='Bias s=%s'%(s_arr3[n]))	
		#ax3.plot(hist_bins,peaks_bias[n]/peaks_true-1,linewidth=1.5, label='Bias s=%s'%(s_arr3[n]))
		#ax4.plot(hist_bins,peaks_bias_noise[n]/peaks_true_noise-1,linewidth=1.5, label='Bias s=%s'%(s_arr3[n]))
	#leg1=ax1.legend(loc=0, ncol=1, labelspacing=.2, prop={'size':14})
	#leg1.get_frame().set_visible(False)
	#leg2=ax2.legend(loc=0, ncol=1, labelspacing=.2, prop={'size':14})
	#leg2.get_frame().set_visible(False)
	
	
	#ax3.set_xlabel('convergence $\kappa$')
	#ax1.set_xlabel('convergence $\kappa$')
	#ax1.set_ylabel('peak counts $N$($\kappa$)')
	#ax3.set_ylabel(r'${\Delta}N/N$')
	
	#ax1.set_xlim(-0.02,0.1)
	#ax2.set_xlim(-0.02,0.1)
	
	#plt.setp(ax1.get_xticklabels(), visible=False)	
	#plt.setp(ax2.get_xticklabels(), visible=False)
	#plt.subplots_adjust(hspace=0.0)
	
	#show()

if seeded100noises:
	fisher_mat=np.genfromtxt('/Users/jia/magbias/fisher_mat.ls')
	f = figure(figsize=(12,4))
	for i in range(3):
		ax = subplot(1,3,i+1) #subplots 2,4,6
		fits = fisher_mat[:,i]
		hist1, edges = histogram(fits, bins=10)
		midpoints = edges[:-1] + 0.5 * (edges[2] - edges[1])
		ax.plot(midpoints, hist1, linewidth=2,drawstyle = 'steps-mid')
		title(label_str[::-1][i])
		#leg=ax.legend(loc=0, ncol=1, labelspacing=.2, prop={'size':10})
		#leg.get_frame().set_visible(False)			
	
	savefig(plot_dir+'seeded100noises.jpg')
	show()	

if noise_jia_vs_xy:
	cov_mat_file='covmat_g6.976_ngal15_bins15.ls'
	
	true_jia = '/Users/jia/magbias/fit_noise_xy/'
	true_xy = '/Users/jia/magbias/fit_noise_xy_usenoisexy60s/'
	
	true_jia = '/Users/jia/magbias/fit_noise_xy/'
	true_xy = '/Users/jia/magbias/fit_noise_xy_usenoisexy60s/'
	
	cov_mat_jia=np.genfromtxt(true_jia+cov_mat_file)
	cov_mat_xy=np.genfromtxt(true_xy+cov_mat_file)
	
	contr_arr = np.array([236.657000, 238.526000, 236.721000, 241.410000, 241.211000, 245.687000, 241.950000, 240.693000, 242.570000, 244.987000, 243.602000, 240.994000, 239.980000, 239.609000, 188.378000])
	lo_si_arr = np.array([225.992000, 235.599000, 236.670000, 244.670000, 245.980000, 251.634000,248.491000, 247.972000, 250.320000, 253.431000, 251.998000, 247.408000,245.414000, 237.749000, 167.363000])
	lo_w_arr = np.array([242.059000, 238.986000, 234.175000, 239.478000, 237.405000, 242.610000,238.618000, 237.780000, 238.333000, 242.242000, 241.476000, 239.155000,241.063000, 243.221000, 196.402000])
	lo_Om_arr = np.array([224.225000, 236.739000, 238.602000, 246.758000, 247.544000, 255.071000,250.867000, 250.408000, 252.222000, 254.013000, 252.042000, 245.739000,242.993000, 234.564000, 164.070000])
	hi_si_arr = np.array([248.993000, 243.004000, 235.649000, 239.312000, 237.058000, 240.182000,234.147000, 232.824000, 233.505000, 235.946000, 235.933000, 232.862000,235.102000, 240.018000, 210.098000])
	hi_w_arr = np.array([231.413000, 237.891000, 237.694000, 243.794000, 243.989000, 250.932000,244.824000, 243.887000, 246.203000, 249.014000, 246.310000, 241.423000,239.337000, 234.988000, 180.593000])
	hi_Om_arr = np.array([249.463000, 240.550000, 233.385000, 236.433000, 234.177000, 237.878000,232.876000, 231.764000, 233.245000, 236.034000, 235.728000, 234.135000,236.809000, 244.649000, 212.716000])
	
	labels=('lo_Om','hi_w','lo_si','contr','hi_si','lo_w','hi_Om')
	peakfiles=('peak_matrix_true_1000m-512b240_Om0.230_Ol0.770_w-1.000_ns0.960_si0.798_g6.976_ngal15_bins15.ls',
'peak_matrix_true_1000m-512b240_Om0.260_Ol0.740_w-0.800_ns0.960_si0.798_g6.976_ngal15_bins15.ls',
'peak_matrix_true_1000m-512b240_Om0.260_Ol0.740_w-1.000_ns0.960_si0.750_g6.976_ngal15_bins15.ls',
'peak_matrix_true_1000m-512b240_Om0.260_Ol0.740_w-1.000_ns0.960_si0.798_g6.976_ngal15_bins15.ls',
'peak_matrix_true_1000m-512b240_Om0.260_Ol0.740_w-1.000_ns0.960_si0.850_g6.976_ngal15_bins15.ls',
'peak_matrix_true_1000m-512b240_Om0.260_Ol0.740_w-1.200_ns0.960_si0.798_g6.976_ngal15_bins15.ls',
'peak_matrix_true_1000m-512b240_Om0.290_Ol0.710_w-1.000_ns0.960_si0.798_g6.976_ngal15_bins15.ls')
	
	true_mat_jia=np.ndarray(shape=(7,1000,15))
	true_mat_xy=np.ndarray(shape=(7,1000,15))
	k=0
	for peak_mat in peakfiles:
		true_mat_jia[k]=np.genfromtxt(true_jia+peakfiles[k])
		true_mat_xy[k]=np.genfromtxt(true_xy+peakfiles[k])
		k+=1
	
	fidu_mat_jia=average(true_mat_jia,axis=1)
	fidu_mat_xy=average(true_mat_xy,axis=1)
	# that calculated by xiuyuan
	fidu_mat_xiuyuan = [lo_Om_arr,hi_w_arr,lo_si_arr,contr_arr,hi_si_arr,lo_w_arr,hi_Om_arr]
	
	his = [4,1,6]
	los = [2,5,0]
	contr = [3]
	
	def fishererr(fidu_mat,cov_mat):
		hi_arr,lo_arr,contr_arr = fidu_mat[his],fidu_mat[los],fidu_mat[contr]
		
		cov_inv = np.mat(cov_mat).I
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
		err = sqrt(mat(Fisher).I)[[0,1,2],[0,1,2]]/sqrt(20000.0/12.0)
		
		return Fisher,err	
	
	
	
	
## replace one dN with dN from jia(xy), use everything else of xy(jia)
	fidu_mat_jia_xydw = np.ndarray(shape=(7,15))
	fidu_mat_jia_xydw[[1,5]]=fidu_mat_xy[[1,5]]
	fidu_mat_jia_xydw[[0,2,3,4,6]]=fidu_mat_jia[[0,2,3,4,6]]
	
	fidu_mat_jia_xyds = np.ndarray(shape=(7,15))
	fidu_mat_jia_xyds[[2,4]]=fidu_mat_xy[[2,4]]
	fidu_mat_jia_xyds[[0,1,3,5,6]]=fidu_mat_jia[[0,1,3,5,6]]
	
	fidu_mat_jia_xydom = np.ndarray(shape=(7,15))
	fidu_mat_jia_xydom[[0,6]]=fidu_mat_xy[[0,6]]
	fidu_mat_jia_xydom[[1,2,3,4,5]]=fidu_mat_jia[[1,2,3,4,5]]
	
	
	fidu_mat_xy_jiadw = np.ndarray(shape=(7,15))
	fidu_mat_xy_jiadw[[1,5]]=fidu_mat_jia[[1,5]]
	fidu_mat_xy_jiadw[[0,2,3,4,6]]=fidu_mat_xy[[0,2,3,4,6]]
	
	fidu_mat_xy_jiads = np.ndarray(shape=(7,15))
	fidu_mat_xy_jiads[[2,4]]=fidu_mat_jia[[2,4]]
	fidu_mat_xy_jiads[[0,1,3,5,6]]=fidu_mat_xy[[0,1,3,5,6]]
	
	fidu_mat_xy_jiadom = np.ndarray(shape=(7,15))
	fidu_mat_xy_jiadom[[0,6]]=fidu_mat_jia[[0,6]]
	fidu_mat_xy_jiadom[[1,2,3,4,5]]=fidu_mat_xy[[1,2,3,4,5]]
	
	
	#############################################
	#######plot fluctuation of dN ###############
	def errorbars_on_peak_distr(a_true_mat):
		a_true_mat = a_true_mat.T
		err_arr = np.ndarray(shape=(a_true_mat.shape[0]))
		i=0
		for ibin in a_true_mat:
			#print stats.bayes_mvs(ibin,0.68)
			err_arr[i]=0.5*(stats.bayes_mvs(ibin,0.68)[2][1][0]+stats.bayes_mvs(ibin,0.68)[2][1][1])
			i+=1
		return err_arr
	
	#f = figure(figsize=(12,9))	
	#for i in range(7):
		#ax=subplot(3,3,i+1)
		
		#yerr_jia = errorbars_on_peak_distr(true_mat_jia[i])
		#yerr_xy = errorbars_on_peak_distr(true_mat_xy[i])
		
		#ax.errorbar(np.arange(15),fidu_mat_jia[i],yerr=yerr_jia, label='Jia (%s)'%(sum(fidu_mat_jia[i])),linewidth=2, drawstyle = 'steps-mid')
		
		#ax.errorbar(np.arange(15),fidu_mat_xy[i], yerr=yerr_xy, label='Xiuyuan (%s)'%(sum(fidu_mat_xiuyuan[i])),linewidth=2, drawstyle = 'steps-mid')
		
		#leg=ax.legend(loc=0, ncol=1, labelspacing=.2, prop={'size':10})
		#leg.get_frame().set_visible(False)
		
		#title(labels[i])
	#savefig(plot_dir+'peak_distr_fluctuations.jpg')
	#close()
	
	f = figure(figsize=(12,6))
	j=1
	for i in (0,1,2,4,5,6):
		ax=subplot(2,3,j)
		
		yerr_jia = errorbars_on_peak_distr(true_mat_jia[i])
		yerr_xy = errorbars_on_peak_distr(true_mat_xy[i])
		
		ax.errorbar(np.arange(15),fidu_mat_jia[i]-fidu_mat_jia[3],yerr=yerr_jia, label='Jia (%s)'%(sum(fidu_mat_jia[i])),linewidth=2, drawstyle = 'steps-mid')
		
		ax.errorbar(np.arange(15),fidu_mat_xy[i]-fidu_mat_xy[3], yerr=yerr_xy, label='Xiuyuan (%s)'%(sum(fidu_mat_xiuyuan[i])),linewidth=2, drawstyle = 'steps-mid')
		
		leg=ax.legend(loc=0, ncol=1, labelspacing=.2, prop={'size':10})
		leg.get_frame().set_visible(False)
		
		title(labels[i])
		j+=1
	savefig(plot_dir+'peak_distr_dN_fluctuations.jpg')
	show()
	
	
	#############################################
	#i=0
	#f = figure(figsize=(12,9))
	#for peak_mat in peakfiles:		
		#true_mat_jia=np.genfromtxt(true_jia+peak_mat)
		#true_mat_xy =np.genfromtxt(true_xy+peak_mat)
		
		#ax=subplot(3,3,i+1)

		#ax.plot(average(true_mat_jia,axis=0),label='Jia noise+count (%s)'%(0.001*sum(true_mat_jia)),drawstyle = 'steps-mid',linewidth=2)
		
		#ax.plot(average(true_mat_xy,axis=0),label='XY noise, Jia count(%s)'%(0.001*sum(true_mat_xy)),drawstyle = 'steps-mid',linewidth=2)
		
		#ax.plot(fidu_mat_xiuyuan[i],label='Xiuyuan noise+count(%s)'%(sum(fidu_mat_xiuyuan[i])),drawstyle = 'steps-mid',linewidth=2)
		
		#leg=ax.legend(loc=0, ncol=1, labelspacing=.2, prop={'size':10})
		#leg.get_frame().set_visible(False)
		
		#xlabel('bins')
		#title = labels[i]
		
		#i+=1
	#savefig(plot_dir+'peak_distr_xybins_noise2.jpg')
	#close()
	########################################################
	#f = figure(figsize=(16,4))
	#subplot(1,3,1)
	#title('Jia Covariance Matrix')
	#imshow(cov_mat_jia, interpolation='nearest', aspect='auto', origin='lower')
	#colorbar()
	#subplot(1,3,2)
	#title('XY Covariance Matrix')
	#imshow(cov_mat_xy, interpolation='nearest', aspect='auto', origin='lower')
	#colorbar()
	#subplot(1,3,3)
	#title('Frac. Diff')
	#imshow(cov_mat_jia/cov_mat_xy-1, interpolation='nearest', aspect='auto', origin='lower')
	#colorbar()
	#savefig(plot_dir+'jia_xy_cov_diff.jpg')
	#close()
	##########################################################
	
	##################################
	#f = figure(figsize=(5,10))
	#subplot(3,1,1)
	#title('Jia Fisher Matrix')
	#imshow(F_jia, interpolation='nearest', aspect='auto', origin='lower')
	#colorbar()
	#subplot(3,1,2)
	#title('XY Fisher Matrix')
	#imshow(F_xy, interpolation='nearest', aspect='auto', origin='lower')
	#colorbar()
	#subplot(3,1,3)
	#title('fra diff. XY from Jia')
	#imshow(F_xy/F_jia-1, interpolation='nearest', aspect='auto', origin='lower')
	#colorbar()
	#savefig(plot_dir+'jia_xy_fisher_diff.jpg')
	#close()	
	
	##################################
	#xx=.5*(fidu_mat_jia[his]-fidu_mat_jia[los])
	#yy=.5*(fidu_mat_xy[his]-fidu_mat_xy[los])
	
	#f = figure(figsize=(14,8))
	#for x in range(3):
		#ax = subplot(2,3,x+1)
		#ax.plot(xx[x],linewidth=2, drawstyle = 'steps-mid',label='Jia')
		#ax.plot(yy[x],linewidth=2, drawstyle = 'steps-mid',label='Xiuyua')
		#title('del_N '+label_str[::-1][x])	
		#leg=ax.legend(loc=0, ncol=1, labelspacing=.2, prop={'size':10})
		#leg.get_frame().set_visible(False)
		##xlabel('ith bin')
		#ylabel('del_N')
		
		#ax2=subplot(2,3,x+4)
		#ax2.plot(xx[x]/yy[x]-1,linewidth=2, drawstyle = 'steps-mid')
		#title('Frac Diff (Jia/XY-1) '+label_str[::-1][x])
		#xlabel('ith bin')
		#ylabel('Frac Diff')
	#savefig(plot_dir+'jia_xy_del_N.jpg')

	##################################
	## test suggestion # 2 by zoltan
	
	#print 'Randomly draw 500 maps from 1000 maps to calculate errors'
	##print 'Jia'
	
	#samplenum=500
	#mapnum=900
	#random_errs_jia=np.ndarray(shape=(samplenum,3))
	#random_errs_xy=np.ndarray(shape=(samplenum,3))
	#for x in range(samplenum):
		#maps=np.array(rand(mapnum)*1000,dtype=int)
		#partial_maps=true_mat_jia[:,maps,:]
		#fidu_matx=average(partial_maps,axis=1)
		#cov_matx=np.cov(partial_maps[3],rowvar=0)
		#fisherx,errx=fishererr(fidu_matx,cov_matx)
		#random_errs_jia[x]=errx
		##print x, errx
		
	##print 'Xiuyuan'
	#for x in range(samplenum):
		#maps=np.array(rand(mapnum)*1000,dtype=int)
		#partial_maps=true_mat_xy[:,maps,:]
		#fidu_matx=average(partial_maps,axis=1)
		#cov_matx=np.cov(partial_maps[3],rowvar=0)
		#fisherx,errx=fishererr(fidu_matx,cov_matx)
		#random_errs_xy[x]=errx
		##print x, errx
	
	#f = figure(figsize=(5,8))
	
	#for x in range(3):
		#a1,b1=histogram(random_errs_jia[:,x],bins=20)
		#a2,b2=histogram(random_errs_xy[:,x],bins=b1)
		#xbin=b1[:-1]+0.5*(b1[1]-b1[0])
		
		#ax = subplot(3,1,x+1)
		#ax.plot(xbin,a1,linewidth=2, drawstyle = 'steps-mid',label='Jia')
		#ax.plot(xbin,a2,linewidth=2, drawstyle = 'steps-mid',label='Xiuyuan')
		##title('error - '+label_str[::-1][x])	
		#leg=ax.legend(loc=0, ncol=1, labelspacing=.2, prop={'size':10})
		#leg.get_frame().set_visible(False)
		#xlabel('err '+label_str[::-1][x])
		##ylabel('del_N')
	#subplot(3,1,1)
	#title('randomly draw %s maps (%s times)'%(mapnum,samplenum))
	#savefig(plot_dir+'jia_xy_%smaps.jpg'%(mapnum))
	#show()
	

if Fisher_vs_Chisq_test:
	from scipy import interpolate
	cov_dir = '/Users/jia/magbias/cov_dir/'
	hi_si,lo_si,hi_w,lo_w,hi_Om,lo_Om,contr=4,2,1,5,6,0,3
	devhi = np.array([0.052,0.2,0.03])
	devlo = np.array([0.048,0.2,0.03])
	
	def fisher_err_and_chisqpair (fidu_mat,cov_inv):
		
		hi_arr =fidu_mat[[hi_si,hi_w,hi_Om]]
		lo_arr = fidu_mat[[lo_si,lo_w,lo_Om]]
		contr_arr = fidu_mat[contr]
		
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

		err = sqrt(mat(Fisher)**-1)[[0,1,2],[0,1,2]]/sqrt(20000.0/12.0)
		
		####### calculating chi^2 pairs
		
		interpOm = lambda ibin: interpolate.UnivariateSpline(Om,fidu_mat[[lo_Om,contr,hi_Om], ibin], k=2, s=0)
		interpw  = lambda ibin: interpolate.UnivariateSpline(w, fidu_mat[[lo_w ,contr,hi_w ], ibin], k=2, s=0)
		interpsi = lambda ibin: interpolate.UnivariateSpline(si,fidu_mat[[lo_si,contr,hi_si], ibin], k=2, s=0)

		Ni = lambda ibin, iOm, iw, isi: interpOm(ibin)(iOm)+interpw(ibin)(iw)+ interpsi(ibin)(isi)-2.0*fidu_mat[contr,ibin]

		Ni_arr = lambda iOm, iw, isi: np.array(map(Ni, np.arange(bins),array([iOm]*bins),array([iw]*bins),array([isi]*bins)))
		
		def chisq (initguess, ibmap):
			iOm, iw, isi = initguess
			del_Ni = np.mat(Ni_arr(iOm, iw, isi) - ibmap)
			chisquare = float64(del_Ni * cov_inv * del_Ni.T)
		#	print chisquare, initguess
			return chisquare
		
		chisqpair=np.ndarray(shape=(6,))
		i=0
		for icosmo in fidu_mat[[hi_si,lo_si,hi_w,lo_w,hi_Om,lo_Om]]:
			chisqpair[i] = chisq((0.26,-1.0,0.798),icosmo)
			i += 1
		
		return err, chisqpair
	
	fisher_noiseless = np.ndarray(shape=(10,3))
	fisher_ngal15    = np.ndarray(shape=(10,3))
	fisher5_noiseless= np.ndarray(shape=(10,3))
	fisher5_ngal15   = np.ndarray(shape=(10,3))
	
	chisqpairs_noiseless  = np.ndarray(shape=(10,6))
	chisqpairs_ngal15     = np.ndarray(shape=(10,6))
	chisqpairs5_noiseless = np.ndarray(shape=(10,6))
	chisqpairs5_ngal15    = np.ndarray(shape=(10,6))
	
	j = 0
	for bins in bins_arr:
		
		cov_inv_noiseless = np.genfromtxt(cov_dir+'Cov_inv_noiseless_bins%s.ls'%(bins))
		fidu_mat_noiseless = np.genfromtxt(cov_dir+'fidu_mat_noiseless_bins%s.ls'%(bins))
		cov_inv_ngal15 = np.genfromtxt(cov_dir+'Cov_inv_ngal15_bins%s.ls'%(bins))
		fidu_mat_ngal15 = np.genfromtxt(cov_dir+'fidu_mat_ngal15_bins%s.ls'%(bins))
		
		fisher_noiseless [j],chisqpairs_noiseless  [j]= fisher_err_and_chisqpair (fidu_mat_noiseless,cov_inv_noiseless)
		fisher_ngal15    [j],chisqpairs_ngal15     [j]= fisher_err_and_chisqpair (fidu_mat_ngal15,cov_inv_ngal15)
		fisher5_noiseless[j],chisqpairs5_noiseless [j]= fisher_err_and_chisqpair (fidu_mat_noiseless*5,cov_inv_noiseless)
		fisher5_ngal15   [j],chisqpairs5_ngal15    [j]= fisher_err_and_chisqpair (fidu_mat_ngal15*5,cov_inv_ngal15)
		j+=1
		
	f = figure(figsize=(14,8))
	for x in range(3):
		ax = subplot(2,3,x+1)
		ax.plot(bins_arr,fisher_noiseless[:,x],label='Noiseless')
		ax.plot(bins_arr,fisher_ngal15[:,x],label='Ngal15')
		ax.plot(bins_arr,fisher5_noiseless[:,x],label='Noiseless dN*5')
		ax.plot(bins_arr,fisher5_ngal15[:,x],label='Ngal15 dN*5')
		title(label_str[::-1][x])
		
		leg=ax.legend(loc=0, ncol=1, labelspacing=.2, prop={'size':10})
		leg.get_frame().set_visible(False)
		xlabel('bins')
		ylabel('err '+label_str[::-1][x])
		
		
		
		ax2 = subplot(2,3,x+4)
		ax2.plot(bins_arr,chisqpairs_noiseless[:,x],label='Noiseless')
		ax2.plot(bins_arr,chisqpairs_ngal15[:,x],label='Ngal15')
		ax2.plot(bins_arr,chisqpairs5_noiseless[:,x],label='Noiseless dN*5')
		ax2.plot(bins_arr,chisqpairs5_ngal15[:,x],label='Ngal15 dN*5')
		title(label_str[::-1][x])
		
		leg2=ax2.legend(loc=0, ncol=1, labelspacing=.2, prop={'size':10})
		leg2.get_frame().set_visible(False)
		xlabel('bins')
		ylabel('err '+label_str[::-1][x])
		
	show()

	

if comp_xiuyuan_distribution:
	contr_arr = np.array([236.657000, 238.526000, 236.721000, 241.410000, 241.211000, 245.687000, 241.950000, 240.693000, 242.570000, 244.987000, 243.602000, 240.994000, 239.980000, 239.609000, 188.378000])
	lo_si_arr = np.array([225.992000, 235.599000, 236.670000, 244.670000, 245.980000, 251.634000,248.491000, 247.972000, 250.320000, 253.431000, 251.998000, 247.408000,245.414000, 237.749000, 167.363000])
	lo_w_arr = np.array([242.059000, 238.986000, 234.175000, 239.478000, 237.405000, 242.610000,238.618000, 237.780000, 238.333000, 242.242000, 241.476000, 239.155000,241.063000, 243.221000, 196.402000])
	lo_Om_arr = np.array([224.225000, 236.739000, 238.602000, 246.758000, 247.544000, 255.071000,250.867000, 250.408000, 252.222000, 254.013000, 252.042000, 245.739000,242.993000, 234.564000, 164.070000])
	hi_si_arr = np.array([248.993000, 243.004000, 235.649000, 239.312000, 237.058000, 240.182000,234.147000, 232.824000, 233.505000, 235.946000, 235.933000, 232.862000,235.102000, 240.018000, 210.098000])
	hi_w_arr = np.array([231.413000, 237.891000, 237.694000, 243.794000, 243.989000, 250.932000,244.824000, 243.887000, 246.203000, 249.014000, 246.310000, 241.423000,239.337000, 234.988000, 180.593000])
	hi_Om_arr = np.array([249.463000, 240.550000, 233.385000, 236.433000, 234.177000, 237.878000,232.876000, 231.764000, 233.245000, 236.034000, 235.728000, 234.135000,236.809000, 244.649000, 212.716000])
	
	
	fidu_mat=np.genfromtxt('/Users/jia/magbias/cov_dir/fidu_mat_ngal15_bins15.ls')
	#fidu_mat_noiseless=np.genfromtxt('/Users/jia/magbias/cov_dir/fidu_mat-bins15.ls')
	labels = ('hi_si','lo_si','hi_w','lo_w','hi_Om','lo_Om','contr')
	hi_si,lo_si,hi_w,lo_w,hi_Om,lo_Om,contr=4,2,1,5,6,0,3
	fidu_mat_jia = fidu_mat[[4,2,1,5,6,0,3]]
	#fidu_mat_jia_noiseless = fidu_mat_noiseless[[4,2,1,5,6,0,3]]
	fidu_mat_xiuyuan = [hi_si_arr,lo_si_arr,hi_w_arr,lo_w_arr,hi_Om_arr,lo_Om_arr,contr_arr]
	f = figure(figsize=(12,9))
	for i in range(7):
		ax=subplot(3,3,i+1)
		#plot(fidu_mat_jia_noiseless[i],label='Jia noiseless (%s)'%(sum(fidu_mat_jia_noiseless[i])),linewidth=2, drawstyle = 'steps-mid')
		
		ax.plot(fidu_mat_jia[i],label='Jia (%s)'%(sum(fidu_mat_jia[i])),linewidth=2, drawstyle = 'steps-mid')
		
		ax.plot(fidu_mat_xiuyuan[i],label='Xiuyuan (%s)'%(sum(fidu_mat_xiuyuan[i])),linewidth=2, drawstyle = 'steps-mid')
		leg=ax.legend(loc=0, ncol=1, labelspacing=.2, prop={'size':10})
		leg.get_frame().set_visible(False)
		title(labels[i])
	savefig(plot_dir+'peak_distr_xybins_noise.jpg')
	close()
	
if all_errs:
	xiuyuan = [0.0078,0.036,0.0057]
	skyscale = sqrt(20000./12)
	err_dir='/Users/jia/magbias/err/'
	#err_xy_g6_fw=np.genfromtxt(err_dir+'errfw_bins_xy_noiseless.ls')
	#err_xy_g6_bw=np.genfromtxt(err_dir+'errbw_bins_xy_noiseless.ls')
	#err_xy_g6=0.5*(err_xy_g6_fw+err_xy_g6_bw)
	err_xy_g6_ngal15=np.genfromtxt(err_dir+'err_g6.976_ngal15_bins_xy.ls')
	err_g6=np.ndarray(shape=(10,3))
	err_g9=np.ndarray(shape=(10,3))
	err_g6_ngal15=np.ndarray(shape=(10,3))
	err_g6_ngal30=np.ndarray(shape=(10,3))
	i=0
	for bins in bins_arr:
		err_g6[i]=np.genfromtxt(err_dir+'err_g6.976_bins%s.ls'%(bins))
		err_g9[i]=np.genfromtxt(err_dir+'err_g9.865_bins%s.ls'%(bins))
		err_g6_ngal15[i]=np.genfromtxt(err_dir+'err_g6.976_ngal15_bins%s.ls'%(bins))
		err_g6_ngal30[i]=np.genfromtxt(err_dir+'err_g6.976_ngal30_bins%s.ls'%(bins))
		i+=1
	f = figure(figsize=(12,4))
	for x in range(3):
		ax = subplot(1,3,x+1)
		ax.plot(bins_arr,err_g6[:,x],label='noisless 6.976 px',linewidth=1.5)
		ax.plot(bins_arr,err_g9[:,x],label='noisless 9.865 px',linewidth=1.5)
		
		ax.plot(bins_arr,err_g6_ngal15[:,x],label='ngal=15 6.976 pix',linewidth=1.5)
		ax.plot(bins_arr,err_g6_ngal30[:,x],label='ngal=30 6.976 pix',linewidth=1.5)
		ax.plot(bins_arr,zeros(len(bins_arr))+err_xy_g6_ngal15[x],label='ngal=15 bins=xiuyuan',linewidth=1.5)
		ax.plot(bins_arr,zeros(len(bins_arr))+xiuyuan[x],label='ngal=15 Yang 2011',linewidth=1.5)		
		
		leg=ax.legend(loc=0, ncol=1, labelspacing=.2, prop={'size':10})
		leg.get_frame().set_visible(False)
		xlabel('bins')
		ylabel('err '+label_str[::-1][x])

		
	savefig(plot_dir+'err_all.jpg')
	show()
	

if fisher_err:
	xiuyuan = [0.0078,0.036,0.0057]
	errfw = np.ndarray(shape=(10,3))
	errbw = np.ndarray(shape=(10,3))
	errfit = np.ndarray(shape=(10,3))
	errfit_bayes_mvs = np.ndarray(shape=(10,3))
	skyscale = sqrt(20000./12)
	i = 0
	for bins in bins_arr:
		errfw_file = fit_dir+'errfw_bins%s.ls'%(bins)
		errbw_file = fit_dir+'errbw_bins%s.ls'%(bins)
		errfw[i] = genfromtxt(errfw_file)
		errbw[i] = genfromtxt(errbw_file)		
		# width of my fits
		for x in range(3):
			fits = Tfit[i,:,x+1]
			hist1, edges = histogram(fits, bins=20)
			step = (edges[2] - edges[1])
			midpoints = edges[:-1] + 0.5 * step
			ind = np.argsort(hist1)[::-1]
			j = 1
			tot = sum(hist1)
			while sum(hist1[ind[:j]]) < 0.683*tot:
				j += 1
			width = 0.5*(np.amax(midpoints[ind[:j+1]]) - np.amin(midpoints[ind[:j+1]]))
			errfit[i,x] = width
			errfit_bayes_mvs[i,x] = 0.5*(stats.bayes_mvs(fits,0.68)[2][1][1]+stats.bayes_mvs(fits,0.68)[2][1][0])
		i += 1
	print 'errfit',errfit
	print 'errfit_bayes_mvs',errfit_bayes_mvs
	f = figure(figsize=(12,4))
	for x in range(3):
		ax = subplot(1,3,x+1)
		ax.plot(bins_arr,errfw[:,x]/skyscale,label='Fisher fw',linewidth=1.5)
		ax.plot(bins_arr,errbw[:,x]/skyscale,label='Fisher bw',linewidth=1.5)
		ax.plot(bins_arr,0.5*(errfw[:,x]+errbw[:,x])/skyscale,label='Fisher avg',linewidth=1.5)
		ax.plot(bins_arr,errfit[:,x]/skyscale,label='Fit 1D width',linewidth=1.5)
		ax.plot(bins_arr,errfit_bayes_mvs[:,x]/skyscale,label='Fit 1D width (baysian)',linewidth=1.5)
		ax.plot(bins_arr,zeros(len(bins_arr))+xiuyuan[x],label='Yang et el. 2011',linewidth=1.5)
		leg=ax.legend(loc=0, ncol=1, labelspacing=.2, prop={'size':10})
		leg.get_frame().set_visible(False)
		xlabel('bins')
		ylabel('err '+label_str[x])
	savefig(plot_dir+'fisher.jpg')
	#show()
	close()


if param_width:	
	true_width = np.ndarray(shape=(10,3))
	bias_width = np.ndarray(shape=(5,10,3))
	for i in range(3):
		ibin = 0
		for bins in bins_arr:	
		
			fits = Tfit[ibin,:,i+1]
			hist1, edges = histogram(fits, bins=20)
			step = (edges[2] - edges[1])
			midpoints = edges[:-1] + 0.5 * step
			ind = np.argsort(hist1)[::-1]
			j = 1
			tot = sum(hist1)
			while sum(hist1[ind[:j]]) < 0.683*tot:
				j += 1
			width = 0.5*(np.amax(midpoints[ind[:j+1]]) - np.amin(midpoints[ind[:j+1]]))
			center_off = midpoints[ind[0]]-params[i]
			print 'bins = %s, true %s width = %.4f, center off = %.5f'%(bins, label_str[i],width, center_off)
			true_width[ibin,i]=width
			
			#for ss in range(len(s_arr)): 
				#fits = Bfit[ss,ibin,:,i+1]
				#hist1, edges = histogram(fits, bins=20)
				#midpoints = edges[:-1] + 0.5 * (edges[2] - edges[1])
				#ind = np.argsort(hist1)[::-1]
				#j = 1
				#tot = sum(hist1)
				#while sum(hist1[ind[:j]]) < 0.683*tot:
					#j += 1
				#width = 0.5*(np.amax(midpoints[ind[:j+1]]) - np.amin(midpoints[ind[:j+1]]))
				#print 'bins = %s, bias s = %s %s width = %.4f'%(bins, s_arr[ss],label_str[i], width)
				#bias_width[ss,ibin,i]=width
	
			ibin += 1
		i+=1
	#savetxt('true_width.ls',true_width)
	#savetxt('bias_width.ls',bias_width.reshape(1,-1))

if err_fisher_vs_fit_ngal30:
	xiuyuan = [0.0078,0.036,0.0057][::-1]
	fisher_err_true = np.ndarray(shape=(10,3))
	fit_err_true = np.ndarray(shape=(10,3))
	fit_err_bias = np.ndarray(shape=(5,10,3))
	i=0
	for bins in bins_arr:
		fisher_err_true[i]=np.genfromtxt(fit_noise_dir+'err_g6.976_ngal30_bins%s.ls'%(bins))
		
		fits = Tfit_ngal30[i,:,1:]
		fit_err_true[i,0] = 0.5*(stats.bayes_mvs(fits[:,0],0.68)[2][1][1]+stats.bayes_mvs(fits[:,0],0.68)[2][1][0])
		
		fit_err_true[i,1] = 0.5*(stats.bayes_mvs(fits[:,1],0.68)[2][1][1]+stats.bayes_mvs(fits[:,1],0.68)[2][1][0])
		
		fit_err_true[i,2] = 0.5*(stats.bayes_mvs(fits[:,2],0.68)[2][1][1]+stats.bayes_mvs(fits[:,2],0.68)[2][1][0])
		
		j=0
		for s in s_arr:
			fits = Bfit_ngal30[j, i,:,1:]
			fit_err_bias[j,i,0] = 0.5*(stats.bayes_mvs(fits[:,0],0.68)[2][1][1]+stats.bayes_mvs(fits[:,0],0.68)[2][1][0])
			
			fit_err_bias[j,i,1] = 0.5*(stats.bayes_mvs(fits[:,1],0.68)[2][1][1]+stats.bayes_mvs(fits[:,1],0.68)[2][1][0])
			
			fit_err_bias[j,i,2] = 0.5*(stats.bayes_mvs(fits[:,2],0.68)[2][1][1]+stats.bayes_mvs(fits[:,2],0.68)[2][1][0])
			
			j+=1
		i+=1

	f = figure(figsize=(8,10))
	label_str=['Om','w','si8']
	for i in range(3): 
		ax = subplot(3,2,i*2+1) #subplots 1,3,5
		y1 = fisher_err_true[:,i]
		y2 = 1.0/40*fit_err_true[:,i]
		ax.plot(bins_arr, y1, linewidth=2, drawstyle = 'steps-mid',label='Fisher err')
		ax.plot(bins_arr, y2, linewidth=2, drawstyle = 'steps-mid',label='Fit err width')
		ax.plot(bins_arr,zeros(len(bins_arr))+xiuyuan[i],label='Yang et el. 2011',linewidth=2)
		leg=ax.legend(loc=0, ncol=1, labelspacing=.2, prop={'size':10})
		leg.get_frame().set_visible(False)
		ax.set_ylabel(label_str[i])
		for ss in range(len(s_arr)):
			ax = subplot(3,2,(i+1)*2) #subplots 2,4,6
			y = 1.0/40*fit_err_bias[ss,:,i]
			ax.plot(bins_arr, y, label='s=%s'%(s_arr[ss]), linewidth=2, drawstyle = 'steps-mid')
			leg=ax.legend(loc=0, ncol=1, labelspacing=.2, prop={'size':10})
			leg.get_frame().set_visible(False)
			ax.set_ylabel(label_str[i])
	subplot(3,2,1)
	title('True Maps (fisher/fit err)')
	subplot(3,2,2)
	title('Bias Maps (fit err only)')
	subplot(3,2,5)
	xlabel('Bins')
	subplot(3,2,6)
	xlabel('Bins')
	savefig(plot_dir+'err_fisher_vs_fit_ngal30.jpg')
	#show()
	close()
	

if margin2params_noise:
	ngal=30
	ibin = 0
	for bins in bins_arr:
		f = figure(figsize=(12, 4))		
		for i in range(3): 
			ax = subplot(1,3,i+1) #subplots 1,3,5
			fits = Tfit_ngal30[ibin,:,i+1]
			hist1, edges = histogram(fits, bins=20)
			midpoints = edges[:-1] + 0.5 * (edges[2] - edges[1])
			ax.plot(midpoints, hist1, label='True', linewidth=2)#, drawstyle = 'steps-mid')
			title('Bins = %s'%(bins))
			ax.set_xlabel(label_str[i])
			
			for ss in range(len(s_arr)):
				ax = subplot(1,3,i+1) #subplots 2,4,6
				fits = Bfit_ngal30[ss,ibin,:,i+1]
				hist1, edges = histogram(fits, bins=20)
				midpoints = edges[:-1] + 0.5 * (edges[2] - edges[1])
				ax.plot(midpoints, hist1, label='s=%s'%(s_arr[ss]), linewidth=2)#, drawstyle = 'steps-mid')
				
				leg=ax.legend(loc=0, ncol=1, labelspacing=.2, prop={'size':10})
				leg.get_frame().set_visible(False)			
			i+=1
		savefig(plot_dir+'margin2params_ngal%s_bins%s.jpg'%(ngal,bins))
		close()
		ibin += 1
		
if margin2params:
	ibin = 0
	for bins in bins_arr:
		f = figure(figsize=(12, 4))		
		for i in range(3): 
			ax = subplot(1,3,i+1) #subplots 1,3,5
			fits = Tfit[ibin,:,i+1]
			hist1, edges = histogram(fits, bins=20)
			midpoints = edges[:-1] + 0.5 * (edges[2] - edges[1])
			ax.plot(midpoints, hist1, label='True', linewidth=2, drawstyle = 'steps-mid')
			title('Bins = %s'%(bins))
			ax.set_xlabel(label_str[i])
			
			for ss in range(len(s_arr)):
				ax = subplot(1,3,i+1) #subplots 2,4,6
				fits = Bfit[ss,ibin,:,i+1]
				hist1, edges = histogram(fits, bins=20)
				midpoints = edges[:-1] + 0.5 * (edges[2] - edges[1])
				ax.plot(midpoints, hist1, label='s=%s'%(s_arr[ss]), linewidth=2, drawstyle = 'steps-mid')
				
				leg=ax.legend(loc=0, ncol=1, labelspacing=.2, prop={'size':10})
				leg.get_frame().set_visible(False)			
			i+=1
		savefig(plot_dir+'margin2params_bins%s_step.jpg'%(bins))
		close()
		ibin += 1

if origin_vs_bins_noise:	
	true_origins = average(Tfit_ngal30[:,:,1:],axis=1) # each bins, averaged 3 params 10*3 matrix
	bias_origins = average(Bfit_ngal30[:,:,:,1:],axis=2) # 5*10*3 matrix

	f = figure(figsize=(8,10))
	label_str=['Om','w','si8']
	for i in range(3): 
		ax = subplot(3,2,i*2+1) #subplots 1,3,5
		y = true_origins[:,i]
		ax.plot(bins_arr, y, linewidth=2, drawstyle = 'steps-mid')
		ax.set_ylabel(label_str[i])
		ax.set_ylim(fidu_params[i]*0.9,fidu_params[i]*1.1)
		for ss in range(len(s_arr)):
			ax = subplot(3,2,(i+1)*2) #subplots 2,4,6
			y = bias_origins[ss,:,i]
			ax.plot(bins_arr, y, label='s=%s'%(s_arr[ss]), linewidth=2, drawstyle = 'steps-mid')
			leg=ax.legend(loc=0, ncol=1, labelspacing=.2, prop={'size':10})
			leg.get_frame().set_visible(False)
			ax.set_ylabel(label_str[i])
	subplot(3,2,1)
	title('True Maps %s ngal%s'%(interpmethod,ngal))
	subplot(3,2,2)
	title('Bias Maps')
	subplot(3,2,5)
	xlabel('Bins')
	subplot(3,2,6)
	xlabel('Bins')
	#show()
	savefig(plot_dir2+'origin_vs_bins_ngal%s_%s.jpg'%(ngal,interpmethod))
	close()
	
if origin_vs_bins:	
	true_origins = average(Tfit[:,:,1:],axis=1) # each bins, averaged 3 params 10*3 matrix
	bias_origins = average(Bfit[:,:,:,1:],axis=2) # 5*10*3 matrix

	f = figure(figsize=(8,10))
	label_str=['Om','w','si8']
	for i in range(3): 
		ax = subplot(3,2,i*2+1) #subplots 1,3,5
		y = true_origins[:,i]
		ax.plot(bins_arr, y, linewidth=2, drawstyle = 'steps-mid')
		ax.set_ylabel(label_str[i])
		for ss in range(len(s_arr)):
			ax = subplot(3,2,(i+1)*2) #subplots 2,4,6
			y = bias_origins[ss,:,i]
			ax.plot(bins_arr, y, label='s=%s'%(s_arr[ss]), linewidth=2, drawstyle = 'steps-mid')
			leg=ax.legend(loc=0, ncol=1, labelspacing=.2, prop={'size':10})
			leg.get_frame().set_visible(False)
			ax.set_ylabel(label_str[i])
	subplot(3,2,1)
	title('True Maps')
	subplot(3,2,2)
	title('Bias Maps')
	subplot(3,2,5)
	xlabel('Bins')
	subplot(3,2,6)
	xlabel('Bins')
	savefig(plot_dir+'origin_vs_bins_spline2.jpg')
	#show()
	close()

##########Method brute: each bins/s/true/bias, fit on points, 3 2D scatter plot####

# true map

if brute_vs_chisqcontour_per_bin:
	for i in range(len(bins_arr)):
		Tpoints = Tfit[i,:,1:]
		bins = bins_arr[i]
		x, y, z = Tpoints.T # 1000 points
		x0, y0, z0 = Tells[i,2:5] # origins
		dx1, dy1, dz1, dth1, dph1 = Tells[i,5:10] # fitting params for brute force
		dx2, dy2, dz2, dth2, dph2 = ChiCont[i,1:6] # fitting params for chisq contour
		
		# project fitted params into the (Om, w, si) axis
		dx1, dy1, dz1 = (RotM(dth1, dph1).I) * np.mat([dx1, dy1, dz1]).T
		dx2, dy2, dz2 = (RotM(dth2, dph2).I) * np.mat([dx2, dy2, dz2]).T
		
		# param pairs for plot loop
		origins = [[x0, y0], [y0, z0], [z0, x0]]
		x_arr = [x, y, z]
		y_arr = [y, z, x]
		
		widths1  = [dx1, dy1, dz1]
		heights1 = [dy1, dz1, dx1]
		widths2  = [dx2, dy2, dz2]
		heights2 = [dy2, dz2, dx2]
		
		angle1 = np.array([dth1, dph1, dph1])*(-pi/180)
		angle2 = np.array([dth2, dph2, dph2])*(-pi/180)
		
		xlims = ((0.15,0.45), (-2.0,0.5), (0.6,1.1))
		ylims = ((-2.0,0.5), (0.6,1.1), (0.15,0.45))
		
		xlabels = ('Om', 'w', 'si8')
		ylabels = ('w', 'si8', 'Om')
		
		legends = ['68.3% pts. enclosed','Chi^2 = 3.53  contour']
		
		f = figure(figsize = (8,10))
		
		for n in range(3):
			title = 'True Map Bins = %s'%(bins)
			ax = f.add_subplot(3,1,n+1) 
			ax.scatter(x_arr[n], y_arr[n], s=1) #s: marker size
			ells12a = Ellipse(xy=origins[n], width=widths1[n], height=heights1[n], angle=angle1[n], label=legends[0], fill = False, linewidth=2, color=mcolors[0:3])
			ells12b = Ellipse(xy=origins[n], width=widths2[n], height=heights2[n], angle=angle2[n], label=legends[1], fill = False, linewidth=2, color=mcolors[3:6])
			
			ax.add_patch(ells12a)
			ax.add_patch(ells12b)
			ells12a.set_clip_box(ax.bbox)
			ells12b.set_clip_box(ax.bbox)
			ax.set_xlabel(xlabels[n])
			ax.set_ylabel(ylabels[n])
			
			leg=ax.legend(title = title, loc=0, ncol=1, labelspacing=.2, prop={'size':10})
			leg.get_frame().set_visible(False) 
			
			ax.set_xlim(xlims[n])
			ax.set_ylim(ylims[n])
			
			grid(True)
		savefig(plot_dir+'true_brute_vs_chisqcontour_bins%s.jpg'%(bins))
		close()
		

if brute_vs_s_per_bin:
	for i in range(len(bins_arr)):
		bins = bins_arr[i]
		x0, y0, z0 = Tells[i,2:5] # origins
		dx1, dy1, dz1, dth1, dph1 = Tells[i,5:10] # fitting params for true points
		
		# project fitted params into the (Om, w, si) axis
		dx1, dy1, dz1 = (RotM(dth1, dph1).I) * np.mat([dx1, dy1, dz1]).T
		
		# param pairs for plot loop
		origins = [[x0, y0], [y0, z0], [z0, x0]]	
		widths1  = [dx1, dy1, dz1]
		heights1 = [dy1, dz1, dx1]
		angle1 = np.array([dth1, dph1, dph1])*(-pi/180)
		
		xlims = ((0.15,0.45), (-2.0,0.5), (0.6,1.1))
		ylims = ((-2.0,0.5), (0.6,1.1), (0.15,0.45))
		
		xlabels = ('Om', 'w', 'si8')
		ylabels = ('w', 'si8', 'Om')
			
		f = figure(figsize = (8,10))

		for j in range(len(s_arr)):
			x0b, y0b, z0b = Bells[i::10,2:5][j]
			dx2, dy2, dz2, dth2, dph2 = Bells[i::10,5:10][j] # fitting params for bias points
			dx2, dy2, dz2 = (RotM(dth2, dph2).I) * np.mat([dx2, dy2, dz2]).T
			widths2  = [dx2, dy2, dz2]
			heights2 = [dy2, dz2, dx2]
			angle2 = np.array([dth2, dph2, dph2])*(-pi/180)
			originsb = [[x0b, y0b], [y0b, z0b], [z0b, x0b]]
			
			for n in range(3):
				ax = f.add_subplot(3,1,n+1) 
				
				ells12b = Ellipse(xy=originsb[n], width=widths2[n], height=heights2[n], angle=angle2[n], label='s = %s'%(s_arr[j]), fill = False, linewidth=1.5, color=mcolors[3*(j+1):3*(j+2)])						
				ax.add_patch(ells12b)
				ells12b.set_clip_box(ax.bbox)
				
				if s_arr[j] == s_arr[-1]: # add true ellipse only once
					title = 'Bins = %s'%(bins)
					
					ells12a = Ellipse(xy=origins[n], width=widths1[n], height=heights1[n], angle=angle1[n], label='True Map', fill = False, linewidth=1.5, color=mcolors[0:3])
					ax.add_patch(ells12a)
					ells12a.set_clip_box(ax.bbox)
					
					ax.set_xlabel(xlabels[n])
					ax.set_ylabel(ylabels[n])
					
					leg=ax.legend(title = title, loc=0, ncol=1, labelspacing=.2, prop={'size':10})
					leg.get_frame().set_visible(False) 
					
					ax.set_xlim(xlims[n])
					ax.set_ylim(ylims[n])
					
					grid(True)
					
		savefig(plot_dir+'bias_brute_s_bins%s.jpg'%(bins))
		#show()
		close()	

if density_sum:
	a=Tfit[-1,:,1:].T

	x,y=a[0],a[1]

	# plot scatter plot
	scatter(x,y,s=1)
	xlim(np.amin(x),np.amax(x))
	ylim(np.amin(y),np.amax(y))
	savefig(plot_dir+'contour_scatter_test.jpg')
	close()

	# plot histogram2d
	H, xedges, yedges = np.histogram2d(y,x, bins=100)
	imshow(H, interpolation='nearest', extent=[yedges[0], yedges[-1], xedges[0], xedges[-1]],aspect='auto', origin='lower')
	colorbar()
	savefig(plot_dir+'contour_hist2d_bins100.jpg')
	close()

	# plot smoothed histogram2d
	import scipy.ndimage as snd
	smooth = lambda kmap, sigma: snd.filters.gaussian_filter(kmap,sigma)
	H10 = smooth(H, 10)

	imshow(H10, interpolation='nearest', extent=[yedges[0], yedges[-1], xedges[0], xedges[-1]],aspect='auto', origin='lower')
	colorbar()
	savefig(plot_dir+'contour_hist2d_bins100_smooth10.jpg')
	close()

	#sort by density, flatten the array first
	H10f = H10.flatten()

	ind = np.argsort(H10f)
	ind = ind[::-1]

	sigma1 = 0.683
	i = 0
	while i < len(ind):
		indsum = ind[:i] # the index of cells to sum
		percent = sum(H10f[indsum])/sum(H10f)
		if percent > sigma1:
			indx = i-1
			break
		i+=1

	cell_ind_flat = ind[:indx]
	cell_ind = np.unravel_index(cell_ind_flat,(100,100))

	# make contour
	mask = np.zeros(shape=(100,100))
	mask[cell_ind] = 1 #mask.nonzero()[0].shape = (1096,)
	H10c = H10*mask

	imshow(H10c, interpolation='nearest', extent=[yedges[0], yedges[-1], xedges[0], xedges[-1]],aspect='auto', origin='lower')
	colorbar()
	savefig('hist2d_bins100_smooth10_contour.jpg')
	close()

if chisq_density_comparison_noise:
	labels = ('Om', 'w', 'si8')
	ngal=30
	
	for ibin in range(len(bins_arr)):
				
		bins = bins_arr[ibin]
		a=Tfit_ngal30[ibin,:,1:].T
		points = a.shape[1]
		
		f = figure(figsize=(10,4))
		for i in ((0,1,2),(1,2,0),(2,0,1)):
			subplot(1,3,i[0]+1)
			# scatter points
			xx, yy = labels[i[0]], labels[i[1]]
			x, y = a[i[0]], a[i[1]]
			pp = 3 # fomr which I computed chisq contour
			xmin, xmax = np.sort(x)[[pp, -pp]]
			ymin, ymax = np.sort(y)[[pp, -pp]]
			
			############ density contour ############
			db = 100 # density bins
			H, xedges, yedges = np.histogram2d(x, y, bins=db, range = [[xmin, xmax],[ymin, ymax]])
			X = (xedges[1:] + xedges[:-1])/2
			Y = (yedges[1:] + yedges[:-1])/2
			H10 = smooth(H, 10).T # changed - added .T
			H10f = H10.flatten()
			ind = np.argsort(H10f)
			ind = ind[::-1]

			sigma1 = 0.683
			ii=0
			while ii < len(ind):
				indsum = ind[:ii] # the index of cells to sum
				percent = sum(H10f[indsum])/points
				if percent > sigma1:
					break
				ii += 1

			cell_ind_flat = ind[:ii]
			cell_ind = np.unravel_index(cell_ind_flat,(db,db))
			# mask the region outside of interest
			H10c = np.zeros(shape=(db,db))
			H10c[cell_ind] = 1
			
			
			# chi^2 contour
			C = genfromtxt(fit_noise_dir+'chisqgrid_true_ngal%s_g6.976_bins%s.ls'%(ngal,bins))
			res = round(C.shape[0]**(1.0/3.0))
			C = C.reshape(res,res,res)
			C2d = np.amin(C,axis=i[-1])
			X2 = linspace(xmin,xmax,res)
			Y2 = linspace(ymin,ymax,res)
			
			# plot scatter & 2 contours
			
			scatter(x,y,s=1,label = 'scatter')
			contour(X,Y,H10c,colors='b')
			##contour(X2,Y2,C2d,(3.53,),colors='r')
			contour(X2,Y2,C2d,(2.3,),colors='r')
			
			title('Bins = %s'%(bins))
			xlabel(xx)
			ylabel(yy)
		savefig(plot_dir+'contour_true_ngal%s_bins%s.jpg'%(ngal,bins))
		close()

if chisq_density_comparison_noise_bias:
	labels = ('Om', 'w', 'si8')
	ngal=30
	iis=2
	s=s_arr[iis]
	for ibin in range(len(bins_arr)):
				
		bins = bins_arr[ibin]
		a=Bfit_ngal30[iis,ibin,:,1:].T
		points = a.shape[1]
		
		f = figure(figsize=(10,4))
		for i in ((0,1,2),(1,2,0),(2,0,1)):
			subplot(1,3,i[0]+1)
			# scatter points
			xx, yy = labels[i[0]], labels[i[1]]
			x, y = a[i[0]], a[i[1]]
			pp = 3 # fomr which I computed chisq contour
			xmin, xmax = np.sort(x)[[pp, -pp]]
			ymin, ymax = np.sort(y)[[pp, -pp]]
			
			############ density contour ############
			db = 100 # density bins
			H, xedges, yedges = np.histogram2d(x, y, bins=db, range = [[xmin, xmax],[ymin, ymax]])
			X = (xedges[1:] + xedges[:-1])/2
			Y = (yedges[1:] + yedges[:-1])/2
			H10 = smooth(H, 10).T # changed - added .T
			H10f = H10.flatten()
			ind = np.argsort(H10f)
			ind = ind[::-1]

			sigma1 = 0.683
			ii=0
			while ii < len(ind):
				indsum = ind[:ii] # the index of cells to sum
				percent = sum(H10f[indsum])/points
				if percent > sigma1:
					break
				ii += 1

			cell_ind_flat = ind[:ii]
			cell_ind = np.unravel_index(cell_ind_flat,(db,db))
			# mask the region outside of interest
			H10c = np.zeros(shape=(db,db))
			H10c[cell_ind] = 1
			
			
			# chi^2 contour
			#C = genfromtxt(fit_noise_dir+'chisqgrid_true_ngal%s_g6.976_bins%s.ls'%(ngal,bins))
			#res = round(C.shape[0]**(1.0/3.0))
			#C = C.reshape(res,res,res)
			#C2d = np.amin(C,axis=i[-1])
			#X2 = linspace(xmin,xmax,res)
			#Y2 = linspace(ymin,ymax,res)
			
			# plot scatter & 2 contours
			
			scatter(x,y,s=1,label = 'scatter')
			contour(X,Y,H10c,colors='b')
			##contour(X2,Y2,C2d,(3.53,),colors='r')
			#contour(X2,Y2,C2d,(2.3,),colors='r')
			
			title('Bins = %s'%(bins))
			xlabel(xx)
			ylabel(yy)
		savefig(plot_dir+'contour_bias_s%s_ngal%s_bins%s.jpg'%(s,ngal,bins))
		close()
		
if chisq_density_comparison:
	labels = ('Om', 'w', 'si8')
	
	for ibin in range(len(bins_arr)):
				
		bins = bins_arr[ibin]
		a=Tfit[ibin,:,1:].T
		points = a.shape[1]
		
		f = figure(figsize=(10,4))
		for i in ((0,1,2),(1,2,0),(2,0,1)):
			subplot(1,3,i[0]+1)
			# scatter points
			xx, yy = labels[i[0]], labels[i[1]]
			x, y = a[i[0]], a[i[1]]
			pp = 3 # fomr which I computed chisq contour
			xmin, xmax = np.sort(x)[[pp, -pp]]
			ymin, ymax = np.sort(y)[[pp, -pp]]
			
			############ density contour ############
			db = 100 # density bins
			H, xedges, yedges = np.histogram2d(x, y, bins=db, range = [[xmin, xmax],[ymin, ymax]])
			X = (xedges[1:] + xedges[:-1])/2
			Y = (yedges[1:] + yedges[:-1])/2
			H10 = smooth(H, 10).T # changed - added .T
			H10f = H10.flatten()
			ind = np.argsort(H10f)
			ind = ind[::-1]

			sigma1 = 0.683
			ii=0
			while ii < len(ind):
				indsum = ind[:ii] # the index of cells to sum
				percent = sum(H10f[indsum])/points
				if percent > sigma1:
					break
				ii += 1

			cell_ind_flat = ind[:ii]
			cell_ind = np.unravel_index(cell_ind_flat,(db,db))
			# mask the region outside of interest
			H10c = np.zeros(shape=(db,db))
			H10c[cell_ind] = 1
			
			
			# chi^2 contour
			C = genfromtxt(fit_dir+'chisqgrid_true_bins%s.ls'%(bins))
			res = round(C.shape[0]**(1.0/3.0))
			C = C.reshape(res,res,res)
			C2d = np.amin(C,axis=i[-1])
			X2 = linspace(xmin,xmax,res)
			Y2 = linspace(ymin,ymax,res)
			
			# plot scatter & 2 contours
			
			scatter(x,y,s=1,label = 'scatter')
			contour(X,Y,H10c,colors='b')
			#contour(X2,Y2,C2d,(3.53,),colors='r')
			contour(X2,Y2,C2d,(2.3,),colors='r')
			
			
			title('Bins = %s'%(bins))
			xlabel(xx)
			ylabel(yy)
			
			##xlim(xmin, xmax)
			##ylim(ymin, ymax)
			
			####savefig(plot_dir+'contour_true_bins%s_%s_%s.jpg'%(bins, xx, yy))
			####close()
		#show()
		savefig(plot_dir+'contour_true_spline2_bins%s.jpg'%(bins))
		close()

if chisq_density_comparison_xy:
	plot_dir = '/Users/jia/Documents/weaklensing/magbias/plot/xy_z1/' # folder to save plots
	fit_dir = '/Users/jia/Documents/weaklensing/magbias/fit_xy_z1/' # folder to save fit files
	
	labels = ('Om', 'w', 'si8')
	Tfit = genfromtxt(fit_dir+'true_fit_g9.865_bins_xy.ls')
	a = Tfit[:,1:].T
	points = a.shape[1]
	
	f = figure(figsize=(10,4))
	for i in ((0,1,2),(1,2,0),(2,0,1)):
		subplot(1,3,i[0]+1)
		# scatter points
		xx, yy = labels[i[0]], labels[i[1]]
		x, y = a[i[0]], a[i[1]]
		pp = 3 # fomr which I computed chisq contour
		xmin, xmax = np.sort(x)[[pp, -pp]]
		ymin, ymax = np.sort(y)[[pp, -pp]]
		
		############ density contour ############
		db = 100 # density bins
		H, xedges, yedges = np.histogram2d(x, y, bins=db, range = [[xmin, xmax],[ymin, ymax]])
		X = (xedges[1:] + xedges[:-1])/2
		Y = (yedges[1:] + yedges[:-1])/2
		H10 = smooth(H, 10).T # changed - added .T
		H10f = H10.flatten()
		ind = np.argsort(H10f)
		ind = ind[::-1]

		sigma1 = 0.683
		ii=0
		while ii < len(ind):
			indsum = ind[:ii] # the index of cells to sum
			percent = sum(H10f[indsum])/points
			if percent > sigma1:
				break
			ii += 1

		cell_ind_flat = ind[:ii]
		cell_ind = np.unravel_index(cell_ind_flat,(db,db))
		# mask the region outside of interest
		H10c = np.zeros(shape=(db,db))
		H10c[cell_ind] = 1
		
		
		# chi^2 contour
		C = genfromtxt(fit_dir+'chisqgrid_true_bins_xy.ls')
		res = round(C.shape[0]**(1.0/3.0))
		C = C.reshape(res,res,res)
		C2d = np.amin(C,axis=i[-1])
		X2 = linspace(xmin,xmax,res)
		Y2 = linspace(ymin,ymax,res)
		
		# plot scatter & 2 contours
		
		scatter(x,y,s=1,label = 'scatter')
		contour(X,Y,H10c,colors='b')
		#contour(X2,Y2,C2d,(3.53,),colors='r')
		contour(X2,Y2,C2d,(2.3,),colors='r')
		
		xlabel(xx)
		ylabel(yy)
		
	savefig(plot_dir+'contour_true_spline_bins_xy.jpg')
	show()
	#close()	
	

if smoothcompxy:
	smoothdir = '/Users/jia/magbias/nonoise/fit_xy_g6/smoothcomp/'
	peaksxy=genfromtxt(smoothdir+'peaksxy.ls')
	peaksxynoise=genfromtxt(smoothdir+'peaksxynoise.ls')
	peaksjl=genfromtxt(smoothdir+'peaksjl.ls')
	peaksjlnoise=genfromtxt(smoothdir+'peaksjlnoise.ls')
	peaksjlnoise_sqrt=genfromtxt(smoothdir+'peaksjlnoise_sqrt.ls')
	binn=linspace(-0.01859393,0.16237037,100)
	b1,a1=histogram(peaksjl,bins=binn)
	b2,a2=histogram(peaksxy,bins=binn)
	b3,a3=histogram(peaksxynoise,bins=binn)
	b4,a4=histogram(peaksjlnoise,bins=binn)
	b5,a5=histogram(peaksjlnoise_sqrt,bins=binn)
	figure()
	subplot(2,1,1)
	plot(a1[1:],b1,linewidth=2, drawstyle = 'steps',label='Jia noiseless 1000r (%s)'%(len(peaksjl)))
	plot(a4[1:],b4,linewidth=2, drawstyle = 'steps',label='Jia with noise 1000r (%s)'%(len(peaksjlnoise)))
	plot(a5[1:],b5,linewidth=2, drawstyle = 'steps',label='Jia with sqrt(noise) 1000r (%s)'%(len(peaksjlnoise_sqrt)))
	plot(a2[1:],b2,linewidth=2, drawstyle = 'steps',label='Xiuyuan noiseless 1000r (%s)'%(len(peaksxy)))
	plot(a3[1:],b3,linewidth=2, drawstyle = 'steps',label='Xiuyuan with noise 1000r (%s)'%(len(peaksxynoise)))
	title('Peak distribution')
	legend(loc=0, ncol=1, labelspacing=.2, prop={'size':10})
	
	subplot(2,1,2)
	plot(a4[1:],b4-b1,linewidth=2, drawstyle = 'steps',label='Jia 1000r (%s)'%(len(peaksjlnoise)))
	plot(a4[1:],b5-b1,linewidth=2, drawstyle = 'steps',label='Jia (sqrt) 1000r (%s)'%(len(peaksjlnoise_sqrt)))
	plot(a3[1:],b3-b2,linewidth=2, drawstyle = 'steps',label='Xiuyuan 1000r (%s)'%(len(peaksxynoise)))
	
	title('Noise - Noiseless distribution')
	xlabel('kappa')
	
	legend(loc=0, ncol=1, labelspacing=.2, prop={'size':10})
	savefig(smoothdir+'compsmooth_r1000_bins100_c.jpg')
	close()


if density_contour:
	# test if density contour is truely surrounding the points
	labels = ('Om', 'w', 'si8')
	
	for ibin in range(len(bins_arr)):
				
		bins = bins_arr[ibin]
		a=Tfit[ibin,:,1:].T
		points = a.shape[1]
		
		f = figure(figsize=(10,4))
		for i in ((0,1,2),(1,2,0),(2,0,1)):
			subplot(1,3,i[0]+1)
			# scatter points
			xx, yy = labels[i[0]], labels[i[1]]
			x, y = a[i[0]], a[i[1]]
			px = 10
			xmin, xmax = np.sort(x)[[px,-px]]
			ymin, ymax = np.sort(y)[[px,-px]]
			#xmin, xmax = np.sort(x)[[0,-1]]
			#ymin, ymax = np.sort(y)[[0,-1]]
			
			############ density contour ############
			db = 100 # density bins
			smoothscale = 10
			
			H, xedges, yedges = np.histogram2d(x, y, bins=db)#, range = [[xmin, xmax],[ymin, ymax]])
			X = (xedges[1:] + xedges[:-1])/2
			Y = (yedges[1:] + yedges[:-1])/2
			H10 = smooth(H, smoothscale).T ## important, to transpose, somehow imshow is transposed..

			# plot scatter & 2 contours
			imshow(H10, interpolation='nearest', extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]],aspect='auto', origin='lower')
			contour(X,Y,H10)
			#colorbar()
			scatter(x,y,s=1)
						
			title('Bins = %s'%(bins))
			xlabel(xx)
			ylabel(yy)
			
			#xlim(xmin, xmax)
			#ylim(ymin, ymax)
			
			#savefig(plot_dir+'contour_true_bins%s_%s_%s.jpg'%(bins, xx, yy))
			#close()
		
		savefig(plot_dir+'contour_bins%s.jpg'%(bins))
		close()
		
		
if test_kinks:
	# test if there's points preferentially select some value, mostly due to extrapolation & interpolation methods
	w_fits=Tfit[:,:,2]
	Om_fits=Tfit[:,:,1]
	si_fits=Tfit[:,:,3]
	#w2=w_fits.flatten()

	i = 0
	for x in w_fits:
		hist(x,bins=100,range = [-1.3, -0.7])
		title ('Bins = %s'%(bins_arr[i]))
		xlabel('w')
		savefig(plot_dir+'w_kink_bins%s.jpg'%(bins_arr[i]))
		close()
		i+=1

	i = 0
	for x in Om_fits:
		hist(x,bins=100,range = [0.22, 0.3])
		title ('Bins = %s'%(bins_arr[i]))
		xlabel('Om')
		savefig(plot_dir+'Om_kink_bins%s.jpg'%(bins_arr[i]))
		close()
		i+=1

	i = 0
	for x in si_fits:
		hist(x,bins=100,range = [0.74, 0.81])
		title ('Bins = %s'%(bins_arr[i]))
		xlabel('si')
		savefig(plot_dir+'si_kink_bins%s.jpg'%(bins_arr[i]))
		close()
		i+=1


########## 3d plot #######
def plot3d (s, bins):
	from mpl_toolkits.mplot3d import Axes3D
	import animate3dplot
	bothfits = fitdistrib (s,bins)
	test = ('True','Bias')
	for n in (0,1):
		fit = bothfits[n]
		chi, x, y, z = fit.T
		
		fig = plt.figure()
		ax = fig.add_subplot(111, projection='3d',xlabel='Om', ylabel='w', zlabel='Sigma8')
		ax.scatter(x, y, z, c=chi)
		ax.text2D(0.7, 0.89,'s = %s, bins = %s'%(s,bins),color='black', transform=ax.transAxes)
		if n == 0:
			ax.text2D(0.7, 0.85,'True (0.26, -1.0, 0.798)',color='red', transform=ax.transAxes)#color='red'
		if n == 1:
			ax.text2D(0.7, 0.85,'Bias (%.3f, %.3f, %.4f)'%(avg_bias_fit[1],avg_bias_fit[2],avg_bias_fit[3]),color='blue', transform=ax.transAxes) # 0.7, 0.85
		
		angles = np.linspace(0,360,21)[:-1]
		animate3dplot.rotanimate(ax, angles,'%s3d_s%s_g%s_bins%s.gif'%(test[n],s,sigma,bins),delay=60,width=8,height=6)	
		