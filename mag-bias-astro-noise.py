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
#python s bins ngal interpmethod iseed=73000
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

s_arr = (0.0, 0.5, 1.0, 1.5, 2.0)
bins_arr =  (10, 20, 30, 40, 50, 60, 70, 80, 90, 100)
initguess = (0.27, -1.1, 0.76) # initial guesses for fitter
sigma = 6.976
#sigma = 9.865

low, high = -0.02, 0.19 # -0.15 , 0.25 kappa range
kmapnum=1000 # number of maps per cosmology
fidu_params = np.array([0.26, -1.0, 0.798])

if len(sys.argv)>5:
	iseed=int(sys.argv[5])
	peak_dir = '/direct/astro+astronfs01/workarea/jia/peaks_noise_seed%s/'%(iseed) # folder to save peak files
	fit_dir = '/direct/astro+astronfs01/workarea/jia/fit_noise_seed%s/'%(iseed)
	peaks_mat_dir = fit_dir
else:
	iseed=0
	peak_dir = '/direct/astro+u/jia/magbias/peaks_noise/' # folder to save peak files
	peaks_mat_dir = '/direct/astro+u/jia/magbias/peaks_mat_noise/'
	fit_dir = '/direct/astro+u/jia/magbias/fit_noise/%s/'%(interpmethod) # folder to save fit files


		
kmap_dir = '/direct/astro+u/jia/kmap/' # where all the map directories locate
plot_dir = '/direct/astro+u/jia/magbias/plot/' # folder to save plots
analytical_dir = '/direct/astro+u/jia/magbias/fit_noise_analytical/%s/'%(interpmethod)

#peak_dir = '/direct/astro+u/jia/magbias/peaks_noise_noseed/' # folder to save peak files
#fit_dir = '/direct/astro+u/jia/magbias/fit_noise_noseed/' # folder to save fit files


#kmap_dir = '/Users/jia/Documents/weaklensing/map2/' # where all the map directories locate
#peak_dir = '/Users/jia/Documents/weaklensing/magbias/peaks/' # folder to save peak files
#plot_dir = '/Users/jia/Documents/weaklensing/magbias/plot/' # folder to save plots
#fit_dir = '/Users/jia/Documents/weaklensing/magbias/fit_noise/' # folder to save fit files

cosmo_all = os.listdir(kmap_dir) #a list of 7 cosmologies, suppose all the map folders are in kmap_dir

########### equations & functions that takes (bins, s, sigma) ########
#bin_size = (high - low) / bins
#hist_bins = np.linspace (low+bin_size/2, high-bin_size/2, bins) # center of each bin

smooth = lambda kmap, sigma: snd.filters.gaussian_filter(kmap,sigma)

biasnoisesmooth = lambda kmap, noise, sigma, s: snd.filters.gaussian_filter(kmap*(1+(5*s-2)*kmap)+noise,sigma)/snd.filters.gaussian_filter(1+(5*s-2)*kmap,sigma)

#biasmooth = lambda kmap, sigma, s: snd.filters.gaussian_filter(kmap*(1+(5*s-2)*kmap),sigma)/snd.filters.gaussian_filter(1+(5*s-2)*kmap,sigma)


def peaklist_noise(fitsfile, s=s, bins=bins, sigma=sigma, bias=0):
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


def batch_kmap_noise(kmap_dir,cosmo_num, s=s, bins=bins, sigma=sigma, bias=0):
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
			peaks = peaklist_noise(fitsfile, s=s, bins=bins, sigma=sigma, bias=bias)
			peak_matrix[n] = np.histogram(peaks, range=(low,high), bins=bins)[0]
			print 'Now making histogram for fitsfile # bias=',str(bias),str(n)
		np.savetxt(file_matrix, peak_matrix)
	return peak_matrix
	

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


################ chi-square fitting ##############################
#cov_mat = np.cov(true_mat[contr],rowvar=0)
cov_mat = np.cov(true_mat_noise,rowvar=0)
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
def fitdistrib_noise (mat_noise, s=s, bins=bins, sigma=sigma, bias=0):
	if bias:
		fit_filename = fit_dir+'bias_fit_s%s_g%s_ngal%s_bins%s.ls'%(s,sigma,int(ngal_arcmin),bins)
	else:
		fit_filename = fit_dir+'true_fit_g%s_ngal%s_bins%s.ls'%(sigma,int(ngal_arcmin),bins)
	if os.path.isfile(fit_filename):
		fit = np.genfromtxt(fit_filename)
	else:
		fit = np.ndarray(shape=(kmapnum,4))
		n=0
		while n < kmapnum:
			ikmap = mat_noise[n] #changed from bias_mat[contr,n]
			a = minimizechisq (ikmap)
			fit[n] = a
			print 'fitting to map #',n,'bias=',bias
			n+=1		
		savetxt(fit_filename,fit)
	return fit

#true_fit_noise = fitdistrib_noise (true_mat_noise, s=s, bins=bins, sigma=sigma, bias=0)
#bias_fit_noise = fitdistrib_noise (bias_mat_noise, s=s, bins=bins, sigma=sigma, bias=1)


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
X=np.mat([dNdOm,dNdw,dNdsi])

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
			Y=np.mat(ibmap-fidu_mat[contr])
			del_p=((X*cov_inv*X.T).I)*(X*cov_inv*Y.T)		
			initguess=np.squeeze(np.array(del_p.T))+fidu_params
			fit_analytical[i,0]=chisq(initguess, ibmap)
			fit_analytical[i,1:]=np.array(initguess)
		savetxt(fit_filename,fit_analytical)
	return fit_analytical

true_fit_analytical=analytical_fits(true_mat_noise, s=s, bins=bins, sigma=sigma, bias=0)
bias_fit_analytical=analytical_fits(bias_mat_noise, s=s, bins=bins, sigma=sigma, bias=1)

###########################################
def chisqgrid(chisq0=3.53,Ns=60):
	lower = np.sort(true_fit_noise,axis=0)[3,1:]
	upper = np.sort(true_fit_noise,axis=0)[-3,1:]
	filename = fit_dir+'chisqgrid_true_ngal%s_g%s_bins%s.ls'%(int(ngal_arcmin),sigma,bins)
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

#print 'Now searching in chisq grid, Ns = 60'
#chisqgrid,lower,upper = chisqgrid(chisq0=3.53,Ns=60)

############# calculate fisher matrix #########

#cov_mat = np.cov(true_mat_noise,rowvar=0)
#cov_inv = np.mat(cov_mat).I
#savetxt(fit_dir+'covmat_g%s_ngal%s_bins%s.ls'%(sigma,int(ngal_arcmin),bins),cov_mat)

#errfile=fit_dir+'err_g%s_ngal%s_bins%s.ls'%(sigma,int(ngal_arcmin),bins)
#if os.path.isfile(errfile):
	#err = genfromtxt(errfile)
#else:
	#hi_arr = avg_noise_all[[hi_si,hi_w,hi_Om]]
	#lo_arr = avg_noise_all[[lo_si,lo_w,lo_Om]]
	#contr_arr = avg_noise_all[contr]

	#devhi = np.array([0.052,0.2,0.03])
	#devlo = np.array([0.048,0.2,0.03])

	#Fisher = np.ndarray(shape=(3,3))
	#for x in range(3):
		#for y in range(3):
			##Na = np.mat(arr[x]/dev[x])
			##Nb = np.mat(arr[y]/dev[y])
			#Nahi = np.mat((hi_arr[x]-contr_arr)/devhi[x])
			#Nbhi = np.mat((hi_arr[y]-contr_arr)/devhi[y])
			#Nalo = -np.mat((lo_arr[x]-contr_arr)/devlo[x])
			#Nblo = -np.mat((lo_arr[y]-contr_arr)/devlo[y])
			#Na = 0.5* (Nahi+Nalo)
			#Nb = 0.5* (Nbhi+Nblo)
			#M = Na.T*Nb+Nb.T*Na
			#Fisher[x,y] = 0.5 * trace(cov_inv*M)

	#err = sqrt(mat(Fisher)**-1)[[0,1,2],[0,1,2]]/sqrt(20000./3.46**2)
	#savetxt(errfile,err)



#xiuyuan = np.array([0.0078,0.036,.0057])
#print 'jia err',err
#print 'xiuyuan err', xiuyuan

print 'DONE-DONE-DONE-DONE-DONE-DONE'