import numpy as np
from scipy import *
from scipy import interpolate
from pylab import *
import pyfits
import scipy.ndimage as snd
import os #math
import sys
from matplotlib.patches import Ellipse
import matplotlib.pyplot as plt
import matplotlib.patches as patches
iseed=1

description='''
X: Masked 1000 maps
C: 5sim 0-499
'''#%(iseed)#X 5sim random 500 (seed%s)
#dNtitle='MaskXsameMaps'
#dNtitle='ScaleXsameMaps'
dNtitle='MaskXdiffMapsSeed%s'%(iseed)
#dNtitle='ScaleXdiffMaps1000'
#dNtitle='MaskXdiffMaps'
#dNtitle='ScaleXdiffMaps'
s = 0.4
bins = 20
num=bins
ngal_arcmin = 30
interpmethod = 'fw'
ngal=12.0*60**2*ngal_arcmin/2048.0**2
sigma = 6.976
mask = 0.7

low, high = -0.02, 0.19 # -0.15 , 0.25 kappa range
kappa_arr=np.linspace(low,high,bins)[:num]	
kmapnum=1000 # number of maps per cosmology
fidu_params = np.array([0.26, -1.0, 0.798])
label_latex=['\Omega_m','w','\sigma_8']

fit_dir = '/Users/jia/magbias/fit_noise_sim45_mask2_NoScaleX/'
peaks_mat_dir_45sim = '/Users/jia/magbias/peaks_mat_noise_sim45_mask2/'
peaks_mat_dir_5sim = '/Users/jia/magbias/peaks_mat_noise_mask2/'
kmap_dir = '/Users/jia/weaklensing/map2/'
cosmo_all = os.listdir(kmap_dir)
plot_dir = '/Users/jia/magbias/dNplot/'

equalwidthbins = 0

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
##########################################################

if equalwidthbins:
	fidu500bins=genfromtxt('/Users/jia/magbias/peaks_mat_noise_mask/peak_matrix_true_1000m-512b240_Om0.260_Ol0.740_w-1.000_ns0.960_si0.798_g6.976_ngal30_bins500.ls')
	fidu500bins=average(fidu500bins,axis=0)
	fiduCumSum=cumsum(fidu500bins)
	Count1=sum(fidu500bins)/bins
	cuts=Count1*(arange(bins)+1)

	def find_nearest(array, value):
		n = [abs(i-value) for i in array]
		idx = n.index(min(n))
		return idx#array[idx]	
	idx_arr=np.zeros(bins,dtype=int)
	j=0
	for i in cuts:
		idx=find_nearest(fiduCumSum,i)
		idx_arr[j]=idx
		print idx
		j+=1
	
	def rebin(array, Redges,axis=1):
		#always bin from batchmatnoise anyways, so it's 1000x500 arrays
		shape=np.array(array.shape)
		shape[axis]=bins
		newarray=np.ndarray(shape=shape)
		for i in range(len(Redges)):
			if i==0:
				l=0
			else:
				l=Redges[i-1]+1
			r=Redges[i]+1
			print l,r
			newarray[:,i]=sum(array[:,l:r],axis=axis)
		return newarray
	
	def batch_kmap_noise(peaks_mat_dir = peaks_mat_dir_45sim, cosmo_num = contr, bins=500, mask=mask):#add mask=0
		cosmo = cosmo_all[cosmo_num]
		kmaps = os.listdir(kmap_dir+cosmo)
		if mask == 0:
			file_matrix = peaks_mat_dir+'peak_matrix_true_'+str(kmapnum)+cosmo+'_g%s_ngal%s_bins%s.ls'%(sigma,int(ngal_arcmin),bins)
		else:#add mask batch maps
			file_matrix = peaks_mat_dir+'peak_matrix_true_'+str(kmapnum)+cosmo+'_g%s_ngal%s_bins%s_mask%s.ls'%(sigma,int(ngal_arcmin),bins,mask)
		peak_matrix = np.genfromtxt(file_matrix)
		peak_matrix = rebin(peak_matrix,idx_arr)

		return peak_matrix
	kappa_arr=np.linspace(low,high,500)[idx_arr]
	
else:
	def batch_kmap_noise(peaks_mat_dir = peaks_mat_dir_45sim, cosmo_num = contr, bins=bins, mask=mask):#add mask=0
		cosmo = cosmo_all[cosmo_num]
		kmaps = os.listdir(kmap_dir+cosmo)
		if mask == 0:
			file_matrix = peaks_mat_dir+'peak_matrix_true_'+str(kmapnum)+cosmo+'_g%s_ngal%s_bins%s.ls'%(sigma,int(ngal_arcmin),bins)
		else:#add mask batch maps
			file_matrix = peaks_mat_dir+'peak_matrix_true_'+str(kmapnum)+cosmo+'_g%s_ngal%s_bins%s_mask%s.ls'%(sigma,int(ngal_arcmin),bins,mask)
		peak_matrix = np.genfromtxt(file_matrix)

		return peak_matrix



def FisherMat (X, cov_inv):
	Fisher = np.ndarray(shape=(3,3))
	for x in range(3):
		for y in range(3):
			#Nahi = np.mat((hi_arr[x]-contr_arr)/devhi[x])
			#Nbhi = np.mat((hi_arr[y]-contr_arr)/devhi[y])
			#Nalo = -np.mat((lo_arr[x]-contr_arr)/devlo[x])
			#Nblo = -np.mat((lo_arr[y]-contr_arr)/devlo[y])
			#Na = 0.5* (Nahi+Nalo)
			#Nb = 0.5* (Nbhi+Nblo)
			Na=mat(X[x])
			Nb=mat(X[y])
			M = Na.T*Nb+Nb.T*Na
			Fisher[x,y] = 0.5 * trace(cov_inv*M)
	err = diagonal(sqrt(mat(Fisher).I))#/sqrt(20000.0/12.0)
	return Fisher,err
###### here we get model fidu map that's already calculated ###
all_fits = np.ndarray(shape=(4,1000,4))
X_arr = np.ndarray(shape=(4,3,num))
dN_arr =  np.ndarray(shape=(4,3,num))
cov_inv_arr = np.ndarray(shape=(4,num,num))
FisherMat_arr = np.ndarray(shape=(4,3,3))
err_arr = np.ndarray(shape=(4,3))
j = 0
for mask in (0, 0.284, 0.531, 0.7):

	mapcutoff=500

	true_mat_noise_all = np.ndarray(shape=(len(cosmo_all), kmapnum, bins))
	for i in np.arange(len(cosmo_all)):
		# (1) 
		##### mask dN / X
		true_mat_noise_all[i] = batch_kmap_noise(peaks_mat_dir = peaks_mat_dir_5sim, cosmo_num = i, mask=mask)
		##### scale dN / X
		#true_mat_noise_all[i] = batch_kmap_noise(peaks_mat_dir = peaks_mat_dir_5sim, cosmo_num = i, mask=0)
		#true_mat_noise_all[i]*=(1.0-mask)

	# (2) covmat
	#5sim 500-999
	#cov_mat = np.cov(true_mat_noise_all[contr][mapcutoff:,:num],rowvar=0)
	#5sim 0-499
	#cov_mat = np.cov(true_mat_noise_all[contr][:mapcutoff,:num],rowvar=0)
	cov_mat = np.cov(true_mat_noise_all[contr][:,:num],rowvar=0)
	
	cov_inv = np.mat(cov_mat).I
	
	# X which maps
	#sim5, use random 500
	seed(iseed)
	maps=randint(0,999,500)
	#fidu_mat = average(true_mat_noise_all[:,maps],axis=1)
		
	#sim5, use #500-999
	#fidu_mat = average(true_mat_noise_all[:,mapcutoff:],axis=1)
	#sim5, use #1000 maps
	fidu_mat = average(true_mat_noise_all,axis=1)

	# 45sim for footprint & MC
	true_mat_noise = batch_kmap_noise(peaks_mat_dir = peaks_mat_dir_45sim, cosmo_num = contr, mask=mask)#[:mapcutoff]
	fidu_mat_sim45 = average(true_mat_noise,axis=0)

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
	
	


	def analytical_fits(peak_mat, bins=bins,mask=mask):
		#fit_filename = analytical_dir+'true_analytical_g%s_ngal%s_bins%s_mask%s.ls'%(sigma,int(ngal_arcmin),bins,mask)
		#if os.path.isfile(fit_filename):
			#fit_analytical = np.genfromtxt(fit_filename)
		#else:	
		matlen=len(peak_mat)
		fit_analytical=np.ndarray(shape=(matlen,4))
		for i in range(matlen):
			#print 'analytical fitting to map #, cosmo',i
			ibmap=peak_mat[i]
			Y=np.mat(ibmap[:num]-fidu_mat_sim45[:num])
			del_p=((X*cov_inv*X.T).I)*(X*cov_inv*Y.T)		
			initguess=np.squeeze(np.array(del_p.T))+fidu_params
			del_N=Y-del_p.T*X
			fit_analytical[i,0]=(matlen-2.0-bins)/(matlen-1.0)*del_N*cov_inv*del_N.T#chisq
			fit_analytical[i,1:]=np.array(initguess)
		#savetxt(fit_filename,fit_analytical)
		return fit_analytical

	all_fits[j]=analytical_fits(true_mat_noise)
	X_arr[j]=X
	dN_arr[j]=np.mat([dNOm,dNw,dNsi])[:,:num]
	cov_inv_arr[j]=cov_inv
	FisherMat_arr[j], err_arr[j]=FisherMat(X,cov_inv)
	j+=1


def plotEllipse(pos,P,edge='r',ls='solid',label=None,mask=0):
	Area=sqrt(1.0-mask)#1.0#sqrt(20000./12)
	#print Area
	#print Area
	U, s, Vh = svd(P)
	orient = math.atan2(U[1,0],U[0,0])*180/pi
	ellipsePlot = Ellipse(xy=pos, width=2.0*math.sqrt(s[0])*Area, height=2.0*math.sqrt(s[1])*Area, angle=orient, fill = False,linewidth=1.5,edgecolor=edge,ls=ls,label=label)
	ax = gca()
	ax.add_patch(ellipsePlot)
	return ellipsePlot

seed(10)
colors=rand(30).reshape(10,3)

def eigen(P): 
	a = P[0,0]
	b = P[0,1]
	c = P[1,0]
	d = P[1,1]
	eig1 = 0.5*(a+d+sqrt((a+d)**2-4*(a*d-b*c)))
	eig2 = 0.5*(a+d-sqrt((a+d)**2-4*(a*d-b*c)))
	tr = a+d
	det = a*d-b*c
	return eig1, eig2, tr, det

###########plot covariance matrix #########
#f=figure(figsize=(5,8))
#ax1=f.add_subplot(211)
#ax2=f.add_subplot(212)
#im1=ax1.imshow(cov_mat,interpolation='nearest')#, origin='lower')
#cax1 = f.add_axes([0.85, 0.6, 0.03, 0.35])#from left, bottom, width, height
#f.colorbar(im1,cax=cax1)
#im2=ax2.imshow(cov_inv,interpolation='nearest')#, origin='lower')
#cax2 = f.add_axes([0.85, 0.1, 0.03, 0.35])
#f.colorbar(im2,cax=cax2)

#ax1.set_title('$C_{ij}$')
#ax2.set_title('$C^{-1}$')
##show()
#savefig(plot_dir+'cov_mat.jpg')
#close()


#####looking at diagnoal terms #######
#chisq = lambda cov_inv, dN: dN*mat(cov_inv)*dN.T
#I=identity(num)
#cov_inv_D_arr = np.array(cov_inv_arr)*np.array([I,I,I,I])

#Fchisq = np.ndarray(shape=(4,3))
#Dchisq = np.ndarray(shape=(4,3))
#for i in range (4):
	#for j in range(3):
		#dN = mat(dN_arr[i,j])
		#Fchisq[i,j]=chisq(cov_inv_arr[i],dN)
		#Dchisq[i,j]=chisq(cov_inv_D_arr[i],dN)

a=[0.0,0.284,0.531,0.7]
#masks=np.array(a*3).reshape(4,3)

#chisqmat = np.ndarray(shape=(4,3,num,num))
#for i in range(4):
	#for j in range(3):
		#for k in range(20):
			#for l in range(20):
				#chisqmat[i,j,k,l]=dN_arr[i,j,k]*cov_inv_arr[i,k,l]*dN_arr[i,j,l]

#left=[0.05,0.92,0.05,0.92]
#bottom=[0.6,0.6,0.1,0.1]
#for j in range(3):
	#f=figure(figsize=(10,10))
	
	#for i in range(4):
		#ax=f.add_subplot(2,2,i+1)
		#im=ax.imshow(chisqmat[i,j],interpolation='nearest')
		##im=ax.imshow(chisqmat[i,j]/(1-a[i]),interpolation='nearest')
		#cax2 = f.add_axes([left[i], bottom[i], 0.02, 0.3])
		#f.colorbar(im,cax=cax2)	
		#ax.set_title(r'${\Delta}N_iC^{-1}_{ij}{\Delta}N_j (%s)$ mask %s'%(label_latex[j],a[i]))
	
	#figtext(.02, .02, 'Diagn(mask 0, 0.284, 0.531, 0.7) = %.1f,%.1f,%.1f,%.1f\nChisq(mask 0, 0.284, 0.531, 0.7) = %.1f,%.1f,%.1f,%.1f'%(Dchisq[0,j],Dchisq[1,j],Dchisq[2,j],Dchisq[3,j],Fchisq[0,j],Fchisq[1,j],Fchisq[2,j],Fchisq[3,j]))
	##savefig(plot_dir+'chisqelem%s_scaleArea.jpg'%(j))
	#savefig(plot_dir+'chisqelem%s.jpg'%(j))
	#close()
######################

###### plot single dN and chisq ########
#f=figure(figsize=(12,6))
#for j in range(3):
	#ax=f.add_subplot(2,2,j+1)
	#ax2=f.add_subplot(2,2,4)
	#for i in range(4):
		#ax.plot(kappa_arr,dN_arr[i,j]/(1-a[i]),label=a[i])
		#ax2.plot(kappa_arr,log10(diagonal(cov_inv_arr[i])))
	#ax.set_title(r'$%s$'%(label_latex[j]))
	#leg=ax.legend(loc=0, ncol=1, labelspacing=.2, prop={'size':10})
	#leg.get_frame().set_visible(False)
	
##show()
#savefig(plot_dir+'dN_inv.jpg')	
###############################################
		
########### plot error ellipse#####
#f=figure(figsize=(10,12))
#k=0
#for imask in (0,0.284,0.531,0.7):
	
	
	#fits = (all_fits[k].T)[1:]
	#fits_std = np.std(fits,axis=1)
	#for i in range(3):#
		#ax = f.add_subplot(4,3,i+1)
		
		#if i == 2:
			#j = 0
		#else:
			#j = i+1
		#x, y = fits[i], fits[j]
		#x0, y0 = average(x), average(y)
		#P = np.cov(x, y)
		#plotEllipse((x0,y0),P,edge=colors[k],label=str(imask))
		#ax.set_xlim(x0-2*fits_std[i],x0+3*fits_std[i])
		#ax.set_ylim(y0-2*fits_std[j],y0+3*fits_std[j])
		
		#ax.set_xlabel(r'$%s$'%(label_latex[i]))
		
		#leg=ax.legend(loc=0, ncol=1, labelspacing=.2, prop={'size':10})
		#leg.get_frame().set_visible(False)
		#ax.set_title('%sbins'%(bins))
		
		#ax2 = f.add_subplot(4,3,i+4)
		#plotEllipse((x0,y0),P,edge=colors[k],label=str(imask),mask=imask)
		#ax2.set_xlabel(r'$%s$'%(label_latex[i]))
		
		#leg2=ax2.legend(loc=0, ncol=1, labelspacing=.2, prop={'size':10})
		#leg2.get_frame().set_visible(False)
		#ax2.set_xlim(x0-2*fits_std[i],x0+3*fits_std[i])
		#ax2.set_ylim(y0-2*fits_std[j],y0+3*fits_std[j])
		#ax2.set_title('*sqrt(area)')
		
		#ax3 = f.add_subplot(4,3,i+7)
		#ax3.plot(kappa_arr,X_arr[k,i]*dp[i]/(1.0-imask),color=colors[k],label=str(imask))
		#ax3.set_xlabel(r'$\kappa$')
		#leg3=ax3.legend(loc=0, ncol=1, labelspacing=.2, prop={'size':10})
		#leg3.get_frame().set_visible(False)
		
		
		#ax4 = f.add_subplot(4,3,i+10)
		#ax4.plot(kappa_arr,X_arr[k,i]/(1.0-imask)/(X_arr[0,i])-1,color=colors[k],label=str(imask))
		#ax4.set_xlabel(r'$\kappa$')
		#leg4=ax4.legend(loc=0, ncol=1, labelspacing=.2, prop={'size':10})
		#leg4.get_frame().set_visible(False)
		#if i ==0:
			#ax.set_ylabel(r'$%s$'%(label_latex[j]))
			#ax2.set_ylabel(r'$%s$'%(label_latex[j]))
			#ax3.set_ylabel('dN')
			#ax4.set_ylabel('$\Delta$dN/dN')
	#k+=1
##show()
#figtext(.02, .02, description)
#savefig(plot_dir+dNtitle+'bins%s_kappa024.jpg'%(bins))
#close()
##################################


########### plot Fisher error ellipse#####
#f=figure(figsize=(10,12))
#k=0
#for imask in (0,0.284,0.531,0.7):
	
	
	#Fisher = FisherMat_arr[k]
	#fits_std = diagonal(sqrt(mat(Fisher).I))
	#F = mat(Fisher).I
	#for i in range(3):#
		#ax = f.add_subplot(4,3,i+1)
		
		#if i == 2:
			#j = 0
		#else:
			#j = i+1
		##x0, y0 = average(x), average(y)
		#x0, y0 = fidu_params[i],fidu_params[j]
		##P = np.cov(x, y)
		#P = F[[i,i,j,j],[i,j,i,j]].reshape(2,2)
		#plotEllipse((x0,y0),P,edge=colors[k],label=str(imask))
		#ax.set_xlim(x0-3*fits_std[i],x0+3*fits_std[i])
		#ax.set_ylim(y0-3*fits_std[j],y0+3*fits_std[j])
		
		#ax.set_xlabel(r'$%s$'%(label_latex[i]))
		
		#leg=ax.legend(loc=0, ncol=1, labelspacing=.2, prop={'size':10})
		#leg.get_frame().set_visible(False)
		#ax.set_title('%sbins (Fisher)'%(bins))
		
		#ax2 = f.add_subplot(4,3,i+4)
		#plotEllipse((x0,y0),P,edge=colors[k],label=str(imask),mask=imask)
		#ax2.set_xlabel(r'$%s$'%(label_latex[i]))
		
		#leg2=ax2.legend(loc=0, ncol=1, labelspacing=.2, prop={'size':10})
		#leg2.get_frame().set_visible(False)
		#ax2.set_xlim(x0-1.2*fits_std[i],x0+1.2*fits_std[i])
		#ax2.set_ylim(y0-1.2*fits_std[j],y0+1.2*fits_std[j])
		#ax2.set_title('axis*sqrt(area)')
		
		#ax3 = f.add_subplot(4,3,i+7)
		#ax3.plot(kappa_arr,X_arr[k,i]*dp[i]/(1.0-imask),color=colors[k],label=str(imask))
		#ax3.set_xlabel(r'$\kappa$')
		#leg3=ax3.legend(loc=0, ncol=1, labelspacing=.2, prop={'size':10})
		#leg3.get_frame().set_visible(False)
		
		
		#ax4 = f.add_subplot(4,3,i+10)
		#ax4.plot(kappa_arr,X_arr[k,i]/(1.0-imask)/(X_arr[0,i])-1,color=colors[k],label=str(imask))
		#ax4.set_xlabel(r'$\kappa$')
		#leg4=ax4.legend(loc=0, ncol=1, labelspacing=.2, prop={'size':10})
		#leg4.get_frame().set_visible(False)
		#if i ==0:
			#ax.set_ylabel(r'$%s$'%(label_latex[j]))
			#ax2.set_ylabel(r'$%s$'%(label_latex[j]))
			#ax3.set_ylabel('dN')
			#ax4.set_ylabel('$\Delta$dN/dN')
	#k+=1
##show()
#figtext(.02, .02, description)
#savefig(plot_dir+dNtitle+'bins%s_Fisher_equalbins.jpg'%(bins))
#close()
##################################

######### plot Fisher error ellipse, omsi combined#####

def prepareplot(fits_fidu,alpha=0.48):
	#fits_fidu = fits_fidu.T
	w_arr=fits_fidu.T[2]
	w_origin=average(w_arr)
	os_origin=average(fits_fidu.T[1])**alpha*average(fits_fidu.T[3])
	os_arr=fits_fidu.T[1]**alpha*fits_fidu.T[3]

	#w_os_mat=np.array([w_arr[~np.isnan(os_arr)],os_arr[~np.isnan(os_arr)]])
	P=np.cov(w_arr[~np.isnan(os_arr)],os_arr[~np.isnan(os_arr)])
	print len(w_arr[~np.isnan(os_arr)])
	return w_origin, os_origin, P

f=figure(figsize=(8,4))
k=0
for imask in (0,0.284,0.531,0.7):
	
	
	#Fisher = FisherMat_arr[k]
	#fits_std = diagonal(sqrt(mat(Fisher).I))
	#F = mat(Fisher).I
	
	fits = all_fits[k]#.T)#[1:]	
	ax = f.add_subplot(1,2,1)
	x, y, P=prepareplot(fits)
	
	x0, y0 = -1, 0.26**.48*0.798
	
	#P = F[1:,1:]
	U, s, Vh = svd(P)
	
	fits_std=[sqrt(s[0]),sqrt(s[1])]
	
	plotEllipse((x0,y0),P,edge=colors[k],label=str(imask))
	ax.set_xlim(x0-3*fits_std[0],x0+3*fits_std[0])
	ax.set_ylim(y0-3*fits_std[1],y0+3*fits_std[1])
	
	ax.set_xlabel('w')
	ax.set_ylabel(r'$%s^{0.48}%s$'%(label_latex[0],label_latex[2]))
	
	leg=ax.legend(loc=0, ncol=1, labelspacing=.2, prop={'size':10})
	leg.get_frame().set_visible(False)
	ax.set_title('%sbins'%(bins))
	
	ax2 = f.add_subplot(122)
	plotEllipse((x0,y0),P,edge=colors[k],label=str(imask),mask=imask)
	
	ax2.set_xlabel('w')
	ax2.set_ylabel(r'$%s^{0.48}%s$'%(label_latex[0],label_latex[2]))
	leg2=ax2.legend(loc=0, ncol=1, labelspacing=.2, prop={'size':10})
	leg2.get_frame().set_visible(False)
	ax2.set_xlim(x0-1.2*fits_std[0],x0+1.2*fits_std[0])
	ax2.set_ylim(y0-1.2*fits_std[1],y0+1.2*fits_std[1])
	ax2.set_title('axis*sqrt(area)')
	
	k+=1
#show()
#figtext(.02, .02, description)
savefig(plot_dir+dNtitle+'bins%s_Fisher_omsi_randommask.jpg'%(bins))
close()
################################
bins=200
getfit = 1
if getfit:
	k=0
	MB_arr=np.zeros(shape=(4,2))
	std_arr=np.zeros(shape=(4,2))
	for imask in a:
		bias_fits = genfromtxt(fit_dir+'bias_analytical_s0.8_g6.976_ngal30_bins%s_mask%s.ls'%(bins,imask))[:,1:].T
		true_fit = genfromtxt(fit_dir+'true_analytical_g6.976_ngal30_bins%s_mask%s.ls'%(bins,imask))[:,1:].T
		
		#MB = average(bias_fits,axis=0) - fidu_params
		#true_std = std(true_fit,axis=0)/sqrt(20000./12)
		
		MB_w = average(bias_fits[1]) - fidu_params[1]
		MB_omsi = average(bias_fits[0])**.48*average(bias_fits[2])) - 0.26**.48*.798
		
		w_std=std(true_fit[1])/sqrt(20000./12)
		omsi_std = std(true_fit[0]**.48*true_fit[2])/sqrt(20000./12)
		
		MB_arr[k,0]=MB_w
		MB_arr[k,1]=MB_omsi
		std_arr[k]=w_std,omsi_std
		
		#std_arr=true_std
		#print 'mask',imask,MB,(MB/true_std)
		k+=1