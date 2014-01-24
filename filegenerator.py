import numpy as np

s_arr = (0.0, 0.5, 1.0, 1.5, 2.0)
bins_arr =  (10, 20, 30, 40, 50, 60, 70, 80, 90, 100)


for s in (0.1,0.4):#np.arange(-4,4,0.1):
	s='%.1f'%(s)
	for bins in (100,):#bins_arr:
		for ngal in (15,30,45):
			for imethod in ('spline','fw','bw'):
		#			for iseed in np.arange(0,100000,1000):
				f=open('jobfilenoise_s%s_bins%s_ngal%s_imethod%s'%(s,bins,ngal,imethod),'w')
				string = 'command: python /direct/astro+u/jia/magbias/mag-bias-astro-noise.py %s %s %s %s'%(s,bins,ngal,imethod)
				f.write(string)
				f.close()  

### clean up fit files, center, width of bias vs s ##
fit_noise_dir = '/Users/jia/Documents/weaklensing/magbias/fit_noise_seed73000/%s/'%(interpmethod)

bins=30
ngal=45
bias_fit_mat_file = fit_noise_dir+'bias_fit_mat_bins30_ngal%s.ls'%(ngal)
if os.path.isfile(bias_fit_mat_file):
	bias_fit_mat=np.genfromtxt(bias_fit_mat_file)
else:
	bias_fit_mat = np.ndarray(shape=(80,6))
	i=0
	for s in np.arange(-4,4,0.1):
		s='%.1f'%(s)
		f=np.genfromtxt(fit_noise_dir+'bias_fit_s%s_g6.976_ngal%s_bins%s.ls'%(s,ngal,bins))
		bias_fits=f[:,1:]
		Oms	=	stats.bayes_mvs(bias_fits[:,0])
		ws	=	stats.bayes_mvs(bias_fits[:,1])
		sis	=	stats.bayes_mvs(bias_fits[:,2])
		iOm	=	Oms[0][0]
		iOm_wid	=	0.5*(Oms[2][1][0]+Oms[2][1][1])
		iw	=	ws[0][0]
		iw_wid	=	0.5*(ws[2][1][0]+ws[2][1][1])
		isi	=	sis[0][0]
		isi_wid	=	0.5*(sis[2][1][0]+sis[2][1][1])	
		print iOm,iOm_wid,iw,iw_wid,isi,isi_wid
		bias_fit_mat[i]=np.array([iOm,iw,isi,iOm_wid,iw_wid,isi_wid])
		i+=1
		
	savetxt(bias_fit_mat_file,bias_fit_mat)