s_arr = (0.0, 0.5, 1.0, 1.5, 2.0)
bins_arr =  (10, 20, 30, 40, 50, 60, 70, 80, 90, 100)

n=1
for s in s_arr:
	for bins in bins_arr:
		txt='''command: python /direct/astro+u/jia/magbias/mag-bias-astro.py %s %s'''%(s,bins)
		f=open('jobfile%s'%(n),'w')
		f.write(txt)
		f.close
		n+=1
		
for s in s_arr:
	for bins in bins_arr:
		print '''wq sub -c \'nohup python /direct/astro+u/jia/magbias/mag-bias-astro.py %s %s > nohup.out\' &'''%(s,bins)