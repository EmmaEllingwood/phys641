import numpy as np
import matplotlib.pyplot as plt
from scipy.linalg import eigh

ord=1000
a=0.1
sig=5
z=1
for a in [0.1,0.5,0.9]:
	for sig in [5,50,500]:
		N=np.zeros([ord,ord])
		for i in range(ord):
			for j in range(ord):
				if i==j:
					N[i,j]=a*np.exp(-0.5*(i-j)**2/sig**2)+(1-a)
				else:
					N[i,j]=a*np.exp(-0.5*(i-j)**2/sig**2)

		#print N
		x=np.array(range(1000))
		x0=np.mean(x)
		amp_true=1.0
		sigt=50
		template=np.exp(-0.5*(x-x0)**2/sigt**2)


		#Create correlated noise signal in data
		r=np.random.randn(ord)
		eigval, eigvec = eigh(N)
		c = np.dot(eigvec, np.diag(np.sqrt(eigval)))
		#c = np.linalg.cholesky(N, lower=True)
		print r[:3]
		y_correlated = np.dot(c, r)		
		dat=template+y_correlated
		plt.subplot(3,3,z)
		plt.plot(x,dat,'k',label='With Signal',alpha=0.8)
		
		Ninv=np.linalg.inv(N)
		
		ninvd=np.dot(Ninv,dat)
		atninvd=np.dot(template.transpose(),ninvd)
			
		d=(np.dot(Ninv,template)) #A^T N^(-1) A
		denom=np.dot(template.transpose(),d)
		###If I want to use a matched filter
		#amp=np.zeros(len(x))
		#for i in range(len(x)):
		#    template_shifted=np.exp(-0.5*(x-x[i])**2/sig**2)
		#    rhs=np.dot(template_shifted,np.dot(Ninv,dat))
		#    amp[i]=rhs/denom
		#plt.plot(x,amp,label='matched')		
		
		y_pred=np.dot(template,atninvd/denom)
		print 'a=',a, '  ,    sig=',sig, '     Error = %.3f'%(1/np.sqrt(denom))
		print '---------------'		
		
		#For part b
		plt.plot(x,y_correlated,'r',label='Without Signal',alpha=0.3)
			
		####Plotting option
		#plt.plot(x,y_pred,label='Fit')
		plt.title('a='+str(a)+' , sigma='+str(sig))
		plt.ylabel('Amplitude')
		#plt.xlabel('X')
		#plt.plot(x,template,label='Template')
		#plt.legend()
		
		z+=1

plt.show()


####Answer to part b
#Large 'a' means that the noise is mostly correlated and the amplitude of the noise will be larger, this results in larger error bars as 'a' increases. When sigma is very small, especially compared to the signal width, the variations are very pronounced and they vary on shorter scales so with a larger width. Similar to when sigma is high compared to the signal width because the deviations are on much longer scales so it is still possible to see the signal. When sigma is on the same scale as the signal width because it is very hard to distinguish any deviation as being part of a signal or just from the noise.

#The a=0.9, sig=50 combination has the worst error bar. The sigma is on the scale of the signal width and it is very correlated which makes it hard to distinguish the signal. The a=0.9,sig=500 has the best error bar, high correlation and sigma makes a smoothish noise, makes it easier to see the signal

#Noise with sigma=50 is the type that I should worry about more, specifically those with higher correlation because the variations in the noise is on the same scale as the width of the signal added in. So the variations are able to hide the signal in a way, to where it is less obvious. The high correlation also contributes to this to make the wiggles in the noise more pronounced.


