import numpy as np
from matplotlib import pyplot as plt


x=np.linspace(-1,1,10001)
y=np.exp(x)
plt.plot(x,y,'k',label='Exponential to be fit')
ord=50 #max order of chebyshev
p=np.zeros([len(x),ord+1])
p[:,0]=1.0
p[:,1]=x
for n in range(1,ord):
	p[:,n+1]=(2*x*p[:,n])-p[:,n-1]

q,r=np.linalg.qr(p)
r_inv=np.linalg.inv(r)
qt=q.transpose()
qtd=np.dot(qt,y)
y_qr=np.dot(r_inv,qtd)
y_qr[7:]=0   ####Comment out to do untruncated version of 50th order

#y_pred=np.dot(p,y_qr)
y_pred=np.dot(p,y_qr)
plt.plot(x,y_pred,label=str(ord))
#ee,vv=np.linalg.eig(ata)
#cond_pred=(2.0/(2*0+1))/(2.0/(2*ord+1)) #which simplifies to 2*ord+1
#print 'condition number is',ee.max()/ee.min(),' compared to expected ',cond_pred

rms=np.sqrt(sum((y_pred-y)**2)/len(y))
print 'Order: ',ord
print 'RMS:  ',rms
max_error=max(abs(y_pred-y))
print 'Max Error:  ',max_error
plt.plot(x,y_pred)
plt.legend()
#plt.show()
#plt.ion()
#plt.clf()
#plt.plot(x,p)
#plt.savefig('legendre_polys.png')
