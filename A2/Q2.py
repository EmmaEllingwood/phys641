import numpy as np
from matplotlib import pyplot as plt


x=np.linspace(-1,1,10001)
y=np.exp(x)
plt.plot(x,y,'k',label='Exponential to be fit')
ord=2 #max order of chebyshev
for ord in [2,3,4,10]:
	p=np.zeros([len(x),ord+1])
	p[:,0]=1.0
	p[:,1]=x
	for n in range(1,ord):
    		p[:,n+1]=(2*n+1)*x*p[:,n]-p[:,n-1]
	atd=np.dot(p.transpose(),y)
	ata=np.dot(p.transpose(),p)
	m=np.dot(np.linalg.inv(ata),atd)
	print len(m)
	y_pred=np.dot(p,m)
	plt.plot(x,y_pred,label=str(ord))
	ee,vv=np.linalg.eig(ata)
	cond_pred=(2.0/(2*0+1))/(2.0/(2*ord+1)) #which simplifies to 2*ord+1
	print 'condition number is',ee.max()/ee.min(),' compared to expected ',cond_pred


plt.legend()
plt.show()
#plt.ion()
#plt.clf()
#plt.plot(x,p)
#plt.savefig('legendre_polys.png')






