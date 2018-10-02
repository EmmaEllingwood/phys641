import numpy as np
from matplotlib import pyplot as plt
x=np.linspace(-5,5,1000)
y_true=5.0*x**3-3*x**2+9*x
y=y_true

ord=25
a=np.zeros([len(x),ord+1])
a[:,0]=1.0
for i in range(ord):
    a[:,i+1]=a[:,i]*x


#Classical
ata=np.dot(a.transpose(),a)
rhs=np.dot(a.transpose(),y) #Skipping the noise matrix A^T N^-1 d = A^T d
ata_inv=np.linalg.inv(ata)
fitp=np.dot(ata_inv,rhs)
q,r=np.linalg.qr(a)
r_inv=np.linalg.inv(r)
qt=q.transpose()
qtd=np.dot(qt,y)
y_pred_qr=np.dot(r_inv,qtd)
y_pred=np.dot(a,fitp)
ma=np.dot(a,y_pred_qr)


plt.scatter(x,y,label='Data Points',color='k',marker='.')
plt.plot(x,y_pred,label='Classical')
plt.plot(x,ma,label='QR')
plt.ylabel('y')
plt.xlabel('x')
plt.title('Question 1: Classical vs. QR Fit')
plt.legend()
plt.show()
