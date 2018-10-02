import numpy as np
from matplotlib import pyplot as plt
from scipy.linalg import eigh,cholesky

ord=6 #Order of the matrix
##Creating the noise matrix
N=np.ones([ord,ord])
for i in range(ord):
	N[i,i]=2

print 'Noise Matrix:  '
print N

num = 100000

ddt_list=np.zeros([ord,ord])
for i in range(num):
	x=np.random.randn(ord)
	#eigval, eigvec = eigh(N)
	#c = np.dot(eigvec, np.diag(np.sqrt(eigval)))
	c = cholesky(N, lower=True)

	d = np.dot(c, x)
	ddt = np.outer(d.transpose(),d)
	ddt_list+=ddt
print 'Converged to Noise Matrix:'
print ddt_list.round(-2)/num #Rounded to make the matrix easier to look at

