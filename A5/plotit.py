import numpy as np
import matplotlib.pyplot as plt
import camb

order=['ombh2','omch2','H0','tau','As','ns']
maxi=10000
mini=1000

#chains_c=np.loadtxt('chains_c_newlong10000.txt') #chain from part c no source
#chains_c=np.loadtxt('chains_d_newlong10000.txt')  #chain from a_src=0.00075
chains_c=np.loadtxt('chains_d_newlong10000_point003src.txt') #chain from asrc=0.003


plt.figure(3)
means=np.zeros(6)
errors=np.zeros(6)
for x in range(6):
	plt.subplot(3,2,x+1)
	plt.plot(chains_c[mini:maxi,x+1])
	plt.title(order[x])
	values=chains_c[mini:maxi,x+1]
	means[x]=np.mean(chains_c[mini:maxi,x+1])
	errors[x]=np.std(chains_c[mini:maxi,x+1])


for xp in range(0,6,1):
	plt.figure(xp+10)
	a=1
	for yp in range(0,6,1):
		plt.subplot(2,3,a)
		#xp=4
		#yp=5
		if xp!=yp:
			plt.plot(chains_c[mini:maxi,xp],chains_c[mini:maxi,yp],'.')
			plt.xlabel(order[xp])
			plt.ylabel(order[yp])
			plt.tight_layout()
		a+=1

print order
print 'means=',means
print 'errors=',errors



##Importing the data file
data=np.loadtxt('wmap_tt_spectrum_9yr_v5.txt')
ell=data[:,0]
dat=data[:,1]  #l(l+1)/2pi C_l
error=data[:,2]  #From diagonal terms of Fisher matrix


#Using planck best fit from table 1 in the planck 2018 cosmology paper
#https://arxiv.org/pdf/1807.06209.pdf
h=0.6732
H0=h*100
omega_b=0.022383
omega_c=0.12011
tau=0.0543
As=np.exp(3.0448)/(10.**10)
ns=0.96605

source=0.003
#source=0.00075
#source=0

plt.figure(1)
dat=dat+(ell**2*source)
plt.plot(ell,dat,label='WMAP Data')


cos=np.asarray([omega_b,omega_c,H0,tau,As,ns])   #Original planck vlaues

cos=np.asarray([2.18467896e-02,2.58789051e-01,4.12424339e+01,2.43169404e-02,3.61875679e-09,1.17756180e+00]) ##src=0.003

#cos=np.asarray([2.23098702e-02, 1.46855814e-01, 6.01412879e+01, 9.04309095e-02, 2.68832248e-09, 1.02657516e+00]) ##src=0.00075

#cos=np.asarray([2.27115863e-02, 1.12801698e-01, 7.01844481e+01, 7.87862480e-02, 2.16943137e-09, 9.77874207e-01])  ##Part c...no source

#####Calculating the power spectrum using CAMB and the intput parameters
pars=camb.CAMBparams()
pars.set_cosmology(H0=cos[2],ombh2=cos[0],omch2=cos[1],tau=cos[3])
pars.InitPower.set_params(As=cos[4],ns=cos[5])
pars.set_for_lmax(1200, lens_potential_accuracy=0)

results=camb.get_results(pars)

powers=results.get_cmb_power_spectra(pars,CMB_unit='muK')
for name in powers: print name
totCl=powers['total']


plt.plot(totCl[:,0],label='CAMB Results from Planck Parameters')
plt.xlabel('Multipole Moment l')
plt.ylabel('l(l+1)$C_{l}/(2\pi)$') 
plt.legend(loc='lower left')
plt.title('Part (d): WMAP Data with CAMB Power Spectrum')

plt.show()


