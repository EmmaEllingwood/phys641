import numpy as np
import matplotlib.pyplot as plt
import healpy

######Question 1

a=np.loadtxt('example_ps.txt')
data = a[:2501,0]
#plt.plot(data)
#plt.xscale('log')
#plt.title('l(l+1)/(2 pi) Cl')

ll=np.zeros(len(data))
for x in range(len(data)):
	ll[x]=x*(x+1)
ll[0]=1


cl=2*np.pi*data/ll

twol=np.zeros(len(data))
for x in range(len(data)):
	twol[x]=2*x+1
print 'The total variance is the sum of the individual Cl which leaves the variance as ', sum(cl*twol)/(2*np.pi)/2.

##plt.figure(10)
##plt.loglog(cl)
##plt.title('Cl')

##plt.figure(15)
##plt.hist(cl,bins=np.linspace(0,5,100))
##plt.yscale('log')


##Part (b)

#lmax=2500

#nside=1024
#nalm=int((lmax+1)*(lmax+2)/2.)
#	

##print 'Var my data',np.var(alm4)
#####Looking at this variance, seems to be off by approximately 4pi
#i=0
#alms=np.zeros(nalm,dtype='complex')
#for m in range(lmax+1):
#	for l in range(m,lmax+1):
#		if m>0:
#			alms[i]=np.sqrt(cl[l])*np.random.randn()+(1J)*np.sqrt(cl[l])*np.random.randn()
#		else:
#			alms[i]=np.sqrt(cl[l])*np.random.randn()
#		i+=1
#print i
#maps=healpy.alm2map(alms,nside)
#print 'map var ',np.var(maps)
#healpy.mollview(maps)
#plt.title('My Attempt')
#alm_adj=alms
#alm_adj[lmax+1:]=alm_adj[lmax+1:]/np.sqrt(2)
#map_adj=healpy.alm2map(alm_adj,nside)
#new_cl_adj=healpy.anafast(map_adj)
#print 'map adj var ',np.var(map_adj)

#new_cls=healpy.anafast(maps,lmax=lmax)

#plt.figure(15)
#plt.loglog(new_cls,label='Output Cl')
#plt.loglog(cl,label='Input Cl')
#plt.loglog(new_cl_adj,label='Adjusted Output Cl')
#plt.title('My attempt return to Cl')
##plt.xscale('log')
#plt.legend()



##Part (d)

#print 'PART D'




#alm=healpy.synalm(cl,new=True)

#map=healpy.alm2map(alm,nside)
#print 'part d map var ',np.var(map)
#healpy.mollview(map)
#plt.title('Using synalm')

#new_cl=healpy.anafast(map,lmax=lmax)



#plt.figure(4)
#plt.loglog(new_cl,label='Output Cl')
#plt.loglog(cl,label='Input Cl')
#plt.title('synalm map return to Cl using anafast')
##plt.xscale('log')
#plt.legend()

#####QUESTION 2
patchsize_deg=20.
pixelsize_deg=0.02

npix=int(patchsize_deg/pixelsize_deg)
kmax=npix


lmax=2500

nalm=int((lmax+1)*(lmax+2)/2.)
print nalm
n=2*npix+1

x=np.arange(n)
x[n/2:]=x[n/2:]-n
dat=np.random.randn(n,n)
datft=np.fft.fft2(dat)
kx=np.repeat([x],n,axis=0)
ky=np.transpose(kx)
k=np.sqrt(kx**2+ky**2)


k[0,0]=0.0

P=np.zeros([n,n])  #P is the power spectrum in terms of k
for x in k:
	k=np.round(k)


for x in range(len(k)):
	for y in range(len(k)):
		l_value=int(k[x,y]*18)
		if l_value<=lmax:
			P[x,y]=cl[l_value]

nalm2=(18001*18002)/2
f=np.sqrt(nalm2/(2*k+1))
dat_back=np.fft.ifft2(datft*np.sqrt(P)*f)
dat_back=np.real(dat_back)
print 'The variance of the flat-sky map is ',np.var(dat_back)
plt.imshow(dat_back)
plt.colorbar()



plt.show()


#####Question 3
##Part a
##The equation for the Planck function is ( 2*h*v**3/(c**2)*(exp((h*v)/(k*T))-1)**(-1)) and the units are erg s**(-1) cm**(-2) Hz**(-1) sr**(-1) so just multiplying the equation by the area lambda**2 and  1 sr and using T=2.725 gives the ergs per second per Hz requested. 
##The one thing I am not sure about is what is meant by lambda**2, the assumption that I am making is that it is the central wavelength corresponding to the center of the band 150 GHz so that is what I will use for the rest
##This gives a full expression as ( 2*h*v**3/(c**2)*(exp((h*v)/(k*T))-1)**(-1))*(c/150GHz)**2 where the units are now ergs s**(-1) Hz**(-1)

#import scipy.integrate as integrate
#h=6.626068e-27 #erg sec
#k=1.38066e-16 #erg/K
#c=2.997925e10 #cm/sec
##Just defining a function so I could put in any frequency range. To convert ergs to number of photons I just divided E=hv because at every frequency where I calculate, if the energy output is the same that corresponds to different 
#def I(wavelow,wavehigh,T): #wavelength in cm
#	I_list=[]
#	
#	
#	for x in np.linspace(wavelow,wavehigh,1000):
#		B=2*h*x**3/(c**2)/(np.exp((h*x)/(k*T))-1)*(c**2/float(150.e9)**2)
#		I_list.append(B)
#	plt.figure(21)
#	plt.plot(np.linspace(wavelow,wavehigh,1000),I_list)
#	plt.title('Planck Function for T='+str(T))
#	plt.xlabel('Frequency (Hz)')
#	plt.ylabel('Intensity ergs s**(-1) Hz**(-1)')

#I(0,400e9,2.725)
#def I_3(wavelow,wavehigh,T): #wavelength in cm
#	return integrate.quad(lambda x:( 2*h*x**3/(c**2)*(np.exp((h*x)/(k*T))-1)**(-1)*(h*x)**(-1)),wavelow,wavehigh)


#photon_rate=int(I_3(135e9,165e9,2.725)[0]*(c**2)/float(150e9**2))
#print 'The number of photons per second in this band is  ',photon_rate,' or roughly ',float(round(photon_rate/1e9*10))/10. ,' billion photons per second'

#plt.show()
##If I have 4.6 billion photons per second like this equation gives, but I am sampling at 150 billion samples per second then on average I would need to sample 30 times to see a photon. This is not such a drastic difference that I can't say for example that I could see two photons in a single sample which would make it more continuous. Since there is less photons per second than the rate I am sampling at this means it is a bit more like shot-noise but it seems a bit on the edge so it does not seem to be clearly one or the other.



##Part c: Since I want the temperature representation for noise, instead of using delta n / n
##delta T / T = 1/sqrt(nt) with t=1s, delta T = 40.2 uK
##This is better than the best Planck detectors, but the detector we calculated is really for a perfect detector where it is seeing exactly the CMB temperature and there is nothing either raising that temperature and no other sources of noise. There is always ways to improve detectors, but just because we are not in a situation where we could have exactly this ideal case, it could be possible to get deeper maps, but it would probably not reach the level of this ideal detector.
