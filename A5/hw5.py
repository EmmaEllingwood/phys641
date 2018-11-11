import numpy as np
import matplotlib.pyplot as plt
import camb
from camb import model,initialpower
import time

####PART A

###Importing the data file
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

plt.figure(1)
dat=dat+(ell**2*source)
plt.plot(ell,dat,label='WMAP Data')

#####Calculating the power spectrum using CAMB and the intput parameters
pars=camb.CAMBparams()
pars.set_cosmology(H0=H0,ombh2=omega_b,omch2=omega_c,tau=tau)
pars.InitPower.set_params(As=As,ns=ns)
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

####Looking at the plot, the curve from CAMB and the Planck values seems to fit well to the WMAP data, there are no significant deviations, it is just a smoother curve showing the same general trend as the WMAP data.

order=['ombh2','omch2','H0','tau','As','ns']
######PART B
###This is just the cholesky calculation from my longest data run I did during (b), I run it again in (c) when calculating the correlation length
chains_c=np.loadtxt('chains_c_10000events.txt')
chains_c=np.asarray(chains_c)
chains_cor=chains_c[:,1:].copy()
for i in range(chains_cor.shape[1]):
	chains_cor[:,i]=(chains_cor[:,i]-chains_cor[:,i].mean())
cor_matrix=np.dot(chains_cor.transpose(),chains_cor)/chains_cor.shape[0]
g=np.linalg.cholesky(cor_matrix)

t0=time.time()

bad_chisq=1e100




def update(parameters,cosmo):
    p2=parameters.copy()
    p2.set_cosmology(ombh2=cosmo[0],omch2=cosmo[1],H0=cosmo[2],tau=cosmo[3])
    p2.InitPower.set_params(As=cosmo[4],ns=cosmo[5])
    return p2

def calculate_chisq(cosmo,dat,parameters):
    try:
        pars=update(parameters,cosmo)
        results=camb.get_results(pars)
        power=results.get_cmb_power_spectra(pars,CMB_unit='muK')['total']
    except:
        return bad_chisq
    inds=np.asarray(dat[:,0],dtype='int')
    predicted=power[inds,0]
    chisq=np.sum( (predicted-dat[:,1])**2/dat[:,2]**2)
    return chisq

pars=camb.CAMBparams()


cosmology=np.asarray([omega_b,omega_c,H0,tau,As,ns])

#I really had no sense of what a reasonable guess on the error would be to start with so I just started with the values set in the example from class and then only changed one of the parameters (see commented out parts of the code) and then ran for 300 samples. I saved the output standard deviation of the parameter as the new error. I did this separately for each of the parameters then ran for a longer stretch of time using the errors as my step size and doing about 10000 iterations.
#errs=np.asarray([2.4e-4, 7.34e-4,0.7,1.8e-3,7e-12,0.005])/2  #Feed in error when just starting (b)
data=np.loadtxt('WMAP_data.txt')
data[:,1]=data[:,1]+(data[:,0]**2*source)


chisq=calculate_chisq(cosmology,data,pars)
#errs=np.asarray([0.00020133384342450284/2., 0.0005506179258220545/4., 1.3035964051991118*2,0.0009772536508284595/2.,7.739417438688322e-12*2,0.0062575429436165265/4.]) #Rescaled step size at the end of (b)


accept_number=0
##changing_value=5
samples=10000
chains=np.zeros([samples,1+len(cosmology)])

#f=open('aafile.txt','w')

for x in range(samples):
    new_cosmo=cosmology+np.dot(g,np.random.randn(g.shape[1]))#np.random.randn(len(errs))*errs

    #For b, when I did the running of the chain with all other parameters fixed I just changed which value was changing and I replaced the live above with the two lines below. I then took the standard deviation of what came out and used that as the new error estimate. When I was only chaing one parameter at a time I found each error/stepsize and the I did a couple more runs with those values changing the overall scaling which is what is given at errs a few lines above where each is multiplied or divided by 2 or 4 or kept constant. This was part of the attempt to get the acceptance rate to be sensible. With these values I got to about 0.243 acceptance rate.

    #new_cosmo=cosmology
    #new_cosmo[changing_value]=cosmology[changing_value]+np.random.randn()*errs[changing_value]
    new_chisq=calculate_chisq(new_cosmo,data,pars)
    accept=False
    if new_chisq<bad_chisq and new_cosmo[3]>0:
        if np.random.rand()<np.exp(-0.5*(new_chisq-chisq)):
            accept=True
    print x, ' chisq = ',new_chisq, ' and the step status is ',accept
    print new_cosmo
    if accept:
        chisq=new_chisq
        cosmology=new_cosmo
        accept_number+=1
    chains[x,0]=chisq
    chains[x,1:]=cosmology
    line=np.zeros(7)
    line[0]=chisq
    line[1:]=cosmology
    np.savetxt('running_testd.txt',chains)  ###I tried to do it the way you did, but whenever I used the types of values that I had I always got the error 'ValueError: Wrong number of columns at line 2', it resolves itself only if I get rid of a lot of decimals or specifically when I get rid of the value which is at e-09. Either way I switched to using this much less elegant method although it seems to run just as quickly even through it is way less efficient. I am annoyed I have to do it this way instead of writing it line by line.
#    my_str=np.array2string(chains[x,:])
#    my_str=my_str[1:-1]+'\n'
#    print my_str
#    f.write(my_str)
#    f.flush()

#f.close()    

print 'Acceptance Rate: ',accept_number/float(samples)


np.savetxt('chains_d_newlong10000_point003src.txt',chains)
plt.figure(9)
for x in range(len(cosmology)):
	plt.subplot(3,2,x+1)
	plt.plot(chains[:,x+1])
	plt.title(order[x])
t1=time.time()
print 'Time Elapsed: ',t1-t0



####plt.show()
######Save one line at a time, flush between writing, keep file open
######PART C 
####Addition to part b code for correlation, in this case I chose to read from files so I didn't have to calculate anew everytime I wanted to test my code
#chains_c=np.loadtxt('chains_c_cortesthalf.txt')
##plt.figure(3)
##for x in range(6):
##	plt.subplot(3,2,x+1)
##	plt.plot(chains_c[:,x+1])
##	plt.title(order[x])
##chains_c=np.asarray(chains_c)
##plt.figure(4)
##ft1=np.fft.rfft(chains_c[:,1])
##plt.plot(ft1)
##plt.yscale('log')
##chains_cor=chains_c[:,1:].copy()
##for i in range(chains_cor.shape[1]):
##	chains_cor[:,i]=(chains_cor[:,i]-chains_cor[:,i].mean())
##cor_matrix=np.dot(chains_cor.transpose(),chains_cor)/chains_cor.shape[0]
##g=np.linalg.cholesky(cor_matrix)
##c=g*np.random.randn(g.shape[0])
##This is the same as what I copied to the beginning of (b) now, but below from any data set I want now I can calculate the correlation length and the mean and error of the data set.

##print cor_matrix
##print cor_matrix.shape[0]


##Calculating correlation length and finding the mean and error
#for i in range(1,chains_c.shape[1]):
#	dat=chains_c[:,i].copy()
#	meanval=np.mean(dat)
#	stdval=np.std(dat)
#	dat=dat-meanval
#	datft=np.fft.rfft(dat)
#	correlation=np.fft.irfft(datft*np.conj(datft))
#	length=np.min(np.where(correlation<0))
#	print 'The correlation for ',order[i-1],' is ', length, 'with about ',len(dat)/length,' independent samples'
#	print order[i-1],' mean=',meanval,' error=',stdval
#	
#	
##Just plotting each parameter as it changes during the chain	
#maxi=2000
##chains_c=np.loadtxt('running_test.txt')
#chains_c=np.loadtxt('twothousand.txt')
#plt.figure(3)
#for x in range(6):
#	plt.subplot(3,2,x+1)
#	plt.plot(chains_c[:maxi,x+1])
#	plt.title(order[x])

##Producing plots of one parameter vs. the others to see if there are any correlations
#for xp in range(0,6,1):
#	plt.figure(xp+10)
#	a=1
#	for yp in range(0,6,1):
#		plt.subplot(2,3,a)
#		#xp=4
#		#yp=5
#		if xp!=yp:
#			plt.plot(chains_c[:maxi,xp],chains_c[:maxi,yp],'.')
#			plt.xlabel(order[xp])
#			plt.ylabel(order[yp])
#			plt.tight_layout()
#		a+=1
#plt.show()

####PART D
#This section is mainly implemented in earlier parts of the code and the explanation is in the write up
