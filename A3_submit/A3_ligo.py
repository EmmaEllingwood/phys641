import numpy as np
import matplotlib.pyplot as plt
import h5py
import glob


#NOTE FOR PART b
#I will make a separate document that has a table of the SNR and possibly one that has all of the plots that I made since it may take a while to get all of them here

directory='LOSC_Event_tutorial/LOSC_Event_tutorial/'

#This is a simplified version of what was in the read_ligo.py file taking out any unnecessary parts
def read_template(filename):
    dataFile=h5py.File(filename,'r')
    template=dataFile['template']
    tl=template[1]
    return tl
def read_file(filename):
    dataFile=h5py.File(filename,'r')
    strain=dataFile['strain']['Strain'].value
    dataFile.close()
    return strain



#This is the smoothing function we used in glass and honestly I liked it as long I control it when I actually apply it which is what I am going to do later.
def smoothing_function(vec,fwhm):
    n=len(vec)
    x=np.arange(n)
    x[n/2:]=x[n/2:]-n
    sig=fwhm/np.sqrt(8*np.log(2))
    y=np.exp(-0.5*x**2/sig**2)
    y=y/np.sum(y)
    vecft=np.fft.rfft(vec)
    yft=np.fft.rfft(y)
    vec_smooth=np.fft.irfft(yft*vecft,n)
    return vec_smooth


plt.figure(1)
####Events
#To run a specific event, uncomment the event and template for each and then depending on whether it is Hanford or Livingston, uncomment out the file that starts with H-H1 or L-L1 depending on the case and the corresponding detector marker. This is just my system so that the aimges, and values are saved to the correct directories

##GW150914
#event='GW150914'
#template_name=directory+'GW150914_4_template.hdf5'
#file=directory+'H-H1_LOSC_4_V2-1126259446-32.hdf5'
#detector='H'
#file=directory+"L-L1_LOSC_4_V2-1126259446-32.hdf5"
#detector='L'



###LVT151012
#event='LVT151012'
#template_name=directory+"LVT151012_4_template.hdf5"
#file=directory+"H-H1_LOSC_4_V2-1128678884-32.hdf5"
#detector='H'
#file=directory+"L-L1_LOSC_4_V2-1128678884-32.hdf5"
#detector='L'



##GW151226
#event='GW151226'
#template_name=directory+"GW151226_4_template.hdf5"
#file=directory+"H-H1_LOSC_4_V2-1135136334-32.hdf5"
#detector='H'
#file=directory+"L-L1_LOSC_4_V2-1135136334-32.hdf5"
#detector='L'



##GW170104
event='GW170104'
template_name=directory+"GW170104_4_template.hdf5"
#file=directory+"H-H1_LOSC_4_V1-1167559920-32.hdf5"
#detector='H'
file=directory+"L-L1_LOSC_4_V1-1167559920-32.hdf5"
detector='L'


print 'reading file ',file
data=read_file(file)

template=read_template(template_name)


n=len(data)

x=np.linspace(-1,1,n)
#win=0.5*np.cos(np.pi*x)+0.5
#plt.plot(x,win)


#Part(a):   Trying a different window shape, this is a planck-taper window which by some research seems to be what they used in LIGO. I looked up the equation and put it here. I tried using the cosine we used in class, but compared to this window it gave a worse SNR overall. The first plot that appears (Figure 1) is the shape of the window.

def pt_window(e,N):
	win=[]
	for b in range(int(N)):
		if b>0 and b<e*(N-1):
			Z_plus=2*e*((1/(1+(((2*b)/(N-1))-1)))+(1/(1-(2*e)+(((2*b)/(N-1))-1))))
			y=1/(np.exp(Z_plus)+1)
			win.append(y)
		elif b>=e*(N-1) and b<=(1-e)*(N-1):
			y=1
			win.append(y)
		elif b>(1-e)*(N-1) and b<(N-1):
			Z_minus=2*e*((1/(1-(((2*b)/(N-1))-1)))+(1/(1-(2*e)-(((2*b)/(N-1))-1))))
			y=1/(np.exp(Z_minus)+1)
			win.append(y)
		else:
			y=0
			win.append(y)
	return np.array(win)


#Defining the window and plotting it just so I know the shape
win=pt_window(0.5,float(n))	
plt.figure(1)
plt.plot(x,win)

plt.savefig(event+'/'+detector+'/1window.png')

Fk=np.abs(np.fft.rfft(win*data))**2


#Smoothing Reasoning
# I know from our discussions that I can't just smooth simply and that I want to get rid of the times when the signal randomly goes low but I also want to severely downweight/get rid of the skinny high peaks which are not actually data.
#This is the smoothing function we used in class and honestly I liked it as long I control it when I actually apply it which is what I am going to do later.
def smoothing_function(vec,fwhm):
    n=len(vec)
    x=np.arange(n)
    x[n/2:]=x[n/2:]-n
    sig=fwhm/np.sqrt(8*np.log(2))
    y=np.exp(-0.5*x**2/sig**2)
    y=y/np.sum(y)
    vecft=np.fft.rfft(vec)
    yft=np.fft.rfft(y)
    vec_smooth=np.fft.irfft(yft*vecft,n)
    return vec_smooth
    
#So the way I am actually applying this is that at first I do a small smooth (N). Then I define a line where all of the strong thin peaks tend to pass. I keep track of the element of the array where this occurs so I can downweight them later because I still want to smooth more to get rid of that low random noise but I do not want the upper peaks to affect it. So I set it to something that is fairly close to the line so even in smoothing, it will not cause a huge shift in where the smooth should be. After I smooth a second time, then I take the weights and set the weight of the beginning messy part and the ending decreasing baseline part as we did in the class version. I also set the weights of the points that I identified before to zero because I did not want them to be represented at all in the final data set.

plt.figure(3)

N=smoothing_function(Fk,4)

plt.plot(np.log(Fk),label='Raw')
plt.plot(np.log(N),label='First Smooth')

index=[]
for i in range(len(N)):
	if np.log(N[i])>-83:
		N[i]=np.exp(-85)
		index.append(i)

plt.plot(np.log(N),label='Downscaling Peaks')

N=smoothing_function(N,15)


plt.plot(np.log(N),label='Second Smooth')
plt.legend()
Nmhalf=1/np.sqrt(N)
Nmhalf[:200]=0  #Cutting out the edges of the data like we did in class, I kept an eye on this while I was working through all the cases so I know these cuts are okay for all of the data sets.
Nmhalf[53000:]=0
for i in index:
	Nmhalf[i]=0

plt.savefig(event+'/'+detector+'/2Fk.png')

##Working with the data and whitening it
fourier_data=np.fft.rfft(data*win)
fourier_data_whitening=fourier_data*Nmhalf
data_whitening=np.fft.irfft(fourier_data_whitening,len(data))



##Same procedure but with the template
fourier_template=np.fft.rfft(template*win)
fourier_template_whitening=fourier_template*Nmhalf
template_whitening=np.fft.irfft(fourier_template_whitening,len(template))


##Matched filter procedure (part b)
fourier_matched_filter=fourier_data_whitening*np.conj(fourier_template_whitening)
matched_filter=np.fft.irfft(fourier_matched_filter,n)

#I am just shifting the signal so it is near the middle which I am pretty sure is what we did in class so it is not far on the left side. I just took the right side of the data and put it into the left
mflow=matched_filter[0:len(matched_filter)/2]
mfhigh=matched_filter[len(matched_filter)/2:len(matched_filter)]		
mfnew=list(mfhigh)+list(mflow)
mfnew=np.array(mfnew)

plt.figure(2)
plt.plot(mfnew)
plt.savefig(event+'/'+detector+'/3wave.png')

print max(matched_filter)


#Part c
#Yet again this is what we spoke about. I set the amplitude of the signal as the maximum of the absolute values
#For the error I looked close to the important signal as I could which I think I will adjust manually as my python has been acting up a bit and did not want to do something peak based.
#I then calculated the standard deviation right next to the signal which is why I set a low edge and high edge.

absmf=abs(matched_filter)
signal=max(absmf)

lowedge=int(8/15.*len(matched_filter))
highedge=int(9/15.*len(matched_filter))
plt.axvline(lowedge,color='r') #Just plotting the range where I am drawing the error from 
plt.axvline(highedge,color='r')
mfnoise=mfnew[lowedge:highedge]
stdmf=np.std(mfnew[lowedge:highedge])
print 'Signal: ',signal,' , Error:',stdmf,' , SNR:',signal/stdmf



#Part d: calculating the frequency, I am basically creating an array just called 'cumulative' which is just the cumulative weights and then finding the index when the weight equals half the weight of the total. 
#A* A really, but abs(A^2) works as well and is computationally a lot quicker for my computer

d=abs(fourier_template**2*Nmhalf**2)
half=sum(d)/2.
c=0
v_half=0
cumulative=np.zeros(len(Nmhalf))
previous=d[0]
cumulative[0]=previous
for x in range(1,len(Nmhalf)):
	cumulative[x]=previous+d[x]
	if cumulative[x-1]<half and cumulative[x]>half:
		c=cumulative[x-1]
		v_half=(x-1)/2/np.pi
	previous=cumulative[x]


#These are just the plots to check that I did the loops correctly. See the second page of my writeup for the plots
plt.figure(9)
plt.plot(d)
plt.axvline(v_half*2*np.pi,color='r')
plt.xlabel('w')
plt.ylabel('Cumulative $T^{2}/\sigma^{2}$')

plt.figure(10)
plt.plot(cumulative)
plt.axvline(v_half*2*np.pi,color='r')
plt.xlabel('w')
plt.ylabel('$T^{2}/\sigma^{2}$')


print 'Half weight frequency: ',v_half		
#Part d: So k is the wavenumber which should corresponds really to a spatial frequency so I am going to give it in these units because to me it makes more sense than to put it in a temporal frequency. Since I am doing this for multiple values I will just print this to the same file as the signal to noise ratio
#Part d: Note that this is a temporal frequency so k=1/lambda=frequency/c. So frequency = c*k

#Just printing the outputs to a file so that I don't have to run the code everytime to get the values I wanted, printing out the data set will be used in combined_hl.py to get the combined significance of the two events.
with open(event+'/'+detector+'/SNR.txt','w+') as file:
        print >> file, signal,stdmf,signal/stdmf,v_half


with open(event+'/'+detector+'/data_'+event+'_'+detector+'.dat','w+') as file:
    for f1 in mfnew:
        print >> file, f1
        

plt.show(block=True)


