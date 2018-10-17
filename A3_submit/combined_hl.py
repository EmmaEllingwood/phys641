#This code just takes the SNR from each event for each detector calculated in the first part and just outputs the combined X^2 where I calculate it as SNR[H]^2+SNR[L]^2


from numpy import *
import matplotlib.pyplot as plt


def read_file(file_name):
    with open(file_name, 'r') as data:
        x=0
        y=0
        z=0
        for line in data:
            f = line.split()
            x=float(f[0])
            y=float(f[1])
            z=float(f[2])
    return x,y,z

#event='GW150914'
#event='LVT151012'
#event='GW151226'
event='GW170104'

hsignal,hnoise,hsnr=array(read_file(event+'/H/SNR.txt'))
lsignal,lnoise,lsnr=array(read_file(event+'/L/SNR.txt'))

combined=abs(hsnr**2+lsnr**2)

print combined

with open(event+'/combined.txt','w+') as file:
        print >> file, combined

