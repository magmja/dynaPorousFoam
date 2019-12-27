'''
Based on DNV-RP-C205, Section3.2.2
page26
'''
import numpy as np
import os 
from numpy import pi

T= 6.0 # wave period
d= 1.0 # water depth
g=9.80665
a4=[0,0.66,0.445,-1.05,0.272]
for i in range(1,5):
    fomigas=a4[i]*pow(4*pi*pi*d/g/T/T,i)
fomigas+=1
L=T*np.sqrt(g*d*(fomigas/(1+4*pi*pi*d/g/T/T*fomigas)))
print('Wave length is '+str(L))