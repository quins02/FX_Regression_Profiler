import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import gaussian_kde

data = np.genfromtxt("tmp/CVA.dat",delimiter="\t")
data1 =	np.genfromtxt("tmp/CVAMean.dat",delimiter="\t")
data2 =	np.genfromtxt("tmp/CVAPMean.dat",delimiter="\t")

BF1 =	np.genfromtxt("tmp/BucketFit1.dat",delimiter="\t")
BF2 =	np.genfromtxt("tmp/BucketFit2.dat",delimiter="\t")
BF3 =	np.genfromtxt("tmp/BucketFit3.dat",delimiter="\t")
BF4 =	np.genfromtxt("tmp/BucketFit4.dat",delimiter="\t")
BF5 =	np.genfromtxt("tmp/BucketFit5.dat",delimiter="\t")
BF6 =	np.genfromtxt("tmp/BucketFit6.dat",delimiter="\t")
BF7 =	np.genfromtxt("tmp/BucketFit7.dat",delimiter="\t")
BF8 =	np.genfromtxt("tmp/BucketFit8.dat",delimiter="\t")
BF9 =	np.genfromtxt("tmp/BucketFit9.dat",delimiter="\t")
BF10 =	np.genfromtxt("tmp/BucketFit10.dat",delimiter="\t")
BF11 =	np.genfromtxt("tmp/BucketFit11.dat",delimiter="\t")
BF12 =	np.genfromtxt("tmp/BucketFit12.dat",delimiter="\t")

x = data[:,0]
y = data[:,1]

xy = np.vstack([x,y])
z = gaussian_kde(xy)(xy)

Mean = np.mean(data1[:,1])

plt.subplot(2,1,1)
plt.scatter(x,y,c=z,s=5, edgecolor='')
plt.plot(data1[:,0],data1[:,1])
plt.plot(data2[:,0],data2[:,1])
plt.xlabel('Time')
plt.ylabel('Exotic Value ($)')
plt.title('Risk Profile of Exotic Derivative')
plt.annotate(Mean,xy=(1,data[0,1]+0.1))
plt.grid(True)

plt.subplot(2,1,2)
plt.plot(BF1[:,0],BF1[:,1])
plt.plot(BF2[:,0],BF2[:,1])
plt.plot(BF3[:,0],BF3[:,1])
plt.plot(BF4[:,0],BF4[:,1])
plt.plot(BF5[:,0],BF5[:,1])
plt.plot(BF6[:,0],BF6[:,1])
plt.plot(BF7[:,0],BF7[:,1])
plt.plot(BF8[:,0],BF8[:,1])
plt.plot(BF9[:,0],BF9[:,1])
plt.plot(BF10[:,0],BF10[:,1])
plt.plot(BF11[:,0],BF11[:,1])
plt.plot(BF12[:,0],BF12[:,1])
plt.xlabel('FX Rate')
plt.ylabel('Exotic Value ($)')
plt.title('Value')
plt.xlim(0.5,1.5)
plt.ylim(0,1)
plt.grid(True)

plt.show()
