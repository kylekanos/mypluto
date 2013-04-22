import os
import sys
from numpy import *
from matplotlib.pyplot import *
import pyPLUTO as pp

plutodir = os.environ['PLUTO_DIR']
wdir = plutodir+'/Test_Problems/MHD/FARGO/Spherical_Disk/'
nlinf = pp.nlast_info(w_dir=wdir,datatype='float')

D = pp.pload(nlinf['nlast'],w_dir=wdir,datatype='float') # Loading the data into a pload object D.

I = pp.Image()

f1 = figure(figsize=[12,7],num=1)
ax1=f1.add_subplot(111)
I.pltSphData(D,w_dir=wdir,datatype='float',plvar='bx1',logvar=False,rphi=False,x3cut=96)
colorbar(orientation='horizontal')
ax1.set_xlabel(r'Radius')
ax1.set_ylabel(r'Height')
ax1.set_title(r'Magnetic field $B_{\rm x}$')

f2 = figure(figsize=[10,10],num=2)
ax2=f2.add_subplot(111)
I.pltSphData(D,w_dir=wdir,datatype='float',plvar='rho',logvar=True,rphi=True,x2cut=24)
colorbar(orientation='vertical')
ax2.set_xlabel(r'x')
ax2.set_ylabel(r'y')
ax2.set_title(r'Log $\rho$')

# Code to plot arrows. --> Spacing between the arrow can be adjusted by modifying the newdims tuple of conrid function.
T = pp.Tools()
newdims = 2*(20,)
R,Z,SphData = I.getSphData(D,w_dir=wdir,datatype='float',rphi=True,x2cut=24)
xcong = T.congrid(R,newdims,method='linear')
ycong = T.congrid(Z,newdims,method='linear')
vel1 = SphData['v1c']
vel2 = SphData['v3c']
                
xveccong = T.congrid(vel1,newdims,method='linear')
yveccong = T.congrid(vel2,newdims,method='linear')
normVp = sqrt(xveccong**2 + yveccong**2)
xveccong = xveccong/normVp
yveccong = yveccong/normVp
ax2.quiver(xcong, ycong, xveccong, yveccong,color='w')


f1.savefig('sphDisk_Bx.png')
f2.savefig('sphDisk_rho.png')
show()
