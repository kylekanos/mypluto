# -*- coding: utf-8 -*-
import os
import sys
import struct
import numpy as np
import linecache
import scipy as S
import scipy.ndimage
from scipy import integrate
import scipy.interpolate
from matplotlib.pyplot import *
from matplotlib.mlab import *
from scipy.interpolate import UnivariateSpline
import time

def curdir():
	""" Get the current working directory.
	"""
        curdir = os.getcwd()+'/'
        return curdir

def get_nstepstr(ns):
    """ Convert the float input *ns* into a string that would match the data file name.

    **Inputs**:

      ns -- Integer number that represents the time step number. E.g., The ns for data.0001.dbl is 1.\n

    **Outputs**:

      Returns the string that would be used to complete the data file name. E.g., for data.0001.dbl, ns = 1 and pyPLUTO.get_nstepstr(1) returns '0001'
    """
    nstepstr = str(ns)
    while len(nstepstr) < 4:
        nstepstr= '0'+nstepstr
    return nstepstr

def nlast_info(w_dir=None,datatype=None):
	""" Prints the information of the last step of the simulation as obtained from dbl.out or flt.out

	**Inputs**:
	
	  w_dir -- path to the directory which has the dbl.out(or flt.out) and the data\n
	  datatype -- If the data is of 'float' type then datatype = 'float' else by default the datatype is set to 'double'.

        **Outputs**:
	
	  This function returns a dictionary with following keywords - \n

	  nlast -- The ns for the last file saved.\n
	  time -- The simulation time for the last file saved.\n
	  dt -- The time step dt for the last file. \n
	  Nstep -- The Nstep value for the last file saved.


	**Usage**:
	
	  In case the data is 'float'.

	  ``wdir = /path/to/data/directory``\n
	  ``import pyPLUTO as pp``\n
	  ``A = pp.nlast_info(w_dir=wdir,datatype='float')``
	

	
	
	"""
	if w_dir is None: w_dir=curdir()
	if datatype == 'float':
		fname_v = w_dir+"flt.out"
	else:
		fname_v = w_dir+"dbl.out"
	last_line = file(fname_v,"r").readlines()[-1].split()
	nlast = int(last_line[0])
	SimTime =  float(last_line[1])
	Dt = float(last_line[2])
	Nstep = int(last_line[3])
	    
	print "------------TIME INFORMATION--------------"
	print 'nlast =',nlast
	print 'time  =',SimTime
	print 'dt    =', Dt
	print 'Nstep =',Nstep
	print "-------------------------------------------"
	    
	return {'nlast':nlast,'time':SimTime,'dt':Dt,'Nstep':Nstep}
    



class pload(object):
    """
    This Class has all the routines loading the data from the
    binary files output from PLUTO Simulations. Assign an object
    when the data is loaded for some *ns*.

    **Inputs:**

      ns -- Time step of the data, data.ns.dbl or data.ns.flt\n
      w_dir -- path to the directory which has the dbl.out(or flt.out) and the data\n
      datatype -- If the data is of 'float' type then datatype = 'float' else by default the datatype is set to 'double'.
    

    **Outputs:**

      A pyPLUTO.pload object  having all the relevant information of the corresponding data file 

    **Usage**:

      ``import pyPLUTO as pp``\n
      ``wdir = '/path/to/the data files/'``\n
      ``D = pp.pload(1,w_dir=wdir)``\n

      Now D is the pyPLUTO.pload object having all the relevant information
      of the corresponding data file - data.0001.dbl.


    """
    def get_varinfo(self,datatype):
	""" This method reads the dbl.out (or flt.out) and stores the information in a dictionary.

	**Inputs:**

	  datatype -- If the data is of 'float' type then datatype = 'float' else by default the datatype is set to 'double'.

	**Outputs:**
	
	  A dictionary with the following keywords - 
	
	  fltype -- returns the filetype of storing data (single_file or multiple_file)\n
	  nvar -- Number of variables\n
	  allvars -- A list of variables names. Each of the variable name will be the attributes to the pyPLUTO.pload object.\n
	
	"""
	if datatype == 'float':
		fname_v = self.wdir+"flt.out"
	else:
		fname_v = self.wdir+"dbl.out"

	f_var = open(fname_v)
        lnum_v = len(f_var.readlines())
        Var_info = linecache.getline(fname_v,1).split()
        fltype = Var_info[4]
        nvar = len(Var_info[6:])
        allvars=[]
        for i in Var_info[6:]:
            allvars.append(i)
        
            if i == 'bx1s':
                nvar = nvar - 1
                del allvars[-1]
            
            if i == 'bx2s':
                nvar = nvar - 1
                del allvars[-1]

	    if i == 'bx3s':
	        nvar = nvar - 1
                del allvars[-1]
                
        f_var.close()
        return {'fltype':fltype, 'nvar':nvar, 'allvars':allvars}

    def geometry(self):
	""" This method has the geometry information of the problem considered.

	"""
	fname_d = self.wdir+"grid.out"
        f_def = open(fname_d)
	
        lnum_d = len(f_def.readlines())
        Geo_info = linecache.getline(fname_d,4).split()
        f_def.close()
        print "GEOMETRY >> " + Geo_info[2]

    def time_info(self,datatype):
	"""
	This method returns a dictionary that has the time information for the
	step ns.

        **Inputs:**

          datatype -- If the data is of 'float' type then datatype = 'float' else by default the datatype is set to 'double'.

        **Outputs:**

          A dictionary which has a following keys -
	
          time -- Gets the simulation time at step ns.\n
          dt -- Get the time step dt for step ns.\n
          Nstep -- Get the value of nstep for the step ns.\n

	"""
	ns = self.nStep
	if datatype == 'float':
		fname_v = self.wdir+"flt.out"
	else:
		fname_v = self.wdir+"dbl.out"
	
        f_var = open(fname_v)
	tlist = []
	for line in f_var.readlines():
		tlist.append(line.split())
		
	
	SimTime =  float(tlist[ns][1])
	Dt = float(tlist[ns][2])
	Nstep = int(tlist[ns][3])

	return {'time':SimTime,'dt':Dt,'Nstep':Nstep}

	
	    

    def grid(self):
	"""
	This method returns the necessary grid information in form a dictionary after reading the grid.out.

        **Outputs:**

          This function outputs a dictionary with following keys. 
	
	  n1 -- number of grid cells in x1 direction\n
	  n2 -- number of grid cells in x2 direction\n
	  n3 -- number of grid cells in x3 direction\n
	  x1 -- Array x1\n
	  x2 - Array x2\n
	  x3 - Array x3\n
	  dx1 - Array dx1\n
	  dx2 - Array dx2\n
	  dx3 - Array dx3\n

        **Usage:**

          ``import pyPLUTO as pp``\n
          ``D = pp.pload(1)``\n
          ``G = D.grid()``

          Now, G is a dictionary such that G['x1'] will give the x-array. **NOTE** - Even x1 is set as attribute for the pyPLUTO.pload object i.e., the x-array can be obtained as D.x1. Infact all the keys of this routine are set as attributes for the pyPLUTO.pload object D.
	
	"""
	    
        fname_g = self.wdir+"grid.out"
        f_grid=open(fname_g)
        lnum_g = len(f_grid.readlines())
	f_grid.close()
	f_grid=open(fname_g)
	headlines = []
	for line in f_grid.readlines():
		if line.split()[0] == '#' and len(line.split()) > 1:
			headlines.append(line.split()[1:])
	
	offset = len(headlines)+1
	
        n1 = linecache.getline(fname_g,1+offset)
        n2 = linecache.getline(fname_g,int(n1)+2+offset)
        n3 = linecache.getline(fname_g,int(n1)+int(n2)+3+offset)
	
        x1=[]
        x2=[]
        x3=[]
        dx1=[]
        dx2=[]
        dx3=[]
        for i in range(2,int(n1)+2):
            A = linecache.getline(fname_g,i+offset).split()
            x1.append(0.5*(float(A[1])+float(A[2])))
            dx1.append(float(A[2])-float(A[1]))
        
        x1 = np.asarray(x1)
        dx1 = np.asarray(dx1)
	

        for j in range(3+int(n1),int(n1)+int(n2)+3):
            B = linecache.getline(fname_g,j+offset).split()
            x2.append(0.5*(float(B[1])+float(B[2])))
            dx2.append(float(B[2])-float(B[1]))
        
        x2 = np.asarray(x2)
        dx2 = np.asarray(dx2)
	

        for k in range(4+int(n1)+int(n2),lnum_g+1):
            C = linecache.getline(fname_g,k+offset).split()
	    if len(C)>1:
		x3.append(0.5*(float(C[1])+float(C[2])))
		dx3.append(float(C[2])-float(C[1]))

        x3 = np.asarray(x3)
        dx3 = np.asarray(dx3)

        f_grid.close()
	
        grid_dict={'n1':int(n1),'n2':int(n2),'n3':int(n3),'x1':x1,'x2':x2,'x3':x3,'dx1':dx1,'dx2':dx2,'dx3':dx3}

	return grid_dict


    def data(self,datatype):
	""" This method loads the data from the file name "data.ns.dbl"(in case of single_file) or "varname.ns.dbl"(in case of multiple files).

        **Inputs:**

          datatype -- If the data is of 'float' type then datatype = 'float' else by default the datatype is set to 'double'.

        **Outputs:**

	  It returns a dictionary with all the variable names as the keys. **NOTE** All these variable names are set as attributes to the pyPLUTO.pload object like the grid() function. 
  
	"""
	print "Working Dir : %s" % (self.wdir)
        grid_dict = self.grid()
        nstep = get_nstepstr(self.nStep)
        varinf= self.get_varinfo(datatype)
        data_dict={}
        n1 = grid_dict.get('n1')
        n2 = grid_dict.get('n2')
        n3 = grid_dict.get('n3')
        print "<DOMAIN> %d x %d x %d " % (n1,n2,n3)
        
        if varinf.get('fltype') == 'single_file':
	    if datatype == 'float':
		    fname_data = self.wdir+"data."+nstep+".flt"
	    else:
		    fname_data = self.wdir+"data."+nstep+".dbl"
	
            f_data = open(fname_data,'rb')
            datout = f_data.read()
	    if datatype == 'float':
		    D=struct.unpack("<"+str(len(datout)/4)+"f",datout)
	    else:
		    D=struct.unpack("<"+str(len(datout)/8)+"d",datout)
	
            A = np.asarray(D)
            for i in range(varinf.get('nvar')):
		    #print "> Reading %s" % (varinf.get('allvars')[i])
                if varinf.get('allvars')[i] == varinf.get('allvars')[-1] :
			if n3 == 1:
				data_dict[(varinf.get('allvars')[i])]=A[-n2*n1:].reshape(n2,n1).transpose()
			else:
				data_dict[(varinf.get('allvars')[i])]=A[-n3*n2*n1:].reshape(n3,n2,n1).transpose()
                else :
			if n3 == 1:
				data_dict[(varinf.get('allvars')[i])]=A[i*n2*n1:(i+1)*n2*n1].reshape(n2,n1).transpose()
			else:
				data_dict[(varinf.get('allvars')[i])]=A[i*n3*n2*n1:(i+1)*n3*n2*n1].reshape(n3,n2,n1).transpose()
            
        else:
            fname_list = []
            f_list = []
            datout =[]
            Dind=[]


            for j in range(varinf.get('nvar')):
		if datatype == 'float':
			fname_list.append(self.wdir+ varinf.get('allvars')[j]+"."+nstep+".flt")
			f_list.append(open(fname_list[j],'rb'))
			datout.append(f_list[j].read())
			Dind.append(struct.unpack("<"+str(len(datout[j])/4)+"f",datout[j]))
		else:
			fname_list.append(self.wdir+ varinf.get('allvars')[j]+"."+nstep+".dbl")
			f_list.append(open(fname_list[j],'rb'))
			datout.append(f_list[j].read())
			Dind.append(struct.unpack("<"+str(len(datout[j])/8)+"d",datout[j]))
	

	    
	    A = np.asarray(Dind)
	    
	    for j in range(varinf.get('nvar')):
		    if n3 == 1:
			    data_dict[(varinf.get('allvars')[j])]=A[j].reshape(n2,n1).transpose()
		    else:
			    data_dict[(varinf.get('allvars')[j])]=A[j].reshape(n3,n2,n1).transpose()
		    

	
	return data_dict
                
    def __init__(self,ns,w_dir=None,datatype=None):
	    self.nStep = ns
	    if w_dir is None :
		    self.wdir=curdir()
	    else:
		    self.wdir = w_dir

	
	    Grid_dictionary=self.grid()
	    Data_dictionary=self.data(datatype)
	    Time_dictionary = self.time_info(datatype)
	    for keys in Grid_dictionary:
		    object.__setattr__(self,keys,Grid_dictionary.get(keys))
	    for keys in Data_dictionary:
		    object.__setattr__(self,keys,Data_dictionary.get(keys))
	    for keys in Time_dictionary:
		    object.__setattr__(self,keys,Time_dictionary.get(keys))






class Tools(object):
	"""
	
	This Class has all the functions doing basic mathematical
	operations to the vector or scalar fields.
	It is called after pyPLUTO.pload object is defined.
	
	"""

	def deriv(self,Y,X=None):
		"""
		Calculates the derivative of Y with respect to X.

		**Inputs:**

		  Y : 1-D array to be differentiated.\n
		  X : 1-D array with len(X) = len(Y).\n

		  If X is not specified then by default X is chosen to be an equally spaced array having same number of elements
		  as Y.

                **Outputs:**

                  This returns an 1-D array having the same no. of elements as Y (or X) and contains the values of dY/dX.
		
		"""
		n = len(Y)
		n2 = n-2
		if X==None : X = np.arange(n)
		Xarr = np.asarray(X,dtype='float')
		Yarr = np.asarray(Y,dtype='float')
		x12 = Xarr - np.roll(Xarr,-1)   #x1 - x2
		x01 = np.roll(Xarr,1) - Xarr    #x0 - x1
		x02 = np.roll(Xarr,1) - np.roll(Xarr,-1) #x0 - x2
		DfDx = np.roll(Yarr,1) * (x12 / (x01*x02)) + Yarr * (1./x12 - 1./x01) - np.roll(Yarr,-1) * (x01 / (x02 * x12))
		# Formulae for the first and last points:

		DfDx[0] = Yarr[0] * (x01[1]+x02[1])/(x01[1]*x02[1]) - Yarr[1] * x02[1]/(x01[1]*x12[1]) + Yarr[2] * x01[1]/(x02[1]*x12[1])
		DfDx[n-1] = -Yarr[n-3] * x12[n2]/(x01[n2]*x02[n2]) + Yarr[n-2]*x02[n2]/(x01[n2]*x12[n2]) - Yarr[n-1]*(x02[n2]+x12[n2])/(x02[n2]*x12[n2])

		return DfDx
	
	def Grad(self,phi,x1,x2,dx1,dx2,polar=False):
		""" This method calculates the gradient of the 2D scalar phi.

                **Inputs:**

                  phi -- 2D scalar whose gradient is to be determined.\n
		  x1 -- The 'x' array\n
		  x2 -- The 'y' array\n
                  dx1 -- The grid spacing in 'x' direction.\n
                  dx2 -- The grid spacing in 'y' direction.\n
                  polar -- The keyword should be set to True inorder to estimate the Gradient in polar co-ordinates. By default it is set to False.
                
		**Outputs:**

                  This routine outputs a 3D array with shape = (len(x1),len(x2),2), such that [:,:,0] element corresponds to the gradient values of phi wrt to x1 and [:,:,1] are the gradient values of phi wrt to x2.
 
		"""
		(n1, n2) = phi.shape 
		grad_phi = np.zeros(shape=(n1,n2,2))
		h2 = np.ones(shape=(n1,n2))
		if polar == True:
			for j in range(n2):
				h2[:,j] = x1
		
		for i in range(n1):
			scrh1 = phi[i,:]
			grad_phi[i,:,1] = self.deriv(scrh1,x2)/h2[i,:]
		for j in range(n2):
			scrh2 = phi[:,j]
			grad_phi[:,j,0] = self.deriv(scrh2,x1)

		return grad_phi

	def Div(self,u1,u2,x1,x2,dx1,dx2,geometry=None):
		""" This method calculates the divergence of the 2D vector fields u1 and u2.

                **Inputs:**

                  u1 -- 2D vector along x1 whose divergence is to be determined.\n
                  u2 -- 2D vector along x2 whose divergence is to be determined.\n
		  x1 -- The 'x' array\n
		  x2 -- The 'y' array\n
                  dx1 -- The grid spacing in 'x' direction.\n
                  dx2 -- The grid spacing in 'y' direction.\n
                  geometry -- The keyword *geometry* is by default set to 'cartesian'. It can be set to either one of the following : *cartesian*, *cylindrical*, *spherical* or *polar*. To calculate the divergence of the vector fields, respective geometric corrections are taken into account based on the value of this keyword.

                **Outputs:**

		  A 2D array with same shape as u1(or u2) having the values of divergence.

		"""
		(n1, n2) = u1.shape
		Divergence = np.zeros(shape=(n1,n2))
		du1 = np.zeros(shape=(n1,n2))
		du2 = np.zeros(shape=(n1,n2))

		A1 = np.zeros(shape=n1)
		A2 = np.zeros(shape=n2)

		dV1 = np.zeros(shape=(n1,n2))
		dV2 = np.zeros(shape=(n1,n2))

		if geometry == None : geometry = 'cartesian'
		
		#------------------------------------------------
		#  define area and volume elements for the
		#  different coordinate systems
		#------------------------------------------------

		if geometry == 'cartesian' :
			A1[:] = 1.0
			A2[:] = 1.0
			dV1   = np.outer(dx1,A2)
			dV2   = np.outer(A1,dx2)

		if geometry == 'cylindrical' :
			A1 = x1
			A2[:] = 1.0
			dV1 = np.meshgrid(x1*dx1,A2)[0].T*np.meshgrid(x1*dx1,A2)[1].T
			for i in range(n1) : dV2[i,:] = dx2[:]
		
		if geometry == 'polar' :
			A1    = x1
			A2[:] = 1.0
			dV1   = np.meshgrid(x1,A2)[0].T*np.meshgrid(x1,A2)[1].T
			dV2   = np.meshgrid(x1,dx2)[0].T*np.meshgrid(x1,dx2)[1].T

		if geometry == 'spherical' :
			A1 = x1*x1
			A2 = np.sin(x2)
			for j in range(n2): dV1[:,j] = A1*dx1
			dV2   = np.meshgrid(x1,np.sin(x2)*dx2)[0].T*np.meshgrid(x1,np.sin(x2)*dx2)[1].T

		# ------------------------------------------------
		#              Make divergence
		# ------------------------------------------------
		
		
		for i in range(1,n1-1):
			du1[i,:] = 0.5*(A1[i+1]*u1[i+1,:] - A1[i-1]*u1[i-1,:])/dV1[i,:]
		for j in range(1,n2-1):
			du2[:,j] = 0.5*(A2[j+1]*u2[:,j+1] - A2[j-1]*u2[:,j-1])/dV2[:,j]

		Divergence = du1 + du2
		return Divergence



	def RTh2Cyl(self,R,Th,X1,X2):
		""" This method does the transformation from spherical coordinates to cylindrical ones.
		
                **Inputs:**

                  R - 2D array of spherical radius coordinates.\n
	          Th - 2D array of spherical theta-angle coordinates.\n
		  X1 - 2D array of radial component of given vector\n
		  X2 - 2D array of thetoidal component of given vector\n
                
		**Outputs:**

                  This routine outputs two 2D arrays after transformation.
				  
                **Usage:**

                  ``import pyPLUTO as pp``\n
		  ``import numpy as np``\n
                  ``D = pp.pload(0)``\n
		  ``ppt=pp.Tools()``\n
                  ``TH,R=np.meshgrid(D.x2,D.x1)``\n
		  ``Br,Bz=ppt.RTh2Cyl(R,TH,D.bx1,D.bx2)``

                D.bx1 and D.bx2 should be vectors in spherical coordinates. After transformation (Br,Bz) corresponds to vector in cilindrical coordinates.

 
		"""
		Y1=X1*np.sin(Th)+X2*np.cos(Th)
		Y2=X1*np.cos(Th)-X2*np.sin(Th)
		return Y1,Y2




	def myInterpol(self,RR,N):
		""" This method interpolates (linear interpolation) vector 1D vector RR to 1D N-length vector. Useful for stretched grid calculations. 
		
                **Inputs:**

                  RR - 1D array to interpolate.\n
		  N  - Number of grids to interpolate to.\n
				  
		**Outputs:**

                  This routine outputs interpolated 1D array to the new grid (len=N).
				  
                **Usage:**

                  ``import pyPLUTO as pp``\n
		  ``import numpy as np``\n
                  ``D = pp.pload(0)``\n
		  ``ppt=pp.Tools()``\n
                  ``x=linspace(0,1,10) #len(x)=10``\n
		  ``y=x*x``\n
		  ``Ri,Ni=ppt.myInterpol(y,100) #len(Ri)=100``

                  Ri - interpolated numbers;
		  Ni - grid for Ri
 
		"""	
		
		NN=np.linspace(0,len(RR)-1,len(RR))
		spline_fit=UnivariateSpline(RR,NN,k=3,s=0)
		
		RRi=np.linspace(RR[0],RR[-1],N)
		NNi=spline_fit(RRi)
		NNi[0]=NN[0]+0.00001
		NNi[-1]=NN[-1]-0.00001
		return RRi,NNi
		
	def getUniformGrid(self,r,th,rho,Nr,Nth):
		""" This method transforms data with non-uniform grid (stretched) to uniform. Useful for stretched grid calculations. 
		
                **Inputs:**

		  r  - 1D vector of X1 coordinate (could be any, e.g D.x1).\n
		  th - 1D vector of X2 coordinate (could be any, e.g D.x2).\n
		  rho- 2D array of data.\n
		  Nr - new size of X1 vector.\n
		  Nth- new size of X2 vector.\n
				  
		**Outputs:**

                  This routine outputs 2D uniform array Nr x Nth dimension
				  
                **Usage:**

                  ``import pyPLUTO as pp``\n
		  ``import numpy as np``\n
                  ``D = pp.pload(0)``\n
		  ``ppt=pp.Tools()``\n
		  ``X1new, X2new, res = ppt.getUniformGrid(D.x1,D.x2,D.rho,20,30)``

                  X1new - X1 interpolated grid len(X1new)=20
		  X2new - X2 interpolated grid len(X2new)=30
		  res   - 2D array of interpolated variable
 
		"""	

		Ri,NRi=self.myInterpol(r,Nr)
		Ra=np.int32(NRi);Wr=NRi-Ra

		YY=np.ones([Nr,len(th)])
		for i in range(len(th)):
		      YY[:,i]=(1-Wr)*rho[Ra,i] + Wr*rho[Ra+1,i]

		THi,NTHi=self.myInterpol(th,Nth)
		THa=np.int32(NTHi);Wth=NTHi-THa

		ZZ=np.ones([Nr,Nth])
		for i in range(Nr):
		      ZZ[i,:]=(1-Wth)*YY[i,THa] + Wth*YY[i,THa+1]

		return Ri,THi,ZZ
	
        def sph2cyl(self,D,Dx,rphi=None,theta0=None):
		""" This method transforms spherical data into cylindrical applying interpolation. Works for stretched grid as well, transforms poloidal (R-Theta) data by default. Fix theta and set rphi=True to get (R-Phi) transformation.
				
                **Inputs:**

                  D  - structure  from 'pload' method.\n
		  Dx - variable to be transformed (D.rho for example).\n
				  
		**Outputs:**

                  This routine outputs transformed (sph->cyl) variable and grid.
				  
                **Usage:**

                  ``import pyPLUTO as pp``\n
		  ``import numpy as np``\n
                  ``D = pp.pload(0)``\n
		  ``ppt=pp.Tools()``\n
		  ``R,Z,res = ppt.sph2cyl(D,D.rho.transpose())``

                  R - 2D array with cylindrical radius values
		  Z - 2D array with cylindrical Z values
		  res - 2D array of transformed variable
 
		"""	
		
		if rphi is None or rphi == False:
		    rx=D.x1
		    th=D.x2		    
		else:
                    rx=D.x1*np.sin(theta0)
		    th=D.x3
		    
		rx,th,Dx=self.getUniformGrid(rx,th,Dx.T,200,200)
		Dx=Dx.T
		
		if rphi is None or rphi == False:
                    
                    r0=np.min(np.sin(th)*rx[0])
                    rN=rx[-1]
                    dr=rN-r0
                    z0=np.min(np.cos(th)*rN)
                    zN=np.max(np.cos(th)*rN)
                    dz=zN-z0
                    dth=th[-1]-th[0]
                    rl=np.int32(len(rx)*dr/(rx[-1]-rx[0]))  
                    zl=np.int32(rl* dz/dr)
                    thl=len(th)
                    r=np.linspace(r0, rN, rl)
                    z=np.linspace(z0, zN, zl)
		else:
                    r0=np.min([np.sin(th)*rx[0] , np.sin(th)*rx[-1]])
                    rN=np.max([np.sin(th)*rx[0] , np.sin(th)*rx[-1]])
                    dr=rN-r0
                    z0=np.min(np.cos(th)*rN)
                    zN=np.max(np.cos(th)*rN)
                    dz=zN-z0
                    dth=th[-1]-th[0]
                    rl=np.int32(len(rx)*dr/(rx[-1]-rx[0]))  
                    zl=np.int32(rl* dz/dr)
                    thl=len(th)
                    r=np.linspace(r0, rN, rl)
                    z=np.linspace(z0, zN, zl)
                
                R,Z = np.meshgrid(r, z)
		Rs = np.sqrt(R*R + Z*Z)
		
		
		Th = np.arccos(Z/Rs)
		kv_34=find(R<0)
		Th.flat[kv_34]=2*np.pi - Th.flat[kv_34]
		
		
		ddr=rx[1]-rx[0]
		ddth=th[1]-th[0]
		
		Rs_copy=Rs.copy()
		Th_copy=Th.copy()
				
		nR1=find(Rs<rx[0])  
		Rs.flat[nR1]=rx[0] 
		nR2=find(Rs>rN)
		Rs.flat[nR2]=rN
		
		nTh1=find(Th>th[-1])
		Th.flat[nTh1]=th[-1]
		nTh2=find(Th<th[0])
		Th.flat[nTh2]=th[0]
		
		
		ra = ((len(rx)-1.001)/(np.max(Rs.flat)-np.min(Rs.flat)) *(Rs-np.min(Rs.flat)))  
		tha = ((thl-1.001)/dth *(Th-th[0]))  

		rn = np.int32(ra)
		thn = np.int32(tha)
		dra=ra-rn
		dtha=tha-thn
		w1=1-dra
		w2=dra
		w3=1-dtha
		w4=dtha
		lrx=len(rx)
		NN1=np.int32(rn+thn*lrx)
		NN2=np.int32((rn+1)+thn*lrx)
		NN3=np.int32(rn+(thn+1)*lrx)
		NN4=np.int32((rn+1)+(thn+1)*lrx)
		n=np.transpose(np.arange(0,np.product(np.shape(R))))
		DD=Dx.copy()
		F=R.copy()
		F.flat[n]=w1.flat[n]*(w3.flat[n]*Dx.flat[NN1.flat[n]] + w4.flat[n]*Dx.flat[NN3.flat[n]] )+\
		    w2.flat[n]*(w3.flat[n]*Dx.flat[NN2.flat[n]] + w4.flat[n]*Dx.flat[NN4.flat[n]] )
		    
		nR1=find(Rs_copy<rx[0]-ddr/1.5)
		nR2=find(Rs_copy>rN+ddr/1.5)
		nTh1=find(Th_copy>th[-1]+ddth/1.5)
		nTh2=find(Th_copy<th[0]-ddth/1.5)

		nmask=np.concatenate((nR1,nR2,nTh1,nTh2))
		F.flat[nmask]=np.nan;
		return R,Z,F
        

	def congrid(self, a, newdims, method='linear', centre=False, minusone=False):
	    """
	    Arbitrary resampling of source array to new dimension sizes.
	    Currently only supports maintaining the same number of dimensions.
	    To use 1-D arrays, first promote them to shape (x,1).

	    Uses the same parameters and creates the same co-ordinate lookup points
	    as IDL''s congrid routine, which apparently originally came from a VAX/VMS
	    routine of the same name.

	    **Inputs:**

	      a -- The 2D array that needs resampling into new dimensions.\n
	      newdims -- A tuple which represents the shape of the resampled data\n
	      method -- This keyword decides the method used for interpolation.\n
	        neighbour - closest value from original data\n
	        nearest and linear - uses n x 1-D interpolations using scipy.interpolate.interp1d
	        (see Numerical Recipes for validity of use of n 1-D interpolations)\n
	        spline - uses ndimage.map_coordinates\n
	      centre -- This keyword decides the positions of interpolation points.\n
	        True - interpolation points are at the centres of the bins\n
	        False - points are at the front edge of the bin\n
	      minusone -- This prevents extrapolation one element beyond bounds of input array\n
	        For example- inarray.shape = (i,j) & new dimensions = (x,y)\n
	        False - inarray is resampled by factors of (i/x) * (j/y)\n
	        True - inarray is resampled by(i-1)/(x-1) * (j-1)/(y-1)\n

            **Outputs:**

              A 2D array with resampled data having a shape corresponding to newdims. 
	    
	    """
	    if not a.dtype in [np.float64, np.float32]:
		a = np.cast[float](a)

	    m1 = np.cast[int](minusone)
	    ofs = np.cast[int](centre) * 0.5
	    old = np.array( a.shape )
	    ndims = len( a.shape )
	    if len( newdims ) != ndims:
		print "[congrid] dimensions error. " \
		      "This routine currently only support " \
		      "rebinning to the same number of dimensions."
		return None
	    newdims = np.asarray( newdims, dtype=float )
	    dimlist = []

	    if method == 'neighbour':
		for i in range( ndims ):
		    base = np.indices(newdims)[i]
		    dimlist.append( (old[i] - m1) / (newdims[i] - m1) \
				    * (base + ofs) - ofs )
		cd = np.array( dimlist ).round().astype(int)
		newa = a[list( cd )]
		return newa

	    elif method in ['nearest','linear']:
		# calculate new dims
		for i in range( ndims ):
		    base = np.arange( newdims[i] )
		    dimlist.append( (old[i] - m1) / (newdims[i] - m1) \
				    * (base + ofs) - ofs )
		# specify old dims
		olddims = [np.arange(i, dtype = np.float) for i in list( a.shape )]

		# first interpolation - for ndims = any
		mint = scipy.interpolate.interp1d( olddims[-1], a, kind=method )
		newa = mint( dimlist[-1] )

		trorder = [ndims - 1] + range( ndims - 1 )
		for i in range( ndims - 2, -1, -1 ):
		    newa = newa.transpose( trorder )

		    mint = scipy.interpolate.interp1d( olddims[i], newa, kind=method )
		    newa = mint( dimlist[i] )

		if ndims > 1:
		    # need one more transpose to return to original dimensions
		    newa = newa.transpose( trorder )

		return newa
	    elif method in ['spline']:
		oslices = [ slice(0,j) for j in old ]
		oldcoords = np.ogrid[oslices]
		nslices = [ slice(0,j) for j in list(newdims) ]
		newcoords = np.mgrid[nslices]

		newcoords_dims = range(n.rank(newcoords))
		#make first index last
		newcoords_dims.append(newcoords_dims.pop(0))
		newcoords_tr = newcoords.transpose(newcoords_dims)
		# makes a view that affects newcoords

		newcoords_tr += ofs

		deltas = (np.asarray(old) - m1) / (newdims - m1)
		newcoords_tr *= deltas

		newcoords_tr -= ofs

		newa = scipy.ndimage.map_coordinates(a, newcoords)
		return newa
	    else:
		print "Congrid error: Unrecognized interpolation type.\n", \
		      "Currently only \'neighbour\', \'nearest\',\'linear\',", \
		      "and \'spline\' are supported."
		return None


        








class Image(object):
	''' This Class has all the routines for the imaging the data
	and plotting various contours and fieldlines on these images.
	CALLED AFTER pyPLUTO.pload object is defined
	'''
	def pldisplay(self,var,**kwargs):
		""" This method allows the user to display a 2D data using the matplotlib's pcolormesh.

		**Inputs:**
		
		  var -- 2D array that needs to be displayed.

		*Required Keywords:*

		  x1 -- The 'x' array\n
		  x2 -- The 'y' array\n

		*Optional Keywords:*
		
		  vmin -- The minimum value of the 2D array (Default : min(var))\n
                  vmax -- The maximum value of the 2D array (Default : max(var))\n
                  title -- Sets the title of the image.\n
                  label1 -- Sets the X Label (Default: 'XLabel')\n
                  label2 -- Sets the Y Label (Default: 'YLabel')\n
                  cbar -- Its a tuple to set the colorbar on or off. 
                     cbar = (True,'vertical') -- Displays a vertical colorbar\n
                     cbar = (True,'horizontal') -- Displays a horizontal colorbar\n
                     cbar = (False,'') -- Displays no colorbar.\n
		
                **Usage:**
		
                  ``import pyPLUTO as pp``\n
                  ``wdir = '/path/to/the data files/'``\n
                  ``D = pp.pload(1,w_dir=wdir)``\n
                  ``I = pp.Image()``\n
                  ``f1 = figure()``\n
                  ``ax1 = f1.add_subplot(111)``\n
                  ``I.pldisplay(D.v2,x1=D.x1,x2=D.x2,cbar=(True,'vertical'),title='Velocity',label1='Radius',label2='Height')``\n

        	"""
		
		x1 = kwargs.get('x1')
		x2 = kwargs.get('x2')
		var = var.T

		f1 = figure(kwargs.get('fignum',1), figsize=kwargs.get('figsize',[10,10]),dpi=80, facecolor='w', edgecolor='k')
		ax1 = f1.add_subplot(111)
		ax1.set_aspect('equal')
		ax1.axis([np.min(x1),np.max(x1),np.min(x2),np.max(x2)])
		pcolormesh(x1,x2,var,vmin=kwargs.get('vmin',np.min(var)),vmax=kwargs.get('vmax',np.max(var)))
		
		title(kwargs.get('title',"Title"),size=kwargs.get('size'))
		xlabel(kwargs.get('label1',"Xlabel"),size=kwargs.get('size'))
		ylabel(kwargs.get('label2',"Ylabel"),size=kwargs.get('size'))
		if kwargs.get('cbar',(False,''))[0] == True:
			colorbar(orientation=kwargs.get('cbar')[1])

	

	def multi_disp(self,*args,**kwargs):
		mvar = []
		var_cart_list=[]
		for arg in args:
			mvar.append(arg.T)
	
		
		xmin = np.min(kwargs.get('x1'))
		xmax = np.max(kwargs.get('x1'))		
		ymin = np.min(kwargs.get('x2'))
		ymax = np.max(kwargs.get('x2'))
		mfig = figure(kwargs.get('fignum',1),figsize=kwargs.get('figsize',[10,10]))
		Ncols = kwargs.get('Ncols')
		Nrows = len(args)/Ncols
		mprod = Nrows*Ncols
		dictcbar=kwargs.get('cbar',(False,'','each'))
		for j in range(mprod):
			mfig.add_subplot(Nrows,Ncols,j+1)
			pcolormesh(kwargs.get('x1'),kwargs.get('x2'), mvar[j])
			axis([xmin,xmax,ymin,ymax])
			gca().set_aspect('equal')
			
			xlabel(kwargs.get('label1',mprod*['Xlabel'])[j])
			ylabel(kwargs.get('label2',mprod*['Ylabel'])[j])
			title(kwargs.get('title',mprod*['Title'])[j])
			if (dictcbar[0] == True) and (dictcbar[2] =='each'):
				colorbar(orientation=kwargs.get('cbar')[1])
			if dictcbar[0] == True and dictcbar[2]=='last':
					if (j == np.max(range(mprod))):colorbar(orientation=kwargs.get('cbar')[1])
				
	def field_interp(self,var1,var2,x,y,dx,dy,xp,yp):
		""" This method interpolates value of vector fields (var1 and var2) on field points (xp and yp).
		The field points are obtained from the method field_line.

		**Inputs:**
		
		  var1 -- 2D Vector field in X direction\n
		  var2 -- 2D Vector field in Y direction\n
		  x -- 1D X array\n
		  y -- 1D Y array\n
		  dx -- 1D grid spacing array in X direction\n
		  dy -- 1D grid spacing array in Y direction\n
		  xp -- field point in X direction\n
		  yp -- field point in Y direction\n

		**Outputs:**

                  A list with 2 elements where the first element corresponds to the interpolate field point in 'x' direction and the second element is the field point in 'y' direction.  

		"""
		q=[]
		U = var1
		V = var2
		i0 = np.abs(xp-x).argmin()
		j0 = np.abs(yp-y).argmin()
		scrhUx = np.interp(xp,x,U[:,j0])
		scrhUy = np.interp(yp,y,U[i0,:])
		q.append(scrhUx + scrhUy - U[i0,j0])
		scrhVx = np.interp(xp,x,V[:,j0])
		scrhVy = np.interp(yp,y,V[i0,:])
		q.append(scrhVx + scrhVy - V[i0,j0])
		return q

	def field_line(self,var1,var2,x,y,dx,dy,x0,y0):
		""" This method is used to obtain field lines (same as fieldline.pro in PLUTO IDL tools).

		**Inputs:**

		  var1 -- 2D Vector field in X direction\n
		  var2 -- 2D Vector field in Y direction\n
		  x -- 1D X array\n
		  y -- 1D Y array\n
		  dx -- 1D grid spacing array in X direction\n
		  dy -- 1D grid spacing array in Y direction\n
		  x0 -- foot point of the field line in X direction\n
		  y0 -- foot point of the field line in Y direction\n

                **Outputs:**

                  This routine returns a dictionary with keys - \n
                  qx -- list of the field points along the 'x' direction.
                  qy -- list of the field points along the 'y' direction.

                **Usage:**

		  See the myfieldlines routine for the same. 
                 

		"""
		xbeg = x[0] - 0.5*dx[0]
		xend = x[-1] + 0.5*dx[-1]

		ybeg = y[0]  - 0.5*dy[0]
		yend = y[-1] + 0.5*dy[-1]

		inside_domain = x0 > xbeg and x0 < xend and y0 > ybeg and y0 < yend

		MAX_STEPS = 5000
		xln_fwd = [x0]
		yln_fwd = [y0]
		xln_bck = [x0]
		yln_bck = [y0]
		rhs = []
		k = 0

		while inside_domain == True:
		    R1 = self.field_interp(var1,var2,x,y,dx,dy,xln_fwd[k],yln_fwd[k])
		    dl = 0.5*np.max(np.concatenate((dx,dy)))/(np.sqrt(R1[0]*R1[0] + R1[1]*R1[1] + 1.e-14))
		    xscrh = xln_fwd[k] + 0.5*dl*R1[0]
		    yscrh = yln_fwd[k] + 0.5*dl*R1[1]

		    R2 = self.field_interp(var1,var2,x,y,dx,dy,xln_fwd[k],yln_fwd[k])

		    xln_one = xln_fwd[k] + dl*R2[0]
		    yln_one = yln_fwd[k] + dl*R2[1]

		    xln_fwd.append(xln_one)
		    yln_fwd.append(yln_one)
		    inside_domain = xln_one > xbeg and xln_one < xend and yln_one > ybeg and yln_one < yend
		    inside_domain = inside_domain and (k < MAX_STEPS-3)
		    k = k + 1


		k_fwd = k
		qx = np.asarray(xln_fwd[0:k_fwd])
		qy = np.asarray(yln_fwd[0:k_fwd])
		flines={'qx':qx,'qy':qy}
		return flines


	def myfieldlines(self,Data,x0arr,y0arr,stream=False,**kwargs):
		""" This method overplots the magnetic field lines at the footpoints given by (x0arr[i],y0arr[i]).

		**Inputs:**
		
		  Data -- pyPLUTO.pload object\n
		  x0arr -- array of x co-ordinates of the footpoints\n
		  y0arr -- array of y co-ordinates of the footpoints\n
		  stream -- keyword for two different ways of calculating the field lines.\n
		    True -- plots contours of rAphi (needs to store vector potential)\n
		    False -- plots the fieldlines obtained from the field_line routine. (Default option)\n

		*Optional Keywords:*
		
		  colors -- A list of matplotlib colors to represent the lines. The length of this list should be same as that of x0arr.\n
		  lw -- Integer value that determines the linewidth of each line.\n
		  ls -- Determines the linestyle of each line.\n

		**Usage:**

		Assume that the magnetic field is a given as **B** = B0$\hat{y}$. Then to show this field lines we have to define the x and y arrays of field foot points.

		``x0arr = linspace(0.0,10.0,20)``\n
		``y0arr = linspace(0.0,0.0,20)``\n
		``import pyPLUTO as pp``\n
		``D = pp.pload(45)``\n
		``I = pp.Image()``\n
		``I.myfieldlines(D,x0arr,y0arr,colors='k',ls='--',lw=1.0)``

	       """
	       
		if len(x0arr) != len(y0arr) : print "Input Arrays should have same size"
		QxList=[]
		QyList=[]
		StreamFunction = []
		levels =[]
		if stream == True:
			X, Y = np.meshgrid(Data.x1,Data.x2.T)
			StreamFunction = X*(Data.Ax3.T)
			for i in range(len(x0arr)):
				nx = np.abs(X[0,:]-x0arr[i]).argmin()
				ny = np.abs(X[:,0]-y0arr[i]).argmin()
				levels.append(X[ny,nx]*Data.Ax3.T[ny,nx])
			
			contour(X,Y,StreamFunction,levels,colors=kwargs.get('colors'),linewidths=kwargs.get('lw',1),linestyles=kwargs.get('ls','solid'))
		else:
			for i in range(len(x0arr)):
				QxList.append(self.field_line(Data.bx1,Data.bx2,Data.x1,Data.x2,Data.dx1,Data.dx1,x0arr[i],y0arr[i]).get('qx'))
				QyList.append(self.field_line(Data.bx1,Data.bx2,Data.x1,Data.x2,Data.dx1,Data.dx1,x0arr[i],y0arr[i]).get('qy'))
				plot(QxList[i],QyList[i],color=kwargs.get('colors'))
			axis([min(Data.x1),max(Data.x1),min(Data.x2),max(Data.x2)])

        def getSphData(self,Data,w_dir=None,datatype=None,**kwargs):
	    """This method transforms the vector and scalar  fields from Spherical co-ordinates to Cylindrical.

	    **Inputs**:
	    
	      Data -- pyPLUTO.pload object\n
	      w_dir -- /path/to/the/working/directory/\n
	      datatype -- If the data is of 'float' type then datatype = 'float' else by default the datatype is set to 'double'.

	    *Optional Keywords*:
	    
	      rphi -- [Default] is set to False implies that the r-theta plane is transformed. If set True then the r-phi plane is transformed.\n
              x2cut -- Applicable for 3D data and it determines the co-ordinate of the x2 plane while r-phi is set to True.\n
	      x3cut -- Applicable for 3D data and it determines the co-ordinate of the x3 plane while r-phi is set to False.\n


            """
 
            Tool = Tools()
            key_value_pairs = []
            if w_dir is None: w_dir = curdir()		    
            allvars = Data.get_varinfo(datatype).get('allvars')
            
            
            if kwargs.get('rphi',False)==True:
		R,TH = np.meshgrid(Data.x1,Data.x3)
                if Data.n3 != 1:
                    for variable in allvars:
                        key_value_pairs.append([variable,getattr(Data,variable)[:,kwargs.get('x2cut',0),:].T])
                    SphData = dict(key_value_pairs)
		    if ('bx1' in allvars) or ('bx2' in allvars):
                        (SphData['b1c'],SphData['b3c']) = Tool.RTh2Cyl(R,TH,SphData.get('bx1'),SphData.get('bx3'))
			allvars.append('b1c')
			allvars.append('b3c')
                    if ('vx1' in allvars) or ('vx2' in allvars):
                        (SphData['v1c'],SphData['v3c']) = Tool.RTh2Cyl(R,TH,SphData.get('vx1'),SphData.get('vx3'))
			allvars.append('v1c')
			allvars.append('v3c')
                else:
                    print "No x3 plane for 2D data"
            else:
                R,TH = np.meshgrid(Data.x1,Data.x2)
                if Data.n3 != 1:
                    for variable in allvars:
                        key_value_pairs.append([variable,getattr(Data,variable)[:,:,kwargs.get('x3cut',0)].T])
                    SphData = dict(key_value_pairs)
                    if ('bx1' in allvars) or ('bx2' in allvars):
                        (SphData['b1c'],SphData['b2c']) = Tool.RTh2Cyl(R,TH,SphData.get('bx1'),SphData.get('bx2'))
			allvars.append('b1c')
			allvars.append('b2c')
                    if ('vx1' in allvars) or ('vx2' in allvars):
                        (SphData['v1c'],SphData['v2c']) = Tool.RTh2Cyl(R,TH,SphData.get('vx1'),SphData.get('vx2'))
			allvars.append('v1c')
			allvars.append('v2c')
                else:
                    for variable in allvars:
                        key_value_pairs.append([variable,getattr(Data,variable)[:,:].T])
                    SphData = dict(key_value_pairs)
                    if ('bx1' in allvars) or ('bx2' in allvars):
                        (SphData['b1c'],SphData['b2c']) = Tool.RTh2Cyl(R,TH,SphData.get('bx1'),SphData.get('bx2'))
			allvars.append('b1c')
			allvars.append('b2c')
                    if ('vx1' in allvars) or ('vx2' in allvars):
                        (SphData['v1c'],SphData['v2c']) = Tool.RTh2Cyl(R,TH,SphData.get('vx1'),SphData.get('vx2'))
			allvars.append('v1c')
			allvars.append('v2c')
            
            for variable in allvars:
                if kwargs.get('rphi',False)==True:
		    R,Z,SphData[variable]= Tool.sph2cyl(Data,SphData.get(variable),rphi=True,theta0=Data.x2[kwargs.get('x2cut',0)])
                else:
                    if Data.n3 != 1:
                        R,Z,SphData[variable] = Tool.sph2cyl(Data,SphData.get(variable),rphi=False)
                    else:
                        R,Z,SphData[variable] = Tool.sph2cyl(Data,SphData.get(variable),rphi=False)

            return R,Z,SphData



        def pltSphData(self,Data,w_dir=None,datatype=None,**kwargs):
            """This method plots the transformed data obtained from getSphData using the matplotlib's imshow

            **Inputs:**

              Data -- pyPLUTO.pload object\n
              w_dir -- /path/to/the/working/directory/\n
	      datatype -- If the data is of 'float' type then datatype = 'float' else by default the datatype is set to 'double'.

            *Required Keywords*:

              plvar -- A string which represents the plot variable.\n

	    *Optional Keywords*:

	      logvar -- [Default = False] Set it True for plotting the log of a variable.
	      rphi -- [Default = False - for plotting in r-theta plane] Set it True for plotting the variable in r-phi plane. 

           """
              
            if w_dir is None: w_dir=curdir()
            R,Z,SphData = self.getSphData(Data,w_dir=w_dir,datatype=datatype,**kwargs)
    
            extent=(np.min(R.flat),max(R.flat),np.min(Z.flat),max(Z.flat))
            dRR=max(R.flat)-np.min(R.flat)
            dZZ=max(Z.flat)-np.min(Z.flat)


            isnotnan=-np.isnan(SphData[kwargs.get('plvar')])
            maxPl=max(SphData[kwargs.get('plvar')][isnotnan].flat)
            minPl=np.min(SphData[kwargs.get('plvar')][isnotnan].flat)
            normrange=False
            if minPl<0:
                normrange=True
            if maxPl>-minPl:
                minPl=-maxPl
            else:
                maxPl=-minPl	  
            if (normrange and kwargs.get('plvar')!='rho' and kwargs.get('plvar')!='prs'):
                SphData[kwargs.get('plvar')][-1][-1]=maxPl
                SphData[kwargs.get('plvar')][-1][-2]=minPl
	
	    if (kwargs.get('logvar') == True):
		    SphData[kwargs.get('plvar')] = np.log10(SphData[kwargs.get('plvar')])
	
            im1= imshow(SphData[kwargs.get('plvar')], aspect='equal', origin='lower', cmap=cm.jet,extent=extent, interpolation='nearest')



	

	
		
		



	


      




		    

	    
    
	    
	    
	    
	    
	    
	    



	
        
	  






    
         
    







    
   
