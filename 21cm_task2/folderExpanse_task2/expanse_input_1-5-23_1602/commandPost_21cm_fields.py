import numpy as np
import os
import matplotlib.pyplot as plt
from matplotlib.pyplot import contourf 
#from matplotlib.pyplot import imshow 
#from matplotlib.pyplot import pcolormesh 

#--------------------------------------------------------------------------------------------------------------------------------------
#Functions
#--------------------------------------------------------------------------------------------------------------------------------------
def callout(bob):
    print("type:",type(bob)) #data type, should be numpy array  
    print("dim:",bob.ndim) #dimensions 
    print("shape:",bob.shape) #row by column 
    print("size:",bob.size,"\n") #number of elements 

def loadIn(model,rs,N):
    #each data file is a binary file with N=300^3 32-bit floats 
    #import and reshape 
    #np.fromfile(filename,dtype=np.float32).reshape((N,N,N))
    directory = "./"+model+"/t21_field.z="+rs
    rawData = np.fromfile(directory,dtype=np.float32).reshape((N,N,N))
    return rawData 

def getRS(model):
    os.chdir(model+"/") #change to directory of data 
    rs = (os.listdir()) #return a list of all files in that directory
    rs = [(x.replace("t21_field.z=","")) for x in rs] #use comprehension to replace all text save for redshift with blank
    os.chdir("../") #return back to original directory 
    return rs 



#--------------------------------------------------------------------------------------------------------------------------------------
#Program 
#--------------------------------------------------------------------------------------------------------------------------------------

#get each model's redshifts 
#each has 98 elements 
rs_D = np.sort(np.array(getRS("democratic")))
rs_F = np.sort(np.array(getRS("fiducial")))
rs_O = np.sort(np.array(getRS("oligarchic")))
'''
print("rs_D")
callout(rs_D)
print(rs_D)
print("\nrs_F")
callout(rs_F)
print(rs_F)
print("\nrs_O")
callout(rs_O)
print(rs_O)
'''

#select a redshift 
rs_chosen = rs_D[58] 
print("\nChosen redshift:",rs_chosen) #~7
print("Confirm chosen redshift is same for each model:",rs_D[58]==rs_F[58]," ",rs_D[58]==rs_O[58],"\n")

#acquire the data 
#the data is a grid of 300 by 300 by 300 for x,y,z axes 
#with 21 cm brightness temp at each location in units of mK
#the grid is evenly spaced from 0 to 300 Mpc/h for each direction
N = 300 #side-length of box or length of each x,y,or z axis in the data
fields_D = np.array(loadIn("democratic",rs_chosen,N))
fields_F = np.array(loadIn("fiducial",rs_chosen,N))
fields_O = np.array(loadIn("oligarchic",rs_chosen,N))
print("fields_D")
callout(fields_D)

#form the plot (contour map) using pyplot contourf
#chose a z to be fixed for the map
#therefore the map is a cross section, fixing the z value and looking at the temperatures along xy plane
z_chosen = 10 #index value or distance along z axis in Mpc/h
fig, ax = plt.subplots(1,3)
#create the xy grid
x = np.linspace(0,N,N)
y = np.linspace(0,N,N)
xy = fields_D[:,:,z_chosen] #z input to the contourf
print("x:")
callout(x)
print("xy:")
callout(xy)
#create filled contour map
#note len(x) = num of columns of z and len(y) = num of rows of z
con_D = ax[1].contourf(x,y,xy,cmap="Greens")
fig.colorbar(con_D)
plt.show()

#save and output the plot 
outfile = "data_contourMap.pdf"
fig.savefig(outfile,format='pdf')
print("output: ",outfile)

