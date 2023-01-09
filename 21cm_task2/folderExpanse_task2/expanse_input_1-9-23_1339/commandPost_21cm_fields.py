import numpy as np
import os
import matplotlib.pyplot as plt
from matplotlib.pyplot import contourf 
import math
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
print("\nChosen z axis distance:",z_chosen,"Mpc/h\n") #~7
fig, ax = plt.subplots(1,3)#3 rows 1 column
#create the xy grid
x = np.linspace(0,N,N)
y = np.linspace(0,N,N)
xy_D = fields_D[:,:,z_chosen] #z input to the contourf
xy_F = fields_F[:,:,z_chosen]
xy_O = fields_O[:,:,z_chosen]
print("x:")
callout(x)
print("xy_D:")
callout(xy_D)

#find max values for color bar
print("xy_D min: ",xy_D.min()," xy_D max: ",xy_D.max())
print("xy_F min: ",xy_F.min()," xy_F max: ",xy_F.max())
print("xy_O min: ",xy_O.min()," xy_O max: ",xy_O.max())
max_array = np.array([xy_D.max(),xy_F.max(),xy_O.max()])
max_xy = max_array.max()
min_xy = 0
print("maximum:",max_xy)
print("minimum:",min_xy)

#create filled contour map
#REQUIRED: len(x) = num of columns of z and len(y) = num of rows of z, where z = xy
#also let the max and min values be the same for all models such that pnly one color bar is needed (they are all on same scale)
con_D = ax[0].contourf(x,y,xy_D,cmap="jet",vmin=min_xy,vmax=max_xy)
con_F = ax[1].contourf(x,y,xy_F,cmap="jet",vmin=min_xy,vmax=max_xy)
con_O = ax[2].contourf(x,y,xy_O,cmap="jet",vmin=min_xy,vmax=max_xy)

#set aspect ratio of x and y axis to be the same
ax[0].set_aspect('equal',anchor=(0.25,0.7)) #range from 0 (left/bottom) to 1 (right/top) for anchor(cx,cy)
ax[1].set_aspect('equal',anchor=(0.5,0.7))
ax[2].set_aspect('equal',anchor=(0.75,0.7))

#change padding/spacing between axes in the figure
#fig.tight_layout(pad=2.0)
fig.tight_layout(pad=1.0,w_pad=0.05,h_pad=1.0)

#create a color bar for all subplots
num = plt.cm.colors.Normalize(vmin=min_xy, vmax=max_xy)
sm = plt.cm.ScalarMappable(norm=num,cmap='jet')
ax4 = fig.add_axes([0.08,0.2,0.7,0.035]) #x,y,width,height
tickList = list(range(int(min_xy),int(math.floor(max_xy)),5)) #select major ticks for colorbar 
cbar = fig.colorbar(sm,orientation='horizontal',cax=ax4,label="mK",ticks=tickList)
'''
fig.colorbar(con_D,ax=ax[0],orientation='horizontal')
fig.colorbar(con_F,ax=ax[1],orientation='horizontal')
fig.colorbar(con_O,ax=ax[2],orientation='horizontal')
'''

#change figure size 
fig.set_size_inches(11,6)


#add titles
fig.suptitle("21cm Brightness Temperature (mK) at "+str(z_chosen)+" Mpc/h and z="+str(rs_chosen),y=0.93)
    #figure coordinates (bottom left corner is (0,0) and top right corner is (1,1)
    #https://stackoverflow.com/questions/55767312/how-to-position-suptitle
ax[0].set_ylabel("Mpc/h")
ax[0].set_xlabel("Mpc/h")
ax[0].set_title("Democratic")
ax[1].set_xlabel("Mpc/h")
ax[1].set_title("Fiducial")
ax[2].set_xlabel("Mpc/h")
ax[2].set_title("Oligarchic")
'''
#set log scale
[ax[i].set_xscale("log") for i in range(3)] #comprehension 
[ax[i].set_yscale("log") for i in range(3)]
'''
plt.show()
#save and output the plot 
outfile = "data_contourMap.pdf"
fig.savefig(outfile,format='pdf')
print("\noutput: ",outfile)

