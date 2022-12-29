from json import load
import matplotlib.pyplot as plt #for plotting
from scipy import interpolate 
import numpy as np 
import os 

#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Functions 
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

def callout(bob):
    print("type:",type(bob)) #data type, should be numpy array  
    print("dim:",bob.ndim) #dimensions 
    print("shape:",bob.shape) #row by column 
    print("size:",bob.size) #number of elements 
    #print("\n")
    print(bob,"\n")

def loadIn(model,rs):
    directory = ".//t21_power_data/t21_power_" + model + "/auto_t21_z=" + rs 
    rawData = np.loadtxt(directory,skiprows=1)
    return rawData 

def plot21cm_vs_k(model,rs,axis,colorName,i,j):
    rawData = loadIn(model,rs)
    x = rawData[:,1] #rowmin:rowmax, column 2 "k"
    y = rawData[:,2] #rowmin:rowmax, column 1 "21 cm power"
    axis[i,j].plot(x,y,color=colorName)

#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Program 
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

#PART 1: plot 21 cm power vs k at diff redshifts
    #create canvas (window) and pads (graph regions on the canvas)
dim = 3 #dimensions of pads in canvas 
fig,ax = plt.subplots(dim,dim)
    #choose redshifts 
z0 = "04.7077"
z1 = "07.0445"
z2 = "10.0859"
'''
#rawDemo = np.loadtxt(".//t21_power_data/t21_power_democratic/auto_t21_z=07.0445",skiprows=0,max_rows=1,usecols=(2,))
#callout(rawDemo)
Regarding the data ASCII files per z: 
1st row, 3 columns
(0,0) = 'k[unitless]'
(0,1) = 'k[h/Mpc]' 
(0,2) = 'D^2(k)' 
1st row are solely string elements, thus ignore 
(dimless k, dimful k, 21 cm power)
'''
    #democratic sources model
plot21cm_vs_k("democratic",z0,ax,"tab:blue",0,0)
plot21cm_vs_k("democratic",z1,ax,"tab:blue",1,0)
plot21cm_vs_k("democratic",z2,ax,"tab:blue",2,0)
    #fiducial sources model
plot21cm_vs_k("fiducial",z0,ax,"tab:green",0,1)
plot21cm_vs_k("fiducial",z1,ax,"tab:green",1,1)
plot21cm_vs_k("fiducial",z2,ax,"tab:green",2,1)
    #oligarchic sources model
plot21cm_vs_k("oligarchic",z0,ax,"tab:red",0,2)
plot21cm_vs_k("oligarchic",z1,ax,"tab:red",1,2)
plot21cm_vs_k("oligarchic",z2,ax,"tab:red",2,2)
    #set other labels for graph
fig.suptitle("21 cm Power vs k (h/Mpc)")
ax[0,0].set_title("Democratic")
ax[0,1].set_title("Fiducial")
ax[0,2].set_title("Oligarchic")
ax[0,0].set_ylabel("z="+z0)
ax[1,0].set_ylabel("z="+z1)
ax[2,0].set_ylabel("z="+z2)
#plt.show() 

#PART 2: plot 21 cm power vs redshift at diff k 
#use interpolate from scipy (to create function from data)
    #get all redshifts (all models have same redshifts)
os.chdir(".//t21_power_data/t21_power_fiducial/") #change to directory of the data
redshift = (os.listdir()) #listdir returns a list
redshift = [float(x.replace("auto_t21_z=","")) for x in redshift] #use comprehension to remove characters for loop 
redshift = np.array(redshift) #convert list to array
#callout(redshift)
os.chdir("../..") #return back to original directory 
    #get all available k values (all models and all redshifts have same k)
kNum = np.loadtxt(".//t21_power_data/t21_power_democratic/auto_t21_z="+z0,skiprows=1,usecols=(1,))
#callout(kNum)
    #choose a k
index_k = 13
kChosen = kNum[index_k] 
print("Chosen k:",kChosen)
    #democratic sources model

#power21cm = [np.loadtxt(".//t21_power_data/t21_power_democratic/auto_t21_z="+x,skiprows=1,usecols=(2,))[index_k] for x in redshift]
#callout(np.array(power21cm))

 

