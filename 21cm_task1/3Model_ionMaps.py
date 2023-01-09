from json import load
import matplotlib.pyplot as plt #for plotting
from scipy import interpolate 
import numpy as np 
import os 
#to fix bug in part 1 (overlapping minor tick labels) do:
from matplotlib import ticker
formatter1 = ticker.ScalarFormatter(useMathText=True) #create object to scalar formatter class with fancy math text
formatter1.set_scientific(True)
formatter2 = ticker.LogFormatterSciNotation(labelOnlyBase=True)
#now set all x/y axes and their major/minor ticks to use this formatter
#https://atmamani.github.io/cheatsheets/matplotlib/matplotlib_2/#Placement-of-ticks-and-custom-tick-labels 

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

def getRandP(model,index_k):
    #get the redshift for the model
    os.chdir(".//t21_power_data/t21_power_"+model+"/") #change to directory of the data
    rs = (os.listdir()) #listdir returns a list: redshift list 
    rs = [(x.replace("auto_t21_z=","")) for x in rs] #use comprehension to remove characters for loop 
    rs = np.array(rs) #convert list to array
    #callout(rs)
    os.chdir("../..") #return back to original directory 
    #get the powers for the model at every redshift corresponding to the chosen fixed k
    power = [np.loadtxt(".//t21_power_data/t21_power_"+model+"/auto_t21_z="+x,skiprows=1,usecols=(2,))[index_k] for x in rs] #use comprehension to select out powers for the selected k (index) for every redshift
    #i.e. go to every redshift, select the power corresponding to the chosen fixed k
    power = np.array(power) #convert list to array
    return rs, power 

def minMax(model,redshift):
    print(model+" redshifts:\n")
    callout(redshift)
    print(model+" min:",min(redshift))
    print(model+" max:",max(redshift),"\n")

def plot21cm_vs_z(f1,rs,ax,colorName,i,j):
    new_rs = np.linspace(min(rs),max(rs),len(rs)) #create new redshifts array for the graph
    new_power = f1(new_rs) #get the new 21cm power output from the function (previously built from interpolation)
    ax[i,j].plot(new_rs,new_power,color=colorName)

#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Program 
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

#PART 1: plot 21 cm power vs k at diff redshifts
print("PART 1: 21cm Power vs k at different Redshifts")
    #create canvas (window) and pads (graph regions on the canvas)
dim = 3 #dimensions of pads in canvas 
fig,ax = plt.subplots(dim,dim)
    #choose redshifts 
z0 = "04.7077"
z1 = "07.0445"
z2 = "10.0859"
print("Chosen redshifts: ",z0," ",z1," ",z2,"\n")
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
#k of interest: ~0.05 to ~0.5 
#use nested comprehension in lieu of nest for loops
[[ax[i,j].set_xlim([0.05,0.5]) for j in range(3)] for i in range(3) ]
#use log scale
[[ax[i,j].set_xscale("log") for j in range(3)] for i in range(3) ]
[[ax[i,j].set_yscale("log") for j in range(3)] for i in range(3) ]
#fix x axis bug of overlapping minor tick labels 
#first solution: make the x-axis minor tick labels super small
#[[ax[i,j].tick_params(axis='x',which='both',labelsize='xx-small') for j in range(3)] for i in range(3) ]
#second solution: changing x-axis minor tick labels to scientific notation
[[ax[i,j].xaxis.set_minor_formatter(formatter2) for j in range(3)] for i in range(3) ]
#plt.show() 


#PART 2: plot 21 cm power vs redshift at diff k 
print("PART 2: 21cm Power vs Redshift at different k")
    #get all available k values (same for all models at all redshifts k)
kNum = np.loadtxt(".//t21_power_data/t21_power_democratic/auto_t21_z="+z0,skiprows=1,usecols=(1,))
print("\nChoices of k:")
callout(kNum)
    #choose wavenumbers k
index_k = [1,5,24] #3 indices from 0 to 149
kChosen = np.array([kNum[x] for x in index_k]) #append k values corresponding to chosen k and convert list to array
print("Chosen k:",kChosen,"\n")
    #create canvas/figure/window with 9 pads/axes/plots
fig2,ax2 = plt.subplots(3,3)
    #get all redshifts (not same for all models) and corresponding 21cm powers
    #then interpolate original data to create function: power(redshift) at chosen k
    #then using the function from interpolation and new redshift domain, plot the graph
    #democratic sources model
rsD, powerD0 = getRandP("democratic",index_k[0])
powerD1 = getRandP("democratic",index_k[1])[1]
powerD2 = getRandP("democratic",index_k[2])[1]
rsD = rsD.astype(float) #convert redshift string array to float array
f1_D0 = interpolate.interp1d(rsD,powerD0)
f1_D1 = interpolate.interp1d(rsD,powerD1)
f1_D2 = interpolate.interp1d(rsD,powerD2)
#minMax("democratic",rsD)
plot21cm_vs_z(f1_D0,rsD,ax2,"tab:blue",0,0)
plot21cm_vs_z(f1_D1,rsD,ax2,"tab:blue",1,0)
plot21cm_vs_z(f1_D2,rsD,ax2,"tab:blue",2,0)
    #fiducial sources model 
rsF, powerF0 = getRandP("fiducial",index_k[0])
powerF1 = getRandP("fiducial",index_k[1])[1]
powerF2 = getRandP("fiducial",index_k[2])[1]
rsF = rsF.astype(float) #convert redshift string array to float array
f1_F0 = interpolate.interp1d(rsF,powerF0)
f1_F1 = interpolate.interp1d(rsF,powerF1)
f1_F2 = interpolate.interp1d(rsF,powerF2)
#minMax("fiducial",rsF)
plot21cm_vs_z(f1_F0,rsF,ax2,"tab:green",0,1)
plot21cm_vs_z(f1_F1,rsF,ax2,"tab:green",1,1)
plot21cm_vs_z(f1_F2,rsF,ax2,"tab:green",2,1)
    #oligarchic sources model
rsO,powerO0 = getRandP("oligarchic",index_k[0])
powerO1 = getRandP("oligarchic",index_k[1])[1]
powerO2 = getRandP("oligarchic",index_k[2])[1]
rsO = rsO.astype(float) #convert redshift string array to float array
f1_O0 = interpolate.interp1d(rsO,powerO0)
f1_O1 = interpolate.interp1d(rsO,powerO1)
f1_O2 = interpolate.interp1d(rsO,powerO2)
#minMax("oligarchic",rsO)
plot21cm_vs_z(f1_O0,rsO,ax2,"tab:red",0,2)
plot21cm_vs_z(f1_O1,rsO,ax2,"tab:red",1,2)
plot21cm_vs_z(f1_O2,rsO,ax2,"tab:red",2,2)
#print("comparison:\n",rsO==rsF) #proof that fiducial model and oligarchic model redshifts differ
    #set labels for plot
fig2.suptitle("21 cm Power vs Redshift")
ax2[0,0].set_title("Democratic")
ax2[0,1].set_title("Fiducial")
ax2[0,2].set_title("Oligarchic")
ax2[0,0].set_ylabel("k="+str(kChosen[0]))
ax2[1,0].set_ylabel("k="+str(kChosen[1]))
ax2[2,0].set_ylabel("k="+str(kChosen[2]))
#change from linear scale to log scale
[[ax2[i,j].set_xscale("log") for j in range(3)] for i in range(3) ]
[[ax2[i,j].set_yscale("log") for j in range(3)] for i in range(3) ]
plt.show()

