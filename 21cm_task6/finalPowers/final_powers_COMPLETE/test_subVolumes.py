import numpy as np
import os 
from itertools import product #cartesian product
from datetime import datetime 
import matplotlib.pyplot as plt
from matplotlib.pyplot import contourf 
import math

#compare the orientations of the density field and 21cm brightness temperature data by creating contour maps at z=12. If the contour maps appear the same then the orientation matches. Otherwise, may need to swap the x,y, or z axes to get correct orientation. 

#each 21cm field file is a binary file with N=300^3 32-bit floats
#which contains the 21cm brightness temperature at each position in a (reshaped) 300^3 box
#NOTE: number of cells is 300 and physical length of box is 300 Mpc/h

def callout(name,bob):
    print(name)
    print("type:",type(bob)) #data type, should be numpy array  
    print("dim:",bob.ndim) #dimensions 
    print("shape:",bob.shape) #row by column 
    print("size:",bob.size,"\n") #number of elements 

def loadIn1(model,rs,N):
    rsFile = model+"t21_field.z="+rs 
    #import by default in C order, alternatively can do F order
    #to get correct orientation of the field to density plots
    rawdata = np.array(np.fromfile(rsFile,dtype=np.float32).reshape((N,N,N),order='F'))
    return rawdata 

def loadIn2(directory,manifest):
    current = os.getcwd()
    os.chdir(directory)
    data = np.genfromtxt(manifest,skip_header=1,dtype='str',usecols=(0,2))
    #skip the second column, such that 1st column is density file ID and 3rd column is redshift 
    data[:,1] = [str(float(f'{float(data[i,1]):0.4f}')) for i in range(len(data[:,1]))] 
    #convert to float, truncate to 4th decemical place and return to string using string float method
    #alternative method to rounding (yields same result)
    #data[:,1] = (np.round(data[:,1].astype(float),decimals=4)).astype(str)
    os.chdir(current) 
    return data 

def listGrab(directory):
    current = os.getcwd() 
    os.chdir(directory)
    data_list = (os.listdir()) #returns a list
    manifest = os.popen('ls *txt').read() #to output the actual commandline 
    manifest = manifest.strip() #remove extra space 
    os.chdir(current) 
    data_list.remove(manifest)  
    return data_list, manifest  

def getDensityFile(directory,density_list,record,z,N):
    #get the corresponding density file ID
    for i, z_record in enumerate(record[:,1]):
        if float(z_record) == float(z):
            ID = record[i,0]   
            break 
    #print("ID:",ID)
    split_name = density_list[0].split('.')
    file_name = split_name[0]+"."+ID+"."+split_name[2]  
    file_path = directory+file_name 
    print("file path:",file_path)
    rawdata = np.array(np.fromfile(file_path,dtype=np.float32).reshape((N,N,N)))
    rawdata = np.log10(rawdata) #convert density to log base 10 (for log scale)
    return rawdata 

def getRS(model_directory):
    current = os.getcwd() #get current directory
    os.chdir(model_directory)#change to directory of data
    rs_str = (os.listdir()) #return a list of all files in that directory
    rs_str = [(x.replace("t21_field.z=","")) for x in rs_str] #use comprehension to replace all text save for redshift with blank
    os.chdir(current)#return back to original directory
    return rs_str

def getTerms(trans,N):
	terms = []
	now = int(0)
	for i in range(1,trans+1):
		next = int(N*i)
		phrase = str(now)+":"+str(next)
		terms.append(phrase)
		now = next
	return terms 

def permutate(terms):
	positions = 3 
	numTerms = len(terms)
	numPerms = numTerms**positions
	#need to iterate through every permutation with repeats
	perms = list(product(terms,repeat=3))
	return perms  

def dirCreator(output_dir):
    if not os.path.isdir(output_dir): print("MAKING %s"%(output_dir)); os.mkdir(output_dir)
    else: print("%s EXISTS"%(output_dir))
    return 0

def loadOut(data,phrase,outputDir):
    filename = "t21_field."+phrase
    outputFile = outputDir + filename  
    print("CREATING %s"%(outputFile))
    data.tofile(outputFile) #save as binary format 
    return 0 

def redshiftSelector(wanted_Z,all_Z):
    #wanted_z must be a list of integers or floats
    chosen_Z = []
    if isinstance(all_Z,np.ndarray):
        all_Z = all_Z.tolist()
    #print("original:",all_Z,"\n")
    sorted_all_Z = all_Z.copy()
    sorted_all_Z.sort()
    #print("sorted:",sorted_all_Z,"\n")
    for i, wantZ in enumerate(wanted_Z):
        #print("wanted Z:",wantZ,type(wantZ))
        isInt = isinstance(wantZ,int) #determine if integer or not 
        #print("isInt:",isInt) 
        if isInt == 1:
            int_all_Z = [int(round(float(sorted_all_Z[i]),0)) for i in range(len(sorted_all_Z))]
            #print("sorted int:",int_all_Z,"\n")
            for j, Z in enumerate(int_all_Z):
                if wantZ == Z:
                    found_Z = sorted_all_Z[j]
                    #print("found Z:",found_Z,type(found_Z))
                    #print("sorted index:",j)
                    true_index = all_Z.index(found_Z)
                    #print("true index:",true_index)
                    #print("len of all_Z:", len(all_Z))
                    if j < len(all_Z)-1:
                        if (abs(float(found_Z)-float(wantZ)) < abs(float(sorted_all_Z[j+1])-float(wantZ))):
                            chosen_Z.append(found_Z)
                            break
                    elif j == len(all_Z)-1:
                            chosen_Z.append(found_Z)
                            break
                    elif (abs(float(found_Z)-float(wantZ))) > abs(float(sorted_all_Z[j+1])-float(wantZ)):
                            continue 
                    else:
                        break 
        elif isInt == 0:
            splitZ = str(wantZ).split(".")
            num = len(splitZ[1]) #number of decimal points to round to
            #print(num) 
            dec_all_Z = [round((float(sorted_all_Z[i])),num) for i in range(len(sorted_all_Z))]
            #print("sorted float:",dec_all_Z,"\n")
            for j, Z in enumerate(dec_all_Z):
                if wantZ == Z:
                    found_Z = sorted_all_Z[j]
                    #print("found Z:",sorted_all_Z[j])
                    #print("sorted index:",j)
                    true_index = all_Z.index(found_Z)
                    #print("true index:",true_index)
                    #print("len of all_Z:", len(all_Z))
                    if j < len(all_Z)-1:
                        if (abs(float(found_Z)-float(wantZ)) < abs(float(sorted_all_Z[j+1])-float(wantZ))):
                            chosen_Z.append(found_Z)
                            break
                    elif j == len(all_Z)-1:
                            chosen_Z.append(found_Z)
                            break
                    elif (abs(float(found_Z)-float(wantZ))) > abs(float(sorted_all_Z[j+1])-float(wantZ)):
                            continue 
                    else:
                        break 
    return chosen_Z

def contourPlotter(directory,density_list,record,model,rs_d,rs_f,N):
    field = np.array(loadIn1(model,rs_f,N)) #get 21cm field 
    density = getDensityFile(directory,density_list,record,rs_d,N) #get corresponding density 
    z_chosen = 10 #index value or distance along z axis in Mpc/h
    print("Chosen z axis distance:",z_chosen,"Mpc/h") #~7
    #create the xy grid
    x = np.linspace(0,N,N)
    y = np.linspace(0,N,N)
    xy_f = field[:,:,z_chosen] #z input to the contourf 
    xy_d = density[:,:,z_chosen] #z input to the contourf
    #find max and min values for color bar
    max_array = np.array([xy_f.max(),xy_d.max()])
    #max_xy = max_array.max() NEGLECT because two want two different color bars 
    min_array = np.array([xy_f.min(),xy_d.min()])
    #print("max field and density:",max_array)
    #print("min field and density:",min_array)
    return x,y,xy_d,xy_f,min_array,max_array,z_chosen

    

#max volume of box is 300 x 300 x 300
L_default = 1000 #Mpc/h, use for labels
N_default = 1024 #cells, use for calculations  
V_max = N_default**3
#list of chosen cuts of the 300 Mpc length to divide the cube 
L_chosen = [L_default,500,250,125] #Mpc/h CHANGE ME
N_chosen = [N_default,512, 256, 128] #cells CHANGE ME
#NOTE: use L for labeling and N for computation 


#prelims 
emissivity_file = '/expanse/lustre/scratch/ccain002/temp_project/for_andrew/1GPCh_fields/' #directory of the 21cm field data to use. CHANGE ME
outputData_dir = '/expanse/lustre/scratch/andrewcaruso/temp_project/final_powers/plots_densityVSfield_test/' #directory of the output data CHANGE ME
density_dir = '/expanse/lustre/scratch/ccain002/temp_project/for_andrew/1GPCh_fields/density_outputs/' #directory of the corresponding density files for the 21cm field data CHANGE ME 
model_list = ["oligarchic/","extreme/"]
emissivity_models = ["",""]
emissivity_models = [emissivity_file+model_list[i] for i in range(2)]
#REMOVE ME
z_str0 = getRS(emissivity_models[0])
z_str1 = getRS(emissivity_models[1])
z_arr_str = np.array([z_str0,z_str1])
print("z_arr_str",z_arr_str)
current_dir = os.getcwd() #get current directory 
#get the manifest list of density file ID and redshifts  
density_list, manifest = listGrab(density_dir)
record = loadIn2(density_dir,manifest)  
#chose a redshift to plot 
wanted_z = [11] #needs to be a list
#get the corresponding redshift in the manifest 
z_density = redshiftSelector(wanted_z,record[:,1])[0]
#confirm chosen directory for output exists 
dirCreator(outputData_dir)
print("wanted z:",wanted_z[0])
print("z_density:",z_density) 

#loop over every model
for m, model in enumerate(emissivity_models):
    outputModelDir = outputData_dir + model_list[m] 
    model_name = model_list[m].replace("/","")
    print("model:",model_name)
    dirCreator(outputModelDir)
    z_field = redshiftSelector(wanted_z,z_arr_str[m])[0] #convert from list to str  
    #print("z_field:",z_field) 
    outputRSdir = outputModelDir + "z=" + str(wanted_z[0])+ "/"
    dirCreator(outputRSdir)
    #create figure and get x,y,xy (or z) values for each model for each of the redshifts
    fig, ax = plt.subplots(1,2)
    x,y,xy_d,xy_f,min_array,max_array,pos = contourPlotter(density_dir,density_list,record,model,z_density,z_field,N_default)
    con_f = ax[0].contourf(x,y,xy_f,cmap="jet",vmin=min_array[0],vmax=max_array[0])
    con_d = ax[1].contourf(x,y,xy_d,cmap="jet",vmin=min_array[1],vmax=max_array[1])
    #set aspect ratio of x and y axis to be the same
    [ax[i].set_aspect('equal') for i in range(len(ax))] 
    #create separate color bars for field and density subplots
    num_f = plt.cm.colors.Normalize(vmin=min_array[0], vmax=max_array[0])
    num_d = plt.cm.colors.Normalize(vmin=min_array[1], vmax=max_array[1])
    sm_f = plt.cm.ScalarMappable(norm=num_f,cmap='jet')
    sm_d = plt.cm.ScalarMappable(norm=num_d,cmap='jet')
    ax_f = fig.add_axes([0.11,0.15,0.4,0.02]) #xmin,ymin,width,height
    ax_d = fig.add_axes([0.54,0.15,0.4,0.02]) #xmin,ymin,width,height
    tickList_f = list(range(int(min_array[0]),int(math.floor(max_array[0])),5)) #select major ticks for colorbar 
    tickList_d = list(np.arange(min_array[1],max_array[1],0.2)) #select major ticks for colorbar 
    #print("density ticks:",tickList_d)
    cbar_f = fig.colorbar(sm_f,orientation='horizontal',cax=ax_f,label="mK",ticks=tickList_f)
    cbar_d = fig.colorbar(sm_d,orientation='horizontal',cax=ax_d,label="mK",ticks=tickList_d)
    #change figure size 
    fig.set_size_inches(11,8)
    #add titles
    fig.suptitle("21cm Brightness Temperature (mK) vs Density for fixed position "+str(pos)+" Mpc/h")
        #figure coordinates (bottom left corner is (0,0) and top right corner is (1,1)
        #https://stackoverflow.com/questions/55767312/how-to-position-suptitle
    ax[0].set_xlabel("Mpc/h")
    ax[1].set_xlabel("Mpc/h")
    ax[0].set_ylabel("Mpc/h")
    ax[1].set_ylabel("Mpc/h")
    ax[0].set_title("21cm Field at z="+z_field)
    ax[1].set_title("Density at z="+z_density)
    current = os.getcwd()
    os.chdir(outputRSdir)
    plt.show()
    #save and output the plot 
    outfile = "contourMap_"+model_name+".pdf"
    fig.savefig(outfile,format='pdf')
    print("\noutput: ",outfile,"\n")
    os.chdir(current) 

print("\nTASK COMPLETED\n") 

