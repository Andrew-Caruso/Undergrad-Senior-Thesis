import numpy as np
import os 
import math 
from itertools import product #cartesian product
from datetime import datetime 

#take in the 21cm fields and chop them up into smaller volumes
#and store them in another directory to be converted into power spectra
#for each redshift and each model

#each 21cm field file is a binary file with N=1024^3 32-bit floats
#NOTE: WARNING: the number of cells is  N=1024^3 while the box volume is 1000^3 (Mpc/h)^2
#which contains the 21cm brightness temperature at each position in a (reshaped) box

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

def loadIn3(directory,manifest):
    current = os.getcwd()
    os.chdir(directory)
    data = np.genfromtxt(manifest,skip_header=1,dtype='str',usecols=(0,2))
    #skip the second column, such that 1st column is density file ID and 3rd column is redshift 
    #print("data unrounded:\n",data) 
    #truncate the redshifts (2nd column) to 4th decimal place
    data[:,1] = [str(f'{float(data[i,1]):0.6f}')[:-2] for i in range(len(data[:,1]))] 
    #convert to float, round to 6th decimal place and return to string, then remove that last two characters to truncate  
    #need to truncate to 4th decimal place (NOT round because it will not make with z from getRS) 
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

#to get data from density file at chosen redshift 
def loadIn2(directory,density_list,record,z,N):
    kill = 1 
    rawdata = -1
    #get the corresponding density file ID
    for i, z_record in enumerate(record[:,1]):
        if float(z_record) == float(z):
            ID = record[i,0]   
            kill = 0
            break 
    if kill == 0:
        #print("ID:",ID)
        split_name = density_list[0].split('.')
        file_name = split_name[0]+"."+ID+"."+split_name[2]  
        file_path = directory+file_name 
        #print("file path:",file_path)
        rawdata = np.array(np.fromfile(file_path,dtype=np.float32).reshape((N,N,N)))
    return rawdata,kill 

def getRS(model_directory):
    current = os.getcwd() #get current directory
    os.chdir(model_directory)#change to directory of data
    rs_str = (os.listdir()) #return a list of all files in that directory
    rs_str = [(x.replace("t21_field.z=","")) for x in rs_str] #use comprehension to replace all text save for redshift with blank
    os.chdir(current)#return back to original directory
    #convert from str to float back to str, to remove extra 0s, e.g. 05.0297 to 5.0297 
    return rs_str

def getTerms(trans,N,step,N_max):
    terms = []
    first = int(0)
    second = int(N)
    if trans > 1: 
        for i in range(1,trans+1):
            if second <= N_max:
                phrase = str(first)+":"+str(second)
                terms.append(phrase)
                first += step
                second += step
    elif trans == 1:
        for i in range(0,trans+1):
            if second <= N_max:
                phrase = str(first)+":"+str(second)
                terms.append(phrase)
                first += step
                second += step
    return terms 

def getTrans(step,N,N_default):
    if N == N_default:
        step_trans = int(1)
    else:
        step_trans = (N_default-N)/step #number of valid transitions using the step
        step_trans = int(math.ceil(step_trans)) #round up 
    return step_trans


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

def loadOut(data,phrase,outputDir):
    filename = "t21_field."+phrase
    outputFile = outputDir + filename  
    print("CREATING %s"%(outputFile))
    data.tofile(outputFile) #save as binary format 


#max volume of box is 1 Gpc/h x 1 Gpc/h x 1 Gpc/h 
L_default = 1000 #Mpc/H or 1 Gpc physical CHANGE ME 
N_default = 1024 #cells CHANGE ME  
V_max = N_default**3
#list of chosen cuts of the 1 Gpc length to divide the cube 
L_chosen = [L_default,500] #Mpc/h CHANGE ME
#N_chosen = [N_default,102.4,204.8,307.2,409.6,512.0] #cells CHANGE ME
N_chosen = [N_default,512] #cells CHANGE ME
#NOTE: use L for labeling and N for computation; also the N_chosen values were rounded accordingly

#prelims 
emissivity_file = '/expanse/lustre/scratch/ccain002/temp_project/for_andrew/1GPCh_fields/' #directory of the 21cm field data to use. CHANGE ME 
outputData_dir = '/expanse/lustre/scratch/andrewcaruso/temp_project/final_powers/21cm_fields_500subs/' #directory of the output data CHANGE ME
density_dir = '/expanse/lustre/scratch/ccain002/temp_project/for_andrew/1GPCh_fields/density_outputs/' #directory of the corresponding density files for the 21cm field data CHANGE ME 
model_list = ["oligarchic/","extreme/"]
#model_list = ["extreme/"]
emissivity_models = ["",""]
emissivity_models = [emissivity_file+model_list[i] for i in range(2)]
z_str0 = getRS(emissivity_models[0])
z_str1 = getRS(emissivity_models[1])
z_arr_str = np.array([z_str0,z_str1])
print("z_arr_str:",z_arr_str)
current_dir = os.getcwd() #get current directory 
weight = 0.1 #weight on the std 
print("std weight:",weight)
numSteps = 2 #number of both valid and invalid steps depending on the cut size N
step = int(N_default/numSteps) #step size for moving the subvolume in the loop
print("step size:",step)
cosmicDensity = 1 #densities are in units of cosmic density which is defined as unity

#extract record of all density file IDs with each redshift from manifest file
density_list, manifest = listGrab(density_dir)
#print("manifest:",manifest) 
#print("density_list\n",density_list)
record = loadIn3(density_dir,manifest)  
#TEST: redshifts from input files matches redshifts in record from manifest
#print("record:\n",record)
#convert redshifts to correct format to match record
#test_z = np.array(z_str0)
#test_z = [str(f'{float(rs):0.6f}')[:-2] for rs in test_z] 
#print("redshifts:\n",test_z)
#missing last density file with ID 124, thus delete it from record and proceed
record = record.tolist()
del record[-1] 
record = np.array(record) 
#print("record:\n",record)
#confirm chosen directory for output exists 
dirCreator(outputData_dir)

#perform the cuts of the fields 
#loop over every model
for m, model in enumerate(emissivity_models):
    outputModelDir = outputData_dir + model_list[m] 
    dirCreator(outputModelDir)
    #loop over every redshift per model 
    for q, rs in enumerate(z_arr_str[m]):
        field = loadIn1(model,rs,N_default)
        outputRsDir = outputModelDir + "z=" + rs + "/"
        #change format of redshift to match the redshifts in record 
        rs = str(np.round(float(rs),decimals=5))
        rs = str(f'{float(rs):0.6f}')[:-2]
        #print("\nz:",rs)
        #get the density data corresponding to the particular redshift 
        #the density files are in units of cosmic density 
        density,kill = loadIn2(density_dir,density_list,record,rs,N_default)
        #handle any unmatched redshifts 
        if kill == 1:
            print("skipped z: ",rs)
            continue 
        dirCreator(outputRsDir)
        #loop over every chosen cut 
        #NOTE: only save the 500 cut not the default size 
        for f, N in enumerate(N_chosen):
            if N != N_default: 
                step_trans = getTrans(step,N,N_default)
                terms = getTerms(step_trans,N,step,N_default) #all terms based on number of transitions 
                perms = permutate(terms) #all permutations of the terms
                #create a new directory for all the new fields from the chosen cut 
                outputCutDir = outputRsDir + "cut_" + str(L_chosen[f]) + "/" 
                dirCreator(outputCutDir) 
                print("cut: ",L_chosen[f]) 

                #get the standard deviation of average densities from all subvolumes for particular cut
                avg_densityList = [] #to save all average densities to calculate STD
                #loop over every permutation of the terms
                for permTuple in perms:
                    permArray = np.array(list(permTuple)) 
                    CP = [phrase.split(':') for phrase in permArray] 
                    CP = [[int(CP[j][i]) for i in range(len(CP[j]))] for j in range(len(CP)) ]
                    dataSize = N**3  
                    longPhrase = str(CP[0][0])+":"+str(CP[0][1])+","+str(CP[1][0])+":"+str(CP[1][1])+","+str(CP[2][0])+":"+str(CP[2][1])
                    newDensity = (density[CP[0][0]:CP[0][1],CP[1][0]:CP[1][1],CP[2][0]:CP[2][1]]).reshape(dataSize,) 
                    avgDensity = np.average(newDensity)
                    avg_densityList.append(avgDensity)
                    if L_chosen[f] == 500:
                        print(avgDensity)
                        #to double check that it is not the same average each time for the 500 cut
                std = np.std(avg_densityList)
                newSTD = weight * std
                print("CosmicDensity"+u"\u00B1"+"STD:",cosmicDensity+newSTD," to ",cosmicDensity-newSTD)

                #get the field subvolumes and only keep those within the cosmic density +- weighed STD 
                saveCount = 0 #number of saved fields (that are within the cosmic density criterion)
                #loop over every permutation of the terms
                for permTuple in perms:
                    #3 column positions in the fields data  
                    permArray = np.array(list(permTuple)) #convert to array 
                    CP = [phrase.split(':') for phrase in permArray] #split at the commas 
                    CP = [[int(CP[j][i]) for i in range(len(CP[j]))] for j in range(len(CP)) ]
                    #create new fields for each cut and reshape it to raw data shape 
                    dataSize = N**3  
                    newfield = (field[CP[0][0]:CP[0][1],CP[1][0]:CP[1][1],CP[2][0]:CP[2][1]]).reshape(dataSize,) 
                    longPhrase = str(CP[0][0])+":"+str(CP[0][1])+","+str(CP[1][0])+":"+str(CP[1][1])+","+str(CP[2][0])+":"+str(CP[2][1])
                    newDensity = (density[CP[0][0]:CP[0][1],CP[1][0]:CP[1][1],CP[2][0]:CP[2][1]]).reshape(dataSize,) 
                    avgDensity = np.average(newDensity)
                    #NOTE: REMOVE THE DENSITY CRITERION FOR THE 500 MPC/h cuts (i.e. keep all cuts)!!!
                    '''
                    if avgDensity < cosmicDensity+newSTD and avgDensity > cosmicDensity-newSTD:
                        #does NOT save N=N_default since STD = 0 because only 1 permutation (full box size), only one average 
                        #print("Save:",longPhrase,"with avg density:",avgDensity) 
                        loadOut(newfield,longPhrase,outputCutDir) 
                        saveCount+=1
                    elif N == N_default:
                        #save the total box size 
                        #print("Save:",longPhrase,"with avg density:",avgDensity) 
                        loadOut(newfield,longPhrase,outputCutDir) 
                        saveCount+=1
                    '''
                    loadOut(newfield,longPhrase,outputCutDir) 
                    saveCount+=1
                print("For N=",N," Save count:",saveCount, " VS permutations:",len(perms))
                if saveCount == 0:
                    print("WARNING: saved ",saveCount," for Z=",rs," N=",N," with permutations: ",len(perms))
                
print("\nTASK COMPLETED\n") 

