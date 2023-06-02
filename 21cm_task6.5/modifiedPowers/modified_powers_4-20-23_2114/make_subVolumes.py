import numpy as np
import os 
from itertools import product #cartesian product
from datetime import datetime 

#take in the 21cm fields and chop them up into smaller volumes
#and store them in another directory to be converted into power spectra
#for each redshift and each model

#each 21cm field file is a binary file with N=300^3 32-bit floats
#which contains the 21cm brightness temperature at each position in a (reshaped) 300^3 box

def callout(name,bob):
    print(name)
    print("type:",type(bob)) #data type, should be numpy array  
    print("dim:",bob.ndim) #dimensions 
    print("shape:",bob.shape) #row by column 
    print("size:",bob.size,"\n") #number of elements 

def loadIn1(model,rs,N):
    rsFile = model+"t21_field.z="+rs 
    type(np.fromfile(rsFile,dtype=np.float32))
    rawdata = np.array(np.fromfile(rsFile,dtype=np.float32).reshape((N,N,N)))
    return rawdata 

def loadIn2(model,rs,N):
    rsFile = model+"t21_field.z="+rs 
    type(np.fromfile(rsFile,dtype=np.float32))
    rawdata = np.array(np.fromfile(rsFile,dtype=np.float32).reshape((N,N,N)))
    return rawdata 

def listGrab(directory):
    current = os.getcwd() 
    os.chdir(directory)
    data_list = (os.listdir()) #returns a list
    manifest = os.popen('ls *txt').read() #to output the actual commandline 
    manifest = manifest.strip() #remove extra space 
    os.chdir(current) 
    data_list.remove(manifest)  
    return data_list, manifest  

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

def loadOut(data,phrase,outputDir):
    filename = "t21_field."+phrase
    outputFile = outputDir + filename  
    print("CREATING %s"%(outputFile))
    data.tofile(outputFile) #save as binary format 

#max volume of box is 300 x 300 x 300
N_default = 300 #Mpc/H  
V_max = N_default**3
#list of chosen cuts of the 300 Mpc length to divide the cube 
N_chosen = [N_default,75,100,150] #Mpc/h CHANGE ME


#prelims 
emissivity_file = '/expanse/lustre/scratch/ccain002/temp_project/for_andrew/21cm_fields/' #directory of the 21cm field data to use. CHANGE ME
outputData_dir = '/expanse/lustre/scratch/andrewcaruso/temp_project/modified_powers/21cm_fields_subVolumes/' #directory of the output data CHANGE ME
density_dir = '/expanse/lustre/scratch/ccain002/temp_project/for_andrew/real_300_mod/' #directory of the corresponding density files for the 21cm field data CHANGE ME 
model_list = ["democratic/","fiducial/","oligarchic/"]
emissivity_models = ["","",""]
emissivity_models = [emissivity_file+model_list[i] for i in range(3)]
#REMOVE ME
'''
z_str0 = getRS(emissivity_models[0])
z_str1 = getRS(emissivity_models[1])
z_str2 = getRS(emissivity_models[2])
z_arr_str = np.array([z_str0,z_str1,z_str2])
'''
current_dir = os.getcwd() #get current directory 

#plot both density and corresponding field as a contour map
#to confirm they are being read in in the correct orientation 
#test using z = 12, wherein both density and 21cm field look the same 
density_list, manifest = listGrab(density_dir)
print("manifest:",manifest) 
print("density_list\n",density_list)
#need to extra data 


#get standard deviation of the fields 
#how to import the data??? 


#moving the box across the total volume but select the cuts with average densities that fall in the c*std of the cosmic density 
#the density files are in units of cosmic density 


'''
#perform the cuts of the fields 
#confirm chosen directory for output exists 
dirCreator(outputData_dir)
#start loop to perform cuts 
#loop over every model
for m, model in enumerate(emissivity_models):
    outputModelDir = outputData_dir + model_list[m] 
    dirCreator(outputModelDir)
    #loop over every redshift per model 
    for rs in z_arr_str[m]:
        fields = loadIn1(model,rs,N_default)
        outputRsDir = outputModelDir + "z=" + rs + "/"
        dirCreator(outputRsDir)
        #loop over every chosen cut 
        for N in N_chosen:
            trans = int(N_default / N)  
            terms = getTerms(trans,N)
            perms = permutate(terms)
            #create a new directory for all the new fields from the chosen cut 
            outputCutDir = outputRsDir + "cut_" + str(N) + "/" 
            dirCreator(outputCutDir) 
            #loop over every permutation of the terms
            for permTuple in perms:
                #3 column positions in the fields data  
                permArray = np.array(list(permTuple)) #convert to array 
                CP = [phrase.split(':') for phrase in permArray] #split at the commas 
                CP = [[int(CP[j][i]) for i in range(len(CP[j]))] for j in range(len(CP)) ]
                #create new fields for each cut and reshape it to raw data shape 
                dataSize = N**3  
                newfields = (fields[CP[0][0]:CP[0][1],CP[1][0]:CP[1][1],CP[2][0]:CP[2][1]]).reshape(dataSize,) 
                longPhrase = str(CP[0][0])+":"+str(CP[0][1])+","+str(CP[1][0])+":"+str(CP[1][1])+","+str(CP[2][0])+":"+str(CP[2][1])
                loadOut(newfields,longPhrase,outputCutDir) 
'''
print("\nTASK COMPLETED\n") 

