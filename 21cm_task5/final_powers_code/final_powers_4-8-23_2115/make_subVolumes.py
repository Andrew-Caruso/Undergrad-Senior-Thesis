import numpy as np
import os 
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

def loadIn(model,rs,N):
    rsFile = model+"t21_field.z="+rs 
    type(np.fromfile(rsFile,dtype=np.float32))
    rawdata = np.array(np.fromfile(rsFile,dtype=np.float32).reshape((N,N,N)))
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
L_chosen = [L_default,500,250,125] #Mpc/h CHANGE ME
N_chosen = [N_default,512, 256, 128] #cells CHANGE ME
#NOTE: use L for labeling and N for computation 



#prelims 
emissivity_file = '/expanse/lustre/scratch/ccain002/temp_project/for_andrew/1GPCh_fields/' #directory of the 21cm field data to use. CHANGE ME 
outputData_dir = '/expanse/lustre/scratch/andrewcaruso/temp_project/final_powers/21cm_fields_subVolumes/' #directory of the output data CHANGE ME
model_list = ["oligarchic/","intermediate/","extreme/"]
#model_list = ["extreme/"]
emissivity_models = ["","",""]
emissivity_models = [emissivity_file+model_list[i] for i in range(3)]
z_str0 = getRS(emissivity_models[0])
z_str1 = getRS(emissivity_models[1])
z_str2 = getRS(emissivity_models[2])
z_arr_str = np.array([z_str0,z_str1,z_str2])
print("z_arr_str:",z_arr_str)

#confirm chosen directory for output exists 
dirCreator(outputData_dir)
#loop over every model
for m, model in enumerate(emissivity_models):
    outputModelDir = outputData_dir + model_list[m] 
    dirCreator(outputModelDir)
    #loop over every redshift per model 
    for rs in z_arr_str[m]:
        fields = loadIn(model,rs,N_default)
        outputRsDir = outputModelDir + "z=" + rs + "/"
        dirCreator(outputRsDir)
        #loop over every chosen cut 
        for i, N in enumerate(N_chosen):
            L = L_chosen[i] #use only for labeling 
            trans = int(N_default / N)  
            terms = getTerms(trans,N)
            perms = permutate(terms)
            #create a new directory for all the new fields from the chosen cut 
            outputCutDir = outputRsDir + "cut_" + str(L) + "/" 
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
print("\nTASK COMPLETED\n") 

