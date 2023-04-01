import numpy as np
import os, sys, glob
from os import listdir
from os.path import isfile, join

def callout(bob):
    print("type:",type(bob)) #data type, should be numpy array
    print("dim:",bob.ndim) #dimensions
    print("shape:",bob.shape) #row by column
    print("size:",bob.size) #number of elements
    #print("\n")
    #print(bob,"\n")

emissivity_file = '/expanse/lustre/scratch/ccain002/temp_project/for_andrew/21cm_fields/' #directory of the 21cm field data to use. CHANGE ME
model_list = ["democratic/","fiducial/","oligarchic/"]
emissivity_models = ["","",""]
emissivity_models = [emissivity_file+model_list[i] for i in range(3)] 
print("\nmodels:",emissivity_models,"\n")

def getRS(model_directory):
    current = os.getcwd() #get current directory 
    os.chdir(model_directory)#change to directory of data 
    rs_str = (os.listdir()) #return a list of all files in that directory
    rs_str = [(x.replace("t21_field.z=","")) for x in rs_str] #use comprehension to replace all text save for redshift with blank
    os.chdir(current)#return back to original directory 
    rs = [float(rs_str[i]) for i in range(len(rs_str))]
    return rs, rs_str 

#get the redshifts
z_arr0,z_str0 = getRS(emissivity_models[0])
z_arr1,z_str1 = getRS(emissivity_models[1])
z_arr2,z_str2 = getRS(emissivity_models[2])
z_arr = np.array([z_arr0,z_arr1,z_arr2])
z_arr_str = np.array([z_str0,z_str1,z_str2])
current_dir = os.getcwd() #get current directory 

for m, model in enumerate(emissivity_models):
    #print("m:",m) 
    z = z_arr[m]
    z_string = z_arr_str[m]
    print("z_string:",z_string)
    model_dir = emissivity_models[m]
    #model = model.split(',')
    #model_dir, model_name = model, model.split('/')[:-1] 
    #home_dir = model_dir.split('/')[:-1]
    chopped_dir = model_dir.split('/')
    model_name = chopped_dir[-2:-1]
    model_dir2 = '/'.join(chopped_dir[:-2])
    #print("model_dir: ", model_dir,"\nmodel_dir2: ",model_dir2,"\nmodel_name: ",model_name)
    output_dir = model_dir2+"/testing/"+model_name[0]+"_test"
    #print("output_dir:",output_dir,"\n")

    #create directory if non-existent 
    if not os.path.isdir(output_dir): print("MAKING %s"%(output_dir)); os.mkdir(output_dir)
    else: print("%s EXISTS"%(output_dir))
     
    #loop over redshifts per model
    #i = 0 #select index of z_string 
    for i, z_str in enumerate(z_string):
        dT_FILE_1 ="'"+model_dir+"t21_field.z="+z_str+"'\n"
        dT_FILE_2 ="'"+model_dir+"t21_field.z="+z_str+"'\n"
        outfile = "'%s/auto_21cm_z=%s'\n"%(output_dir, z_str)
        #print("dT_FILE_1:",dT_FILE_1)
        #print("dT_FILE_2:",dT_FILE_2)
        #print("outfile:",outfile)

        #open the fortran code and edit it
        filename = "mainbody_t21.f90"
        with open(filename, 'r') as f: 
            lines = f.readlines()

        for i in range(len(lines)): 
            if lines[i].find("dir     = ") != -1: 
                oldfile = lines[i].split('= ')[1]
                lines[i] = lines[i].replace(oldfile, "%s\n"%(output_dir))
            if lines[i].find("outfile = ") != -1: 
                oldfile = lines[i].split('= ')[1]
                #print("REPLACING %s with %s"%(oldfile, outfile))
                lines[i] = lines[i].replace(oldfile, outfile)
            if lines[i].find("datafile1 = ") != -1: 
                oldin1 = lines[i].split('= ')[1]
                #print("REPLACING %s with %s"%(oldin1, dT_FILE_1))
                lines[i] = lines[i].replace(oldin1, dT_FILE_1)
            if lines[i].find("datafile2 = ") != -1:         
                oldin2 = lines[i].split('= ')[1]
                #print("REPLACING %s with %s"%(oldin2, dT_FILE_2))
                lines[i] = lines[i].replace(oldin2, dT_FILE_2)	

        #over-write the fortran code
        with open(filename, 'w+') as f1: 
            f1.writelines(lines)
        
        #access the FFT in home and create executable file, thus executing the fortran code 
        os.system('ifort -O3 -qopenmp -I/home/andrewcaruso/fftw-3.3.9/include/ nrtype.f90 powers_t21.f90 mainbody_t21.f90  -L/home/andrewcaruso/fftw-3.3.9/lib/ -lfftw3_threads -lfftw3 -lm -o powers_t21.x')
        #execute the file
        os.system('./powers_t21.x')
print("\n \n \n ")
