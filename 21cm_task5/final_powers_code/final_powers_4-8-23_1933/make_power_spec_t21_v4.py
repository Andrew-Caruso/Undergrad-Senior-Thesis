import numpy as np
import os, sys, glob
from os import listdir
from os.path import isfile, join
import time

def callout(bob):
    print("type:",type(bob)) #data type, should be numpy array
    print("dim:",bob.ndim) #dimensions
    print("shape:",bob.shape) #row by column
    print("size:",bob.size) #number of elements
    #print("\n")
    #print(bob,"\n")

def getRS(model_directory):
    current = os.getcwd() #get current directory 
    os.chdir(model_directory)#change to directory of data 
    rs_str = (os.listdir()) #return a list of all files in that directory
    rs_str = [(x.replace("z=","")) for x in rs_str] #use comprehension to replace all text save for redshift with blank
    os.chdir(current)#return back to original directory 
    rs = [float(rs_str[i]) for i in range(len(rs_str))]
    return rs, rs_str 

def getCuts(cut_dir,rs):
    current = os.getcwd() #get current directory 
    os.chdir(cut_dir)#change to directory of data 
    cut_str = (os.listdir())
    cut_str = [(x.replace("cut_","")) for x in cut_str]
    os.chdir(current)
    cutValue = [float(cut_str[i]) for i in range(len(cut_str))]
    return cutValue, cut_str 

def dirChecker(directory):
        #create directory if non-existent 
        if not os.path.isdir(directory): print("MAKING %s"%(directory)); os.mkdir(directory)
        else: print("%s EXISTS"%(directory))

def dataGrab(directory):
    current = os.getcwd() 
    os.chdir(directory)
    data_list = (os.listdir())
    os.chdir(current)
    cutRanges = [(x.replace("t21_field.","")) for x in data_list]
    return data_list, cutRanges  
    

emissivity_file = '/expanse/lustre/scratch/andrewcaruso/temp_project/final_powers/21cm_fields_subVolumes/' #directory of the sliced 21cm field data to use. CHANGE ME
model_list = ["oligarchic/","intermediate/","extreme/"]
emissivity_models = ["","",""]
emissivity_models = [emissivity_file+model_list[i] for i in range(3)] 
#print("\nmodels:",emissivity_models,"\n")
output_path = "/expanse/lustre/scratch/andrewcaruso/temp_project/final_powers/" #CHANGE ME
output_mainDir = "21cm_spectra_subVolumes" 
output_dir = output_path+output_mainDir 
dirChecker(output_dir)

#get the redshifts
z_arr0,z_str0 = getRS(emissivity_models[0])
z_arr1,z_str1 = getRS(emissivity_models[1])
z_arr2,z_str2 = getRS(emissivity_models[2])
z_arr = np.array([z_arr0,z_arr1,z_arr2])
z_arr_str = np.array([z_str0,z_str1,z_str2])
current_dir = os.getcwd() #get current directory 

#loop over models 
for m, model_dir in enumerate(emissivity_models):
    #get major input directory 
    z = z_arr[m]
    z_string = z_arr_str[m]
    #print("z_string:",z_string)
    chopped_dir = model_dir.split('/')
    model_name = chopped_dir[-2:-1]
    model_dir2 = '/'.join(chopped_dir[:-3])
    #print("model_dir: ", model_dir,"\nmodel_dir2: ",model_dir2,"\nmodel_name: ",model_name[0])
    output_dir1 = output_dir+"/"+model_name[0] 
    dirChecker(output_dir1) #check output main directory 
    interval_list = [] #create interval list for the redshift loop 

    #loop over redshifts per model
    for n, z_str in enumerate(z_string):
        start = time.time() #begin timer 
        cut_dir = model_dir+"z="+z_str
        #print("cut directory:",cut_dir)
        cuts, cuts_string = getCuts(cut_dir,z_str) 
        #print("cuts:",cuts_string)
        output_dir2 = output_dir1+"/z="+z_str
        dirChecker(output_dir2) #check output 1st sub directory 
        #loop over cuts per redshift 
        for j, cut_str in enumerate(cuts_string):
            model_dir3 = cut_dir+"/cut_"+cut_str  
            #print("model_dir3:",model_dir3) 
            data_list,cut_ranges = dataGrab(model_dir3) 
            #print("data_list",data_list) 
            #print("cut ranges:",cut_ranges)
            output_dir3 = output_dir2+"/cut_"+cut_str
            dirChecker(output_dir3) #check output 2nd sub directory 
            #NOTE: necessary change implemented below: 
            #cut_str = L or physical distance in Mpc/h but
            #the number of cells = N is what needs to be used for computation 
            #while cut_str L is used for labeling 
            #e.g. the 1000 Mpc box is really 1024 cells 
            num = str(1024/(1000/int(cut_str))) #CHANGE ME!  
            boxSize = num +"\n" #new box size code fortran code 
            numCells = num #new number of cells for fortran code
            #loop over every field per cut
            for n, data_name in enumerate(data_list):
                #two input data files for FFT 
                dT_FILE_1 ="'"+model_dir3+"/"+data_name+"'\n"
                dT_FILE_2 ="'"+model_dir3+"/"+data_name+"'\n"
                #form output directory for the data file
                outfile = "'%s/auto_21cm_=%s'\n"%(output_dir3,cut_ranges[n])
                #print("dT_FILE_1:",dT_FILE_1)
                #print("dT_FILE_2:",dT_FILE_2)
                #print("outfile:",outfile)
                #open the fortran code and edit it
                filename = "mainbody_t21.f90"
                with open(filename, 'r') as f: 
                    lines = f.readlines()
                #loop over every line in fortran code 
                for i in range(len(lines)): 
                    if lines[i].find("dir     = ") != -1: 
                        oldfile = lines[i].split('= ')[1]
                        lines[i] = lines[i].replace(oldfile, "%s\n"%(output_dir))
                    if lines[i].find("outfile = ") != -1: 
                        oldfile = lines[i].split('= ')[1]
                        lines[i] = lines[i].replace(oldfile, outfile)
                    if lines[i].find("datafile1 = ") != -1: 
                        oldin1 = lines[i].split('= ')[1]
                        lines[i] = lines[i].replace(oldin1, dT_FILE_1)
                    if lines[i].find("datafile2 = ") != -1:         
                        oldin2 = lines[i].split('= ')[1]
                        lines[i] = lines[i].replace(oldin2, dT_FILE_2)	
                    if lines[i].find("integer, parameter :: n ") != -1: 
                        newline1 = lines[i].split('= ')[1]
                        oldnum = newline1.split(', ')[0]
                        #print("oldnum:",oldnum)
                        lines[i] = lines[i].replace(oldnum, numCells)
                    if lines[i].find("real(sp), parameter :: boxsize ") != -1: 
                        oldbox = lines[i].split('= ')[1]
                        #print("oldbox:",oldbox)
                        lines[i] = lines[i].replace(oldbox, boxSize)
                #over-write the fortran code
                with open(filename, 'w+') as f1: 
                    f1.writelines(lines)
                #access the FFT in home and create executable file, thus executing the fortran code 
                os.system('ifort -O3 -qopenmp -I/home/andrewcaruso/fftw-3.3.9/include/ nrtype.f90 powers_t21.f90 mainbody_t21.f90  -L/home/andrewcaruso/fftw-3.3.9/lib/ -lfftw3_threads -lfftw3 -lm -o powers_t21.x')
                #execute the file
                os.system('./powers_t21.x')
        #in redshift loop
        end = time.time() #end time  
        interval = end-start 
        print("Time to finish redshift:",interval)
        interval_list.append(interval)
    #in model loop
    if len(interval_list) != 0:
        avg_interval = sum(interval_list)/len(interval_list)
        print("Average redshift computation time:",avg_interval)

print("\n \n ")
