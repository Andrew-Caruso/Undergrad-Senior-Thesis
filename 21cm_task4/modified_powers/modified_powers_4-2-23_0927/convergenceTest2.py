import numpy as np
import os
import time 
import matplotlib.pyplot as plt
from operator import itemgetter

#create power P vs redshift z plots for fixed wavenumber k  

#----------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------
#functions 
#----------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------

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

def listGrab(directory):
    current = os.getcwd() 
    os.chdir(directory)
    data_list = (os.listdir())
    os.chdir(current)
    return data_list

def getDefaultCut(cuts_float):
    maxCut = max(cuts_float)
    cuts_float.remove(maxCut)
    maxCut = str(int(maxCut))
    cuts_str = [str(int(cuts_float[i])) for i in range(len(cuts_float))]
    return maxCut, cuts_str

def loadIn(filepath):
    #the power spectra are ASCII text files 
    rawData = np.loadtxt(filepath,skiprows=1)
    #row are strings of the titles (not actual data)
    #column 0 is k[unitless]
    #column 1 is k[h/Mpc]
    #column2 is D^2(k) or power spectrum at the corresponding k value 
    return rawData 

def loadIn2(model,rs,N):
    rsFile = model+"t21_field.z="+rs 
    type(np.fromfile(rsFile,dtype=np.float32))
    rawdata = np.array(np.fromfile(rsFile,dtype=np.float32).reshape((N,N,N)))
    return rawdata 

def isIntChecker(x):
    num = int(abs(float(int(x) - x)*10)) #tenth place of number 
    if num == 0:
        isInt = 1
    else:
        isInt = 0
    return isInt  

def wavenumberSelector(wanted_Z,all_Z):
    chosen_Z = []
    all_Z = all_Z.tolist()
    #print("original:",all_Z,"\n")
    sorted_all_Z = all_Z.copy()
    sorted_all_Z.sort()
    #print("sorted:",sorted_all_Z,"\n")
    int_all_Z = [int(float(sorted_all_Z[i])) for i in range(len(sorted_all_Z))]
    #print("sorted int:",int_all_Z,"\n")
    dec_all_Z = [round((float(sorted_all_Z[i])),1) for i in range(len(sorted_all_Z))]
    #print("sorted float:",dec_all_Z,"\n")
    for i, wantZ in enumerate(wanted_Z):
        #print("wanted Z:",wantZ) 
        isInt = isIntChecker(wantZ) 
        #print(isInt) 
        if isInt == 1:
            for j, Z in enumerate(int_all_Z):
                if wantZ == Z:
                    found_Z = sorted_all_Z[j] 
                    #print("found Z:",sorted_all_Z[j])
                    #print("sorted index:",j)
                    true_index = all_Z.index(found_Z) 
                    #print("true index:",true_index)
                    if true_index < len(all_Z):
                        if (abs(found_Z-wantZ) < abs(sorted_all_Z[j+1]-wantZ)):
                            chosen_Z.append(found_Z)
                            break
                    elif (abs(found_Z-wantZ) > abs(sorted_all_Z[j+1]-wantZ)):
                            continue 
                    else:
                        break 
        elif isInt == 0: 
            for j, Z in enumerate(dec_all_Z):
                if wantZ == Z:
                    found_Z = sorted_all_Z[j] 
                    #print("found Z:",sorted_all_Z[j])
                    #print("sorted index:",j)
                    true_index = all_Z.index(found_Z) 
                    #print("true index:",true_index)
                    if true_index < len(all_Z):
                        if (abs(found_Z-wantZ) < abs(sorted_all_Z[j+1]-wantZ)):
                            chosen_Z.append(found_Z)
                            break
                    elif (abs(found_Z-wantZ) > abs(sorted_all_Z[j+1]-wantZ)):
                            continue 
                    else:
                        break 
    return chosen_Z 

def PowerFinder(data_vec,spec_k):
    listk = data_vec[:,1] 
    power = 0
    for i, k in enumerate(listk):
        if k == spec_k:
            power = data_vec[i:i+1,2] 
            power = power[0] #convert fom list to float
    return power  
    
#----------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------
#program 
#----------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------


#prelims
emissivity_file = '/expanse/lustre/scratch/andrewcaruso/temp_project/modified_powers/21cm_spectra_subVolumes/' #directory of the sliced 21cm power spectra data to use. CHANGE ME
model_list = ["democratic/","fiducial/","oligarchic/"]
emissivity_models = ["","",""]
emissivity_models = [emissivity_file+model_list[i] for i in range(3)] 
output_path = "/expanse/lustre/scratch/andrewcaruso/temp_project/modified_powers/" #CHANGE ME
output_mainDir = "plots_convergenceTest_PvsZ" 
output_dir = output_path+output_mainDir 
dirChecker(output_dir)

#get all redshifts for each model (in case there is a discrepancy)
z_arr0,z_str0 = getRS(emissivity_models[0])
z_arr1,z_str1 = getRS(emissivity_models[1])
z_arr2,z_str2 = getRS(emissivity_models[2])
z_arr = np.array([z_arr0,z_arr1,z_arr2])
z_arr_str = np.array([z_str0,z_str1,z_str2])
current_dir = os.getcwd() #get current directory 
wanted_k = [0.1,0.2,0.3]  
#need to pick out wavenumbers: 0.1, 0.2, and 0.3 

#loop over models 
for m, model_dir in enumerate(emissivity_models):
    #STEP 1: get all the redshifts 
    z_all = z_arr[m] #all redshift values  
    original_z_all = z_all.copy() #make an original copy of the unordered redshifts
    z_string = z_arr_str[m] #string redshift values  
    chopped_dir = model_dir.split('/')
    model_name = chopped_dir[-2:-1]
    model_dir2 = '/'.join(chopped_dir[:-3])
    output_dir1 = output_dir+"/"+model_name[0] 
    dirChecker(output_dir1) #check output main directory 
    #STEP 2: get the cuts from first redshift directory  
    z_strb = z_string[0]
    #print("z_strb:",z_strb)
    start = time.time() #begin timer 
    cut_dirb = model_dir+"z="+z_strb
    cuts, cuts_string = getCuts(cut_dirb,z_strb) 
    #STEP 2.5: get and remove max cut
    #max cut will be the same regardless of model 
    maxCut, cuts_string = getDefaultCut(cuts) 
    #STEP 3: get the allow k values from selected k values 
    #loop over cuts 
    for i, cut_str in enumerate(cuts_string):
        model_dir3b = cut_dirb+"/cut_"+cut_str  
        #print("input dir:",model_dir3b)
        spectra_listb = listGrab(model_dir3b) 
        #get the all k values from first file for each cut
        first_file = spectra_listb[0]
        first_fileDir = model_dir3b+"/"+first_file
        first_data = loadIn(first_fileDir)
        k_vec = np.transpose(np.delete(np.delete(first_data,0,1),1,1)) #obj,number,axis (1 for column 0 for row)
        #k_list = list(k_vec[0]) #extract single list and convert to list
        chosen_k = wavenumberSelector(wanted_k,k_vec[0])
        print("chosen wavenumbers:",chosen_k) 
        output_dir2 = output_dir1+"/cut_"+cut_str
        dirChecker(output_dir2) #check output 1st sub directory 
        interval_list = [] #create time-interval list for the wavenumber loop 
        #STEP 4: get power for each z fixing k, cut, and model
        #loop over chosen k
        for j, k in enumerate(chosen_k):
            start = time.time() #begin timer 
            #create blank canvas/plot and default arrays 
            fig,ax = plt.subplots(1,1)
            P_max = [] #max box power, not from the cut
            P_avg = [] #average power from all files in the cut 
            P_stdp = [] #positive STD from avg power 
            P_stdn = [] #negative STD from avg power
            output_dir3 = output_dir2+"/k="+str(k)
            dirChecker(output_dir3) #check output 2nd sub directory 
            #print("output dir:",output_dir3)
            #loop over redshifts  
            for n, z_str in enumerate(z_string):
                cut_dir = model_dir+"z="+z_str
                model_dir3 = cut_dir+"/cut_"+cut_str  
                spectra_list = listGrab(model_dir3) 
                #print("input dir:",model_dir3)
                #get max power 
                max_dir = cut_dir+"/cut_"+maxCut
                max_file = listGrab(max_dir)[0] #returns a list therefore get first and only element 
                max_dir = max_dir+"/"+max_file
                max_data = loadIn(max_dir)
                max_data.astype(float)
                powerMax = PowerFinder(max_data,k)
                P_max.append(powerMax)
                #get the average power from all data files in the cut at specific z (fixed k, model, and cut) 
                P_temp = [] #temporary list for power data (from the spectra files)
                #loop over all the spectra files
                for g, data_file in enumerate(spectra_list):
                    data_dir = model_dir3+"/"+data_file
                    data = loadIn(data_dir) 
                    powerData = PowerFinder(data,k)
                    P_temp.append(powerData)
                #take the average and STD of all temporary powers for the specific z
                avgPower = sum(P_temp)/len(P_temp)
                stdPower = np.std(P_temp) 
                #update the corresponding list 
                P_avg.append(avgPower)
                P_stdp.append(avgPower+stdPower) #y+error
                P_stdn.append(avgPower-stdPower) #y-error 
            #order the data in ascending order
            z_all = original_z_all.copy() #get unordered redshifts each time
            zP_sorted = sorted(zip(z_all,P_max,P_avg,P_stdp,P_stdn),key=itemgetter(0))
            z_all,P_max,P_avg,P_stdp,P_stdn = zip(*zP_sorted)
            #in wavenumber loop
            end = time.time() #end time  
            interval = end-start 
            interval_list.append(interval)
            #plot the data for the cut
            ax.set_title("21cm Power vs Redshift at k="+str(k))
            ax.plot(z_all,P_max,color='red',label=maxCut+" Mpc/h")
            ax.plot(z_all,P_avg,color='green',linestyle='dashed',label="Average for "+cut_str+" Mpc/h cut")
            ax.fill_between(z_all,P_stdn,P_stdp,alpha=0.5,color='green')
            ax.set_xlabel("Redshift") 
            ax.set_ylabel("Power mk^2") 
            #ax.set_xlim([0.1,0.5]) 
            ax.set_xscale("log")
            ax.set_yscale("log")
            ax.legend(loc='best')
            #output the plot in the output directory
            os.chdir(output_dir3)
            fig.show()
            outfile = "plot_cut"+cut_str+"_k="+str(k)+".pdf"
            fig.savefig(outfile,format='pdf')
            print("output plot: ",outfile)
            os.chdir(current_dir)
            #end wavenumber loop
    #in model loop
    avg_interval = sum(interval_list)/len(interval_list)
    print("Average wavenumber computation time:",avg_interval," seconds")
    #end of model for loop
print("\n \n ")
