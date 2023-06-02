import numpy as np
import os
import time 
import matplotlib.pyplot as plt
import math 

#create power P vs wavenumber k for fixed redshift z

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

def truncate(n,decimals):
    if isinstance(decimals,int) and decimals != 0:
        c = 10**decimals
    else:
        c = 1
    return (math.trunc(n*c)) / c

def redshiftSelector(wanted_Z,all_Z):
    chosen_Z = []
    #all_Z = all_Z.tolist()
    #print("original:",all_Z,"\n")
    sorted_all_Z = all_Z.copy()
    sorted_all_Z.sort()
    #print("sorted:",sorted_all_Z,"\n")
    int_all_Z = [int(float(sorted_all_Z[i])) for i in range(len(sorted_all_Z))]
    #print("sorted int:",int_all_Z,"\n")
    dec_all_Z = [((float(sorted_all_Z[i]))) for i in range(len(sorted_all_Z))]
    #print("sorted float:",dec_all_Z,"\n")
    for i, wantZ in enumerate(wanted_Z):
        #print("wanted Z:",wantZ)
        isInt = isIntChecker(wantZ)
        if isInt == 1:
            for j, Z in enumerate(int_all_Z):
                if wantZ == Z:
                    found_Z = sorted_all_Z[j]
                    #print("found Z:",sorted_all_Z[j])
                    #print("sorted index:",j)
                    true_index = all_Z.index(found_Z)
                    #print("true index:",true_index)
                    chosen_Z.append(found_Z)
                    break
        elif isInt == 0:
            for j, Z in enumerate(dec_all_Z):
                N = int(len(str(wantZ).replace(".",""))-1) 
                #print("number of digits:",N) 
                Z = truncate(Z,N)
                if wantZ == Z:
                    found_Z = sorted_all_Z[j]
                    #print("found Z:",sorted_all_Z[j])
                    #print("sorted index:",j)
                    true_index = all_Z.index(found_Z)
                    #print("true index:",true_index)
                    chosen_Z.append(found_Z)
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
emissivity_file = '/expanse/lustre/scratch/andrewcaruso/temp_project/final_powers/21cm_spectra_subVolumes/' #directory of the sliced 21cm power spectra data to use. CHANGE ME
model_list = ["oligarchic/","intermediate/","extreme/"]
emissivity_models = ["","",""]
emissivity_models = [emissivity_file+model_list[i] for i in range(3)] 
output_path = "/expanse/lustre/scratch/andrewcaruso/temp_project/final_powers/" #CHANGE ME
output_mainDir = "plots_convergenceTest_PvsK" 
output_dir = output_path+output_mainDir 
dirChecker(output_dir)

#get all redshifts for each model (in case there is a discrepancy)
z_arr0,z_str0 = getRS(emissivity_models[0])
z_arr1,z_str1 = getRS(emissivity_models[1])
z_arr2,z_str2 = getRS(emissivity_models[2])
#because ragged nested sequences in z_arr and z_arr_str, need to specify 'dtype=object' to prevent warning message
#due to only having the extreme model (intermediate and oligarchic models are empty 
z_arr = np.array([z_arr0,z_arr1,z_arr2],dtype=object)
#print("z_arr:",z_arr)
z_arr_str = np.array([z_str0,z_str1,z_str2],dtype=object)
current_dir = os.getcwd() #get current directory 
wanted_z = [6.6,6.9]  
print("wanted_z:",wanted_z)
#need to pick out redshifts: 6, 7, and 8.5
#print("z_arr_str:",z_arr_str) 
#print("z_arr_str[2]:",z_arr_str[2]) 
chosen_z = redshiftSelector(wanted_z,z_arr_str[2])
print("chosen_z:",chosen_z)

#loop over models 
for m, model_dir in enumerate(emissivity_models):
    chosen_z = redshiftSelector(wanted_z,z_arr_str[m])
    print("chosen redshifts:",chosen_z) 
    z = z_arr[m]
    #z_string = z_arr_str[m]
    z_string = chosen_z 
    chopped_dir = model_dir.split('/')
    model_name = chopped_dir[-2:-1]
    print("model:",model_name)
    model_dir2 = '/'.join(chopped_dir[:-3])
    output_dir1 = output_dir+"/"+model_name[0] 
    dirChecker(output_dir1) #check output main directory 
    interval_list = [] #create time-interval list for the redshift loop 
    #loop over redshifts per model
    for n, z_str in enumerate(z_string):
        start = time.time() #begin timer 
        cut_dir = model_dir+"z="+z_str
        cuts, cuts_string = getCuts(cut_dir,z_str) 
        #get and remove max cut
        maxCut, cuts_string = getDefaultCut(cuts) 
        output_dir2 = output_dir1+"/z="+z_str
        dirChecker(output_dir2) #check output 1st sub directory 
        #Get the power spectra from maxCut (only 1 data spectrum file)
        max_dir = cut_dir+"/cut_"+maxCut
        max_file = listGrab(max_dir)[0] #returns a list, therefore get first and only element 
        max_dir = max_dir+"/"+max_file
        max_data = loadIn(max_dir)
        max_data.astype(float)
        #print("max_data:\n",max_data)
        #print("max_data:\n")
        #[print(max_data[i]) for i in range(len(max_data))]
        #loop over cuts per redshift 
        for j, cut_str in enumerate(cuts_string):
            model_dir3 = cut_dir+"/cut_"+cut_str  
            output_dir3 = output_dir2+"/cut_"+cut_str
            dirChecker(output_dir3) #check output 2nd sub directory 
            print("input dir:",model_dir3)
            print("output dir:",output_dir3)
            spectra_list = listGrab(model_dir3) 
            #1st: create blank canvas/plot and default arrays 
            fig,ax = plt.subplots(1,1)
            fig2,ax2 = plt.subplots(1,1) #for corresponding error plot
            P_max = [] #max box power, not from the cut
            P_avg = [] #average power from all files in the cut 
            P_stdp = [] #positive STD from avg power 
            P_stdn = [] #negative STD from avg power
            #2nd: get the all k values from first file 
            first_file = spectra_list[0]
            first_fileDir = model_dir3+"/"+first_file
            first_data = loadIn(first_fileDir)
            k_vec = np.transpose(np.delete(np.delete(first_data,0,1),1,1)) #obj,number,axis (1 for column 0 for row)
            k_list = list(k_vec[0]) #extract single list and convert to list
            #Notice that k_vec and max_data may not be the same lenght!!!
            #3nd: loop over k value (which is same for all spectra files in the cut) 
            for i, k in enumerate(k_list):
                #find and store the max box's power
                powerMax = PowerFinder(max_data,k)
                P_max.append(powerMax)
                P_temp =  [] #temporary list for power data (from the spectra files)
                #loop over all the spectra files
                for g, data_file in enumerate(spectra_list):
                    data_dir = model_dir3+"/"+data_file
                    data = loadIn(data_dir) 
                    powerData = PowerFinder(data,k)
                    P_temp.append(powerData)
                #take the average and STD of all temporary powers for the specific k
                avgPower = sum(P_temp)/len(P_temp)
                stdPower = np.std(P_temp) 
                #update the corresponding list 
                P_avg.append(avgPower)
                P_stdp.append(avgPower+stdPower) #y+error
                P_stdn.append(avgPower-stdPower) #y-error 
            #create the power ratios out of the total power data  
            P_max2 = [P_max[i]/P_max[i] for i in range(len(P_max))]
            P_avg2 = [P_avg[i]/P_max[i] for i in range(len(P_max))] 
            P_stdn2 = [P_stdn[i]/P_max[i] for i in range(len(P_max))] 
            P_stdp2 = [P_stdp[i]/P_max[i] for i in range(len(P_max))] 
            #for the corresponding error plot 
            ax2.set_title("Ratios: 21cm Power vs Redshift at z="+z_str)
            ax2.semilogx(k_list,P_max2,color='red',label=maxCut+" Mpc/h")
            ax2.semilogx(k_list,P_avg2,color='green',linestyle='dashed',label="Average for "+cut_str+" Mpc/h cut")
            #ax2.plot(k_list,P_max2,color='red',label=maxCut+" Mpc/h")
            #ax2.plot(k_list,P_avg2,color='green',linestyle='dashed',label="Average for "+cut_str+" Mpc/h cut")
            ax2.fill_between(k_list,P_stdn2,P_stdp2,alpha=0.5,color='green')
            ax2.set_xlabel("Redshift") 
            ax2.set_xlabel("Wave number h/Mpc") 
            ax2.set_ylim([0.25,1.5])
            ax2.legend(loc='best')
            #for the P vs k plot  
            ax.set_title("21cm Power vs Wavenumber at z="+z_str)
            ax.plot(k_list,P_max,color='red',label=maxCut+" Mpc/h")
            ax.plot(k_list,P_avg,color='green',linestyle='dashed',label="Average for "+cut_str+" Mpc/h cut")
            ax.fill_between(k_list,P_stdn,P_stdp,alpha=0.5,color='green')
            ax.set_xlabel("Wave number h/Mpc") 
            ax.set_ylabel("Power mk^2") 
            ax.set_xlim([0.1,0.5]) #limit k values from 0.1 to 0.5
            ax.set_xscale("log")
            ax.set_yscale("log")
            ax.set_ylim([10**1,None])
            ax.legend(loc='best')
            #output the plots in the output directory
            os.chdir(output_dir3) #change to intended output directory 
            #output the P vs K plot 
            fig.show()
            outfile = "plot_cut"+cut_str+"_z="+z_str+".pdf"
            fig.savefig(outfile,format='pdf')
            print("output plot: ",outfile)
            #output the corresponding error plot
            fig2.show()
            outfile2 = "error_plot_cut"+cut_str+"_z="+z_str+".pdf"
            fig2.savefig(outfile2,format='pdf')
            print("output plot: ",outfile2)
            os.chdir(current_dir) #return to current directory 
        #in redshift loop
        end = time.time() #end time  
        interval = end-start 
        #print("Time to finish redshift:",interval)
        interval_list.append(interval)
    #in model loop
    if len(interval_list) != 0:
        avg_interval = sum(interval_list)/len(interval_list)
        print("Average redshift computation time:",avg_interval," seconds")
    #end of model for loop

print("\n \n ")
