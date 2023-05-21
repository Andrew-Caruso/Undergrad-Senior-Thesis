import numpy as np
import os
import time 
import matplotlib.pyplot as plt
from operator import itemgetter
from scipy import interpolate

#create power P vs redshift z plots for fixed wavenumber k  
#plot the old data (without cosmic density criterion) in same plot with new data (with cosmic density criterion)

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
    #interpolate the k and p data to get power at the exact k (instead of selecting closest value)
    listk = data_vec[:,1] 
    listp = data_vec[:,2] 
    f = interpolate.interp1d(listk,listp)
    power = f(spec_k)
    return power  


def z_sorter(original_zall,P_max,P_avg,P_stdp,P_stdn):
    #order the data in ascending order 
    z_all = original_zall.copy() 
    combined = list(zip(z_all,P_max,P_avg,P_stdp,P_stdn))
    combined.sort(key=lambda x: x[0]) 
    z_all,P_max,P_avg,P_stdp,P_stdn = zip(*combined)
    return list(z_all),list(P_max),list(P_avg),list(P_stdp),list(P_stdn)

    
#----------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------
#program 
#----------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------


#prelims
#new data denoted with 1
#old data denoted with 2
emissivity_file1 = '/expanse/lustre/scratch/andrewcaruso/temp_project/final_powers/21cm_spectra_subVolumes/' #new power spectra data with cosmic density constraint 
emissivity_file2 = '/expanse/lustre/scratch/andrewcaruso/temp_project/OLD_final_powers/21cm_spectra_subVolumes/' #old power spectra data WITHOUT cosmic density constraint
model_list = ["oligarchic/","extreme/"]
emissivity_models1 = ["",""]
emissivity_models2 = ["",""]
emissivity_models1 = [emissivity_file1+model_list[i] for i in range(2)] 
emissivity_models2 = [emissivity_file2+model_list[i] for i in range(2)] 
output_path = "/expanse/lustre/scratch/andrewcaruso/temp_project/OLD_final_powers/" #CHANGE ME
output_mainDir = "dualPlots_convergenceTest_PvsZ" 
output_dir = output_path+output_mainDir 
dirChecker(output_dir)

#get all redshifts for each model (in case there is a discrepancy)
#for new data 
z_arr01,z_str01 = getRS(emissivity_models1[0])
z_arr11,z_str11 = getRS(emissivity_models1[1])
#need to specify 'dtype=object' for following numpy arrays because the models differ in lengths (number of z)
z_arr1 = np.array([z_arr01,z_arr11],dtype=object)
z_arr_str1 = np.array([z_str01,z_str11],dtype=object)
#for old data 
z_arr02,z_str02 = getRS(emissivity_models2[0])
z_arr12,z_str12 = getRS(emissivity_models2[1])
#need to specify 'dtype=object' for following numpy arrays because the models differ in lengths (number of z)
z_arr2 = np.array([z_arr02,z_arr12],dtype=object)
z_arr_str2 = np.array([z_str02,z_str12],dtype=object)

current_dir = os.getcwd() #get current directory 
wanted_k = [0.1,0.2,0.3]  
#need to pick out wavenumbers: 0.1, 0.2, and 0.3 

#NOTE: need to ensure that the order of the old and new z,cuts,etc are all in the same order!!! 

#loop over models 
for m, model_dir1 in enumerate(emissivity_models1):
    model_dir2 = emissivity_models2[m]
    #STEP 1: get all the redshifts 
    z_all1 = z_arr1[m] #all redshift values  
    z_all2 = z_arr2[m] #all redshift values  
    original_z_all1 = z_all1.copy() #make an original copy of the unordered redshifts
    z_string1 = z_arr_str1[m] #string redshift values  
    z_string2 = z_arr_str2[m] #string redshift values  
    print("z_string1:",z_string1)
    print("z_string2:",z_string2)
    chopped_dir1 = model_dir1.split('/')
    chopped_dir2 = model_dir2.split('/')
    model_name1 = chopped_dir1[-2:-1]
    model_name2 = chopped_dir2[-2:-1]
    print("models: ",model_name1," ",model_name2)
    #model_dir1 = '/'.join(chopped_dir1[:-3])
    #model_dir2 = '/'.join(chopped_dir2[:-3])
    output_dir1 = output_dir+"/"+model_name1[0] 
    dirChecker(output_dir1) #check output main directory 
    #STEP 2: get the cuts from first redshift directory (any redshift will do, i.e. do not need to coincide)
    z_strb1 = z_string1[0]
    z_strb2 = z_string2[0]
    #print("z_strb:",z_strb)
    start = time.time() #begin timer 
    cut_dirb1 = model_dir1+"z="+z_strb1
    cut_dirb2 = model_dir2+"z="+z_strb2
    cuts1, cuts_string1 = getCuts(cut_dirb1,z_strb1) 
    cuts2, cuts_string2 = getCuts(cut_dirb2,z_strb2) 
    #sort both in ascending order
    cuts1.sort()
    cuts2.sort()
    cuts_string1.sort()
    cuts_string2.sort()
    #STEP 2.5: get and remove max cut
    #max cut will be the same regardless of model 
    maxCut1, cuts_string1 = getDefaultCut(cuts1) 
    maxCut2, cuts_string2 = getDefaultCut(cuts2) 
    #print("cuts_string:",cuts_string1," ",cuts_string2)
    #STEP 3: get the allow k values from selected k values 
    #loop over cuts 
    for i, cut_str1 in enumerate(cuts_string1):
        cut_str2 = cuts_string2[i] 
        model_dir3b1 = cut_dirb1+"/cut_"+cut_str1  
        model_dir3b2 = cut_dirb2+"/cut_"+cut_str2  
        #print("input dir:",model_dir3b)
        spectra_listb1 = listGrab(model_dir3b1) 
        spectra_listb2 = listGrab(model_dir3b2) 
        #get the all k values from first file for each cut
        first_file1 = spectra_listb1[0]
        first_file2 = spectra_listb2[0]
        first_fileDir1 = model_dir3b1+"/"+first_file1
        first_fileDir2 = model_dir3b2+"/"+first_file2
        first_data1 = loadIn(first_fileDir1)
        first_data2 = loadIn(first_fileDir2)
        k_vec1 = np.transpose(np.delete(np.delete(first_data1,0,1),1,1)) #obj,number,axis (1 for column 0 for row)
        k_vec2 = np.transpose(np.delete(np.delete(first_data2,0,1),1,1)) #obj,number,axis (1 for column 0 for row)
        #k_list = list(k_vec[0]) #extract single list and convert to list
        #chosen_k = wavenumberSelector(wanted_k,k_vec[0])
        chosen_k = wanted_k
        print("chosen wavenumbers:",chosen_k) 
        #print("cut_str: ",cut_str1," ",cut_str2)
        output_dir2 = output_dir1+"/cut_"+cut_str1
        dirChecker(output_dir2) #check output 1st sub directory 
        interval_list = [] #create time-interval list for the wavenumber loop 
        #STEP 4: get power for each z fixing k, cut, and model

        #loop over chosen k
        for j, k in enumerate(chosen_k):
            start = time.time() #begin timer 
            #create blank canvas/plot and default arrays 
            fig,ax = plt.subplots(1,1) #for P vs z plot
            fig2,ax2 = plt.subplots(1,1) #for corresponding error plot
            P_max1 = [] #max box power, not from the cut
            P_max2 = [] #max box power, not from the cut
            P_avg1 = [] #average power from all files in the cut 
            P_avg2 = [] #average power from all files in the cut 
            P_stdp1 = [] #positive STD from avg power 
            P_stdp2 = [] #positive STD from avg power 
            P_stdn1 = [] #negative STD from avg power
            P_stdn2 = [] #negative STD from avg power
            output_dir3 = output_dir2+"/k="+str(k)
            dirChecker(output_dir3) #check output 2nd sub directory 
            #print("output dir:",output_dir3)

            #loop over redshifts  
            for n, z_str1 in enumerate(z_string1):
                #need to get SAME redshift for both new and old data 
                b = z_string2.index(z_str1)
                z_str2 = z_string2[b] 
                #NOTE: both z_str1 and z_str_2 are used BUT only z_str1 order will be plotted!!!
                cut_dir1 = model_dir1+"z="+z_str1
                cut_dir2 = model_dir2+"z="+z_str2
                model_dir31 = cut_dir1+"/cut_"+cut_str1  
                model_dir32 = cut_dir2+"/cut_"+cut_str2  
                spectra_list1 = listGrab(model_dir31) 
                spectra_list2 = listGrab(model_dir32) 
                #print("input dir:",model_dir3)
                #get max power 
                max_dir1 = cut_dir1+"/cut_"+maxCut1
                max_dir2 = cut_dir2+"/cut_"+maxCut2
                max_file1 = listGrab(max_dir1)[0] #returns a list therefore get first and only element 
                max_file2 = listGrab(max_dir2)[0] #returns a list therefore get first and only element 
                max_dir1 = max_dir1+"/"+max_file1
                max_dir2 = max_dir2+"/"+max_file2
                max_data1 = loadIn(max_dir1)
                max_data2 = loadIn(max_dir2)
                max_data1.astype(float)
                max_data2.astype(float)
                powerMax1 = PowerFinder(max_data1,k)
                powerMax2 = PowerFinder(max_data2,k)
                P_max1.append(powerMax1)
                P_max2.append(powerMax2)
                #get the average power from all data files in the cut at specific z (fixed k, model, and cut) 
                P_temp1 = [] #temporary list for power data (from the spectra files)
                P_temp2 = [] #temporary list for power data (from the spectra files)
                #for new data
                #loop over all the spectra files
                for g, data_file1 in enumerate(spectra_list1):
                    data_dir1 = model_dir31+"/"+data_file1
                    data1 = loadIn(data_dir1) 
                    powerData1 = PowerFinder(data1,k)
                    P_temp1.append(powerData1)
                #for old data
                #loop over all the spectra files
                for g, data_file2 in enumerate(spectra_list2):
                    data_dir2 = model_dir32+"/"+data_file2
                    data2 = loadIn(data_dir2) 
                    powerData2 = PowerFinder(data2,k)
                    P_temp2.append(powerData2)
                #take the average and STD of all temporary powers for the specific z
                avgPower1 = sum(P_temp1)/len(P_temp1)
                avgPower2 = sum(P_temp2)/len(P_temp2)
                stdPower1 = np.std(P_temp1) 
                stdPower2 = np.std(P_temp2) 
                #update the corresponding list 
                P_avg1.append(avgPower1)
                P_avg2.append(avgPower2)
                P_stdp1.append(avgPower1+stdPower1) #y+error
                P_stdp2.append(avgPower2+stdPower2) #y+error
                P_stdn1.append(avgPower1-stdPower1) #y-error 
                P_stdn2.append(avgPower2-stdPower2) #y-error 
                #end redshift loop
            #order the data in ascending order
            #use original_z_all1 for both, since redshift loop is over the new data z_string1
            z_all1,P_max1,P_avg1,P_stdp1,P_stdn1 = z_sorter(original_z_all1,P_max1,P_avg1,P_stdp1,P_stdn1)
            z_all2,P_max2,P_avg2,P_stdp2,P_stdn2 = z_sorter(original_z_all1,P_max2,P_avg2,P_stdp2,P_stdn2)
            print("Matching z_all 1 vs 2:",z_all1==z_all2)
            z_all = z_all1

            #in wavenumber loop
            end = time.time() #end time  
            interval = end-start 
            interval_list.append(interval)
            #plot the data for the cut
            #for the P vs z plot
            ax.set_title("21cm Power vs Redshift at k="+str(k))
            ax.plot(z_all,P_max1,color='red',label="Constrained:"+maxCut1+" Mpc/h")
            ax.plot(z_all,P_max2,color='tab:red',label="Unconstrained:"+maxCut2+" Mpc/h")
            ax.plot(z_all,P_avg1,color='green',linestyle='dashed',label="Constrained: Average for "+cut_str1+" Mpc/h cut")
            ax.plot(z_all,P_avg2,color='blue',linestyle='dashed',label="Unconstrained: Average for "+cut_str2+" Mpc/h cut")
            ax.fill_between(z_all,P_stdn1,P_stdp1,alpha=0.5,color='green')
            ax.fill_between(z_all,P_stdn2,P_stdp2,alpha=0.5,color='blue')
            ax.set_xlabel("Redshift") 
            ax.set_ylabel("Power mk^2") 
            ax.set_xscale("log")
            ax.set_yscale("log")
            ax.set_ylim([10**(-1),None]) #want to see the peak  
            ax.legend(loc='best')
            #for the error plot
            #create the power ratios out of the total power data  
            P_max31 = [P_max1[i]/P_max1[i] for i in range(len(P_max1))]
            P_max32 = [P_max2[i]/P_max2[i] for i in range(len(P_max2))]
            P_avg31 = [P_avg1[i]/P_max1[i] for i in range(len(P_max1))] 
            P_avg32 = [P_avg2[i]/P_max2[i] for i in range(len(P_max2))] 
            P_stdn31 = [P_stdn1[i]/P_max1[i] for i in range(len(P_max1))] 
            P_stdn32 = [P_stdn2[i]/P_max2[i] for i in range(len(P_max2))] 
            P_stdp31 = [P_stdp1[i]/P_max1[i] for i in range(len(P_max1))] 
            P_stdp32 = [P_stdp2[i]/P_max2[i] for i in range(len(P_max2))] 
            ax2.set_title("Ratios: 21cm Power vs Redshift at k="+str(k))
            ax2.plot(z_all,P_max31,color='red',label="Constrained:"+maxCut1+" Mpc/h")
            ax2.plot(z_all,P_max32,color='tab:red',label="Unconstrained:"+maxCut2+" Mpc/h")
            ax2.plot(z_all,P_avg31,color='green',linestyle='dashed',label="Constrained: Average for "+cut_str1+" Mpc/h cut")
            ax2.plot(z_all,P_avg32,color='blue',linestyle='dashed',label="Unconstrained: Average for "+cut_str2+" Mpc/h cut")
            ax2.fill_between(z_all,P_stdn31,P_stdp31,alpha=0.5,color='green')
            ax2.fill_between(z_all,P_stdn32,P_stdp32,alpha=0.5,color='blue')
            ax2.set_xlabel("Redshift") 
            ax2.set_ylabel("Power mk^2") 
            ax2.set_ylim([0,2])
            #ax2.set_xscale("log")
            #ax2.set_yscale("log")
            #ax2.set_ylim([10**(-1),10**(1)]) 
            ax2.legend(loc='best')
            #output the plot in the output directory
            os.chdir(output_dir3) #change to intended output directory
            #output the P vs z plot
            fig.show()
            outfile = "plot_cut"+cut_str1+"_k="+str(k)+".pdf"
            fig.savefig(outfile,format='pdf')
            print("output plot: ",outfile)
            #output corresponding error plot
            fig2.show()
            outfile2 = "error_plot_cut"+cut_str1+"_k="+str(k)+".pdf"
            fig2.savefig(outfile2,format='pdf')
            print("output plot: ",outfile2)
            os.chdir(current_dir) #return to current directory
            #end wavenumber loop
    #in model loop
    avg_interval = sum(interval_list)/len(interval_list)
    print("Average wavenumber computation time:",avg_interval," seconds")
    #end of model for loop

print("\n \n ")
