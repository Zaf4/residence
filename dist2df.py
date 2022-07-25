import numpy as np
import os
import pandas as pd
import glob
import sys

def read_residence(resident):
    
    """
    Finds the distribution of occurence with respect to time.
    Returns only the occurences.
    """
    
    def resider(val=1):
        """"
        Measures the time of staying bound.
        Appends them to a list and return the list.
        """
        res_time = []
        
        for res in resident: #for each row or the binding molecule
            time = 0 #set occurence time to zero
            for r in res: #binding status 0 or 1,2
                if r:
                    time+=1
                elif r==0: #when the binding status is not satisfied anymore, binding duration is saved and time reset to 0
                    res_time.append(time)
                    time = 0
        return np.array(res_time)   
             
    # unbound = resider(0) #unbound residences
    bound = resider(1) #bound residences
    
    val_max = len(resident[0]) #max time value to decide array length
    length = val_max
    serial = np.arange(1,val_max+1,dtype=int) #time values 
    distribution = np.zeros([length,2],dtype='float32')#empty array to store durations vs occurence rates
    distribution[:,0] = serial #first row is durations
    
    for i in range(length):#counting occurences
    
        distribution[i,1] = np.count_nonzero(bound == distribution[i,0]) #how many occurences for occurences observed for each duration

    
    return distribution[:,1].copy()

def find_npz():
    
    """
    finds .npz files in all dirs or subdirs. 
    Returns their locations.
    """
    file_list = []
    for dirpath, dirnames, filenames in os.walk("."):
        for filename in [f for f in filenames if f.endswith('npz')]:
            file_list.append(os.path.join(dirpath))     
        
    return file_list
        
def retrieve_sub_array(sub_arr='fl_nap'):

    """
    Retrieves given sub arrays of .npz files. 
    Saves them to grandparent directories.
    Returns the grandparent directories.
    """
        
    flies = find_npz()
    locations = []
    cwd = os.getcwd()
    for f in flies:
        
        os.chdir(f)
    
        state = np.load('states.npz')[sub_arr]
        #mol_num, time_num = len(state),len(state[0])
        dist = read_residence(state)
        dist = np.where(dist == 0, np.nan, dist) #zeros removed
        
        #determining the name of the file and saving location
        sepr = "/" #for mac and linux
        if sys.platform == 'win32':
            sepr = '\\'
        
        cloc = os.getcwd().split(sepr)
        
        arr_name = sub_arr +'_'+ cloc[-3]+'_'+cloc[-2]+'_'+cloc[-1]+'.npy'
        save_l = cloc[:-2]
        
        save_loc = sepr.join(save_l)
        
        locations.append(save_loc)
        os.chdir(save_loc)
        
        np.save(arr_name, dist)
        print(f'{arr_name} saved to {save_loc}')
        os.chdir(cwd)
    
    return locations

def unite_arrays(locations):
    
    """
    Finds 1d npy array files in the directory to unite them in a 2d array 
    first columns is time array.
    Array is transformed into pandas dataframe with column names 
    obtained from file names.
    Returns the dataframe
    """
    cwd = os.getcwd()
    for l in locations:
        os.chdir(l)
        
        arr_list = glob.glob('*.npy')
        
        first = np.load(arr_list[0])#to check overall structure of the arrays
        ncol= len(arr_list)+1 #additional row is 'time'
        nrow = len(first) #time num
        
        
        total = np.zeros([nrow,ncol]) #big 2d array to contain all the data
        total[:,0] = np.arange(1,nrow+1) #first cloumns is time
        
        #loading each .npy file ass array 
        for i,arr in enumerate(arr_list):
            total[:,i+1] = np.load(arr)#arrays assigned to their relative columns.
            os.remove(arr)
        
        titles = ['Time']
        cases = []
        
        for a in arr_list:
            parts = a[:-4].split('_') #.npy and work id removed
            conc = parts[-1]+"\u03BCM"
            case = parts[-2]+"kT"
            titles.append(case+conc)
            cases.append(case)
            cases = list(set(cases))
        
        
        df = pd.DataFrame(total,columns = titles)
        
        df_name = parts[0]+parts[1]+parts[-3]+'_df.csv'
        df.to_csv(df_name)
        print(f'Saved {df_name} to {l}.... Arrays deleted')

        
        os.chdir(cwd)
    
    return 

if __name__ == '__main__':
    
    
    locations = retrieve_sub_array('fl_nap')
    locations = list(set(locations))
    unite_arrays(locations)