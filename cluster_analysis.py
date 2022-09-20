import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns
import glob
import os
import imageio
from moviepy.editor import *
import time

#clustering algorithm
def find_cluster(atoms:np.ndarray,threshold:float=2.1,
                 min_size:int=5)->np.ndarray:
    """
    Takes an array(atoms) with the coordinates of atoms
    clusters atoms that are close to each other less
    than threshold.

    Returns cluster_arr 
    row indexes: unique clusters
    column indexes: atoms, marked 1 if associated.
    
    """
    
    num = len(atoms)
    groups = np.zeros((num,num),dtype='float32')
    
    #neighbors lists for each atom is found
    for i in range(num):
        distance = np.linalg.norm(np.subtract(atoms,atoms[i]),axis=1)#euclidian distance of distances in each (x,y,z) direction
        neighbors = distance<threshold
        neighbors[i] = False
        groups[i] = neighbors
    
    #neighbor lists are united repeats are deleted
    arr = groups.copy()
    jump = []#atoms that are assigned to a cluster to be jumped
    for j in range(num):
        if j in jump:
            continue
        index, = np.where(arr[j] > 0)
        
        for k in index:
            if np.sum(arr[j]) > 0 and k>j:
                arr[k] += arr[j]
                arr[j] = np.zeros(num,dtype='float32')
                ind, = np.where(arr[k] != 0)
                if np.sum(arr[ind]) == 0:
                    jump.append(k)
    
    
    c, = np.where(np.sum(arr,axis=1)>0)
    cluster_arr = arr[c].copy()
    
    cluster_arr[cluster_arr>1] = 1
    cluster_size = np.sum(cluster_arr,axis=1)/3 #3 atoms per protein, 
    cluster_arr = cluster_arr[cluster_size>min_size]
    return cluster_arr

def cluster_id(cluster_arr:np.ndarray,n:int=0)->np.ndarray:
    """
    from cluster_arr of clusters and atoms marked to be contained,
    returns a list (numpy array) that marks atoms with their respective
    cluster ids. 

    if there are not-to-be-clustered atoms with number (n) they are 
    initiated as 0 before the cluster atoms.
    """

    width = len(cluster_arr[0])+n #cluster atoms + other atoms
    cluster = np.zeros([width])

    for j,c in enumerate(cluster_arr):
        cluster_num = j+1 #giving ids to the cluster
        cluster_atoms =  np.where(c==1)[0] #indexes of atoms
        index_in_df = n+cluster_atoms  #indexes in the df (n, number of not-to-be-clustered atoms)
        cluster[index_in_df] = cluster_num

    return cluster

def imagify(frame:pd.DataFrame,palette:str='viridis')->plt.Axes:
    """

    Parameters
    ----------
    frame : pandas.DataFrame
          DataFrame from single timestep of the system containing types and
          coordinates of the particles  
    palette : str
        color palette to be used

    Returns
    -------
    ax : seaborn.scatterplot
        plot of cluster positions from x,y perspective.

    """
    #figure plotting
    ax = plt.figure()
    matplotlib.use('agg')
    sns.set(rc = {'figure.figsize':(14,7)})
    sns.set_style(style = 'white')
    sns.scatterplot(data=frame[frame.cluster>0], 
                    x = 'xco',y = 'yco',
                    hue='cluster',
                    palette=palette,
                    legend=None)
    plt.grid(False)
    
    plt.xlim([-95,95])
    plt.ylim([-40,40])

    return ax

def cluster_single_frame(frame:np.ndarray,min_size:int=12,
                         clusterit:bool=True)->pd.DataFrame:
    """
    Takes single time point as argument and returns a pandas.DataFrame with
    cluster IDs for each particle

    Parameters
    ----------
    frame : np.ndarray
        takes a single time point Frame.

    Returns
    -------
    ss
        DataFrame with cluster IDs for each particle.
   

    """

    col_names = ['type','xco','yco','zco']
    ss = pd.DataFrame(frame,columns=col_names)
    n_dna = len(ss[ss.type<3])
    tf = ss[ss.type>2]
    tf_coor = np.array(tf)[:,1:]

    
    if clusterit:
        cluster = find_cluster(tf_coor,min_size=min_size)#returns cluster array
        cluster_index = cluster_id(cluster,n=n_dna)
        ss['cluster'] = cluster_index
    
    return ss

def frameit_self(arr_name:str,skip:int=1)->None:
    """
    Create image for the given Frame based on the current cluster formations.
    Allows the tracking of nascent formations.
    Saves the images.

    Parameters
    ----------
    arr_name : str
        takes arr_name as input to load respective .npy file.

    Returns
    -------
    None.

    """
    
    frames = np.load(arr_name)[::skip]
    if not os.path.isdir('images_self'):
        os.mkdir('images_self')
    os.chdir('./images_self')
    for i,frame in enumerate(frames):
        ss = cluster_single_frame(frame)
        ax = imagify(ss)
        plt.savefig(f'{i}',dpi=100)
        print(f'{i} of {len(frames)}   ', end='\r')
    plt.close('all')
    os.chdir('..')
    return


def frameit_index(arr_name:str,index:int=-1,skip:int=1)->None:
    """
    Create image for the given Frame based on the latest cluster formations.
    Allows the tracking of movements that formed the clusters and latest frame.
    Saves the images.

    Parameters
    ----------
    arr_name : str
        takes arr_name as input to load respective .npy file.

    Returns
    -------
    None.

    """
    frames = np.load(arr_name)[::skip]
    cluster_index = cluster_single_frame(frames[index]).cluster

    
    if not os.path.isdir('images_index'):
        os.mkdir('images_index')
    os.chdir('./images_index')
    
    for i,frame in enumerate(frames):
        ss = cluster_single_frame(frame, clusterit = False)
        #adding cluster indexes of the last instance to the dataframe
        ss['cluster'] = cluster_index
        ax = imagify(ss)
        plt.savefig(f'{i}',dpi=100)
        print(f'{i} of {len(frames)}   ', end='\r')
    plt.close('all')
    os.chdir('..')

    return


def gifit(gif_name:str='clusters'):
    """
    Create a gif from the images in the present location save it with
    the given name.
    
    
    Parameters
    ----------
    gif_name : str, optional
             File name to be saved.

    Returns
    -------
    None.

    """
    names = [f'{i}.png' for i in range(1000)]
    with imageio.get_writer('{gif_name}.gif', mode='I') as movie:
        for name in names:
            image = imageio.imread(name)
            print(name)
            movie.append_data(image)
    return

def movieit(movie_name:str='clusters'):
    """
    Create a movie (.mp4) from the images in the present location save it with
    the given name.
    
    Parameters
    ----------
    movie_name : str, optional
        The default is 'clusters'.

    Returns
    -------
    None.

    """
    names_len = len(glob.glob('*.png'))
    names = [f'{i}.png' for i in range(names_len)]
    clip = ImageSequenceClip(names, fps = 24)
    clip.write_videofile(f'{movie_name}.mp4', fps = 24)
    
    return

def cluster_table(frame:np.ndarray):
    """
    Generates a cluster table (pd.DataFrame) that summarizes each cluster
    within the given frame.

    Parameters
    ----------
    frame : np.ndarray
        Single time point or frame of the 3d array (dump).

    Returns
    -------
    
    cluster_frame : pd.DataFrame
        Summary table of cluster
        
    eigM : np.ndarray
        eigenvalue vector for each cluster in a 2d array.

    """
    
    ss = cluster_single_frame(frame,min_size = 12)
    n_cluster = int(ss.cluster.max())
    liste = []
    
    eigM = np.zeros([n_cluster,3])

    for index in range(1,n_cluster+1):
        cluster_sub_df = ss[ss.cluster == index]
        coor = np.array(cluster_sub_df[['xco','yco','zco']])
        average = np.average(coor,axis=0)#to center the cluster
        coor_norm = np.subtract(coor,average)#cluster centerd
        std_xyz = np.std(coor_norm,axis=0)#
        rg_xyz = np.average(coor_norm**2,axis=0)
        rg_std = np.std(rg_xyz)
        rg_avg = np.average(rg_xyz)
        rg_min = np.min(rg_xyz)
        rg_max = np.max(rg_xyz)
        mm_avg = (rg_max+rg_min)/2
        
        x = average[0]
        y = average[1]
        z = average[2]
        
        #rg x-axis
        rgxx = np.sqrt(np.sum(abs((coor[:,0]-x)*(coor[:,0]-x))))
        rgxy = np.sqrt(np.sum(abs((coor[:,0]-x)*(coor[:,1]-y))))
        rgxz = np.sqrt(np.sum(abs((coor[:,0]-x)*(coor[:,2]-z))))
        #rg y-axis
        rgyy = np.sqrt(np.sum(abs((coor[:,1]-y)*(coor[:,1]-y))))
        rgyz = np.sqrt(np.sum(abs((coor[:,1]-y)*(coor[:,2]-z))))
        #rg z-axis
        rgzz = np.sqrt(np.sum(abs((coor[:,2]-z)*(coor[:,2]-z))))
        
        rgM = [[rgxx,rgxy,rgxz],
               [rgxy,rgyy,rgyz],
               [rgxz,rgyz,rgzz]]
        
        rgM = np.square(rgM)
        
        eigV,eigM1 = np.linalg.eig(rgM)
        eigV = np.linalg.eigvals(rgM)
        
        eigM[index-1] = eigV
        
        shape = 'globular'
        if np.max(eigV)>15*np.min(eigV):
            shape = 'filamentous'
        
        values = [
            index,
            shape,
            x,
            y,
            z,
            int(len(cluster_sub_df)/3),
            std_xyz[0],
            std_xyz[1],
            std_xyz[2],
            np.sum(rg_xyz),
            rg_xyz[0],
            rg_xyz[1],
            rg_xyz[2],
            rg_avg,
            ]
        
        liste.append(values)
        
    col_names = ['id',
                 'conformation',
                 'x',
                 'y',
                 'z',
                 'size', 
                 'std_x',
                 'std_y',
                 'std_z',
                 'Rg',
                 'Rgx',
                 'Rgy',
                 'Rgz',
                 'Rg_avg']
    cluster_frame = pd.DataFrame(data=liste,columns = col_names)
    return cluster_frame,eigM

def labeledGraph(arr_name:str,index:int=-1,savepng:bool=False)->plt.Axes:
    """
    Generate a labeled cluster graph.
    Label consist of clusterID and shape.

    Parameters
    ----------
    arr_name : str
        File to read.
    index : int, optional
        Time point. The default is -1.

    Returns
    -------
    ax : plt.Axes
        Labeled graph.

    """
    frame = np.load(arr_name)[index]
    ss = cluster_single_frame(frame)
    ax = imagify(ss)
    table,matrix = cluster_table(frame)
    
    for i in table.index:
        label =f'{i+1}{table.conformation[i][:1]}'
        plt.text(x=table.x[i]+3, y=table.y[i]+3,
                 s=label,backgroundcolor='red')
    
    plt.savefig(f'{arr_name[:-4]}.png', dpi=100)
    
    return ax

    
if __name__ == '__main__':
    # os.chdir('C:\\Users\\zaf4-PC\\Desktop\\temporary\\clusan')
    start = time.time()
    
    frameit_index('small_280_60.npy',skip=1)
    os.chdir('images_index')
    movieit()
    os.chdir('..')


    print(f'{time.time()-start:.2f} seconds')
    
    
    
    
    
    
    
