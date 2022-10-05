import numpy as np
import os
import sys
import time
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns
import glob


def restate2D(arr1:np.ndarray,arr2:np.ndarray,threshold:float=0.9):
    """ Creates a zeros array by the size of the beads.
    Finds distances between the polymer and the free beads. 
    If a bead is near a polymer, it is considered bound (1) in the state array.
    Returns the state array.
    """
    lenf = len(arr2) #molecule to compare
    state= np.zeros(lenf,dtype=np.float32) #bound or not

    for i in range(lenf): #compare every free bead
            a = np.subtract(arr1,arr2[i]) #distance in each (x,y,z) direction
            distance = np.linalg.norm(a,axis=1)#euclidian distance

            if np.any(np.less(distance,threshold)): # if less than threshold distance
                state[i] = 1  #marked as bound

    return state


def restate3D(arrPol:np.ndarray,arrMol:np.ndarray, work_title:str = "Computing"):
    """
    Uses restate2D for every time depth.
    """
    time_range = len(arrPol) 
    mol_num = len(arrMol[0])
    
    state_whole = np.zeros([time_range,mol_num],dtype=np.float32) #creates a 3D state array to contain state for each bead at each time
    time_earlier = time.perf_counter()
    for i in range(time_range):
        state_whole[i] = restate2D(arrPol[i], arrMol[i])
        
        #performs occasional timing to estimate ending time.
        if i%10 == 1 or i == time_range-1:
            end_time = time.perf_counter()
            time_diff = end_time-time_earlier
            progressbar(i+1, time_range, time_diff,prefix=work_title)
    
    print(f'{work_title} done in {time_diff/60:.1f} minutes')
    
    return state_whole.T #formatting for the distribution analysis

def progressbar(step, length, time_elapsed,prefix = "Computing"):

    per = step/time_elapsed #performance
    remaining = length-step	#remaining steps

    #max bar length
    max_bar = 40
    step_norm = int(step/length*max_bar)
    remaining_norm = max_bar-step_norm

    eta = remaining/per

    unit = 's'
    if eta>60 and eta<3600:
        eta = eta/60
        unit = 'min'
    elif eta>3600:
        eta = eta/3600
        unit = 'h'

    ending = '\r'
    if step == length:
        ending = '\n'

    print(f'{prefix}: [{u"█"*step_norm}{" "*remaining_norm}]{step}/{length} -- Step/s:{per:.1f} -- ETA:{eta:.1f}{unit}     ',flush=True,end=ending)

def state_all(fname:str='dump.npy',ignore:int=1001)->None:
    """
    Creates a '.npz'file (using the dump file)
    storing all the state arrays and saves it.

    Parameters
    ----------
    fname : str, optional
        dump file to be read. The default is 'dump.npy'.
    ignore : int, optional
        DESCRIPTION. The default is 1001.

    Returns
    -------
    None
        DESCRIPTION.

    """
    arr = np.load(fname)
    part = arr[0]
    #slicing part to identify the indexes of dna sites and naps
    # dna_ind = np.where(part[:,0]<3)[0]
    sp_sites = np.where(part[:,0]==2)[0]
    ns_sites = np.where(part[:,0]==1)[0] 
    naps = part[(part[:,0]>2)]
    nap_ind = np.where(naps[:,0]!=5)[0]+len(sp_sites)+len(ns_sites)
    
    arr = arr[ignore:,:,1:]
    
    naps = arr[:,nap_ind].copy()      
    sp_dna = arr[:,sp_sites].copy() # only sp sites
    ns_dna = arr[:,ns_sites].copy() # only ns sites
    arr = 0
    #comparisons for both legs       
    
    sp_state = restate3D(sp_dna, naps, work_title='SP  ') #sp and naps
    ns_state = restate3D(ns_dna, naps, work_title='NS  ') #ns and naps
    nap_state = np.add(ns_state,sp_state)
    nap_state[nap_state==2] = 1
    """
    #legs from whole
    lfoot = nap_state[0::2].copy() #leg1s
    rfoot = nap_state[1::2].copy() #leg2s
    #sp ---legs
    sp_lfoot = sp_state[0::2].copy() #leg1s
    sp_rfoot = sp_state[1::2].copy() #leg2s
    #ns ---legs
    ns_lfoot = ns_state[0::2].copy() #leg1s
    ns_rfoot = ns_state[1::2].copy() #leg2s
     
    #addition of for full view
    fl_nap = np.add(lfoot,rfoot) #partials (foots) are added together 1 means partial 2 means complete binding
    sp_nap = np.add(sp_lfoot,sp_rfoot) #only sp interactions, 0:free, 1: one foot only, 2: both feet
    ns_nap = np.add(ns_lfoot,ns_rfoot) #only ns interactions, 0:free, 1: one foot only, 2: both feet
    
    #addition of sp to ns
    fl_lfoot = np.add(sp_lfoot,ns_lfoot) #left foot,- 0:free,1:only sp or ns, 2:both
    fl_rfoot = np.add(sp_rfoot,ns_rfoot) #right foot, 0:free,1:only sp or ns, 2:both
    #addition of both feet to view number of  interactions
    fl_feet = np.add(fl_rfoot,fl_lfoot) #gives how many sites a nap interacts
    """
    np.savez('states',ns_state=ns_state,sp_state=sp_state,
                    # nap_state=nap_state,
                    # lfoot=lfoot,rfoot=rfoot,
                    # sp_lfoot=sp_lfoot,sp_rfoot=sp_rfoot,
                    # ns_lfoot=ns_lfoot,ns_rfoot=ns_rfoot,
                    # fl_nap=fl_nap,sp_nap=sp_nap,ns_nap=ns_nap,
                    # fl_lfoot=fl_lfoot,fl_rfoot=fl_rfoot,fl_feet=fl_feet
                    )
    
    return



def read_residence(resident:np.ndarray)->np.ndarray:
    
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
        duration = len(resident[0])
        
        for res in resident: #for each row or the binding molecule
            time = 0 #set occurence time to zero
            if np.sum(res) == 0:
                continue
            elif 0 not in res:
                time=duration
                res_time.append(time)
                continue
            for r in res: #binding status 0 or 1+
                if r:
                    time+=1
                elif r==0 and time>0: #when the binding status is not satisfied anymore, binding duration is saved and time reset to 0
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
    
def save_residence():

    state = np.load('states.npz')['fl_nap']
    dist = read_residence(state)
    dist = np.where(dist == 0, np.nan, dist) #zeros removed
    
    #determining the name of the file and saving location
    sepr = "/" #for mac and linux
    if sys.platform == 'win32':
        sepr = '\\'
    
    cloc = os.getcwd().split(sepr)
    arr_name = cloc[-3]+'_'+cloc[-2]+'_'+cloc[-1]+'.npy'
    save_l = cloc[:-2]
    save_loc = sepr.join(save_l)
    os.chdir(save_loc)
    np.save(arr_name, dist)
    print(f'{arr_name} saved to {save_loc}')

    return dist


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
                    x = 'x',y = 'y',
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

    col_names = ['type','x','y','z']
    ss = pd.DataFrame(frame,columns=col_names)
    n_dna = len(ss[ss.type<3])
    tf = ss[ss.type>2]
    tf_coor = np.array(tf)[:,1:]

    
    if clusterit:
        cluster = find_cluster(tf_coor,min_size=min_size)#returns cluster array
        if len(cluster) == 0:
            ss['cluster'] = np.zeros(len(ss))
        else:
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
        coor = np.array(cluster_sub_df[['x','y','z']])
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
    
    plt.savefig(f'{arr_name[:-4]}.png',dpi=400)
    
    return ax


def makecircle(r:float=1, rotate:float=0, step:int=100)->pd.DataFrame:
    """

    Parameters
    ----------
    r : float, optional
        Radius of the circle. The default is 1.
    rotate : float, optional
        Radian angle to rotate the circle. The default is 0.
    step : int, optional
        Number of point in the circle. The default is 100.

    Returns
    -------
    coor : pd.DataFrame
        DataFrame with XYZ coordinates of the circle.

    """
    
    incre = (2*np.pi/step) #incrementation for desired steps
    angles =  np.arange(-np.pi,np.pi,incre)#len(angles) == step
    #initializing the array
    x = np.zeros(len(angles))
    y = np.zeros(len(angles))
    z = np.zeros(len(angles))
    
    #for every angle a point (or particle) is created
    for i,a in enumerate(angles):
        x[i] = r*np.cos(a)*np.sin(rotate)
        y[i] = r*np.sin(a)*np.sin(rotate)
        z[i] = r*np.cos(rotate)
        
    coor = pd.DataFrame()#empty DataFrame
    #coordinates are assigned as columns
    coor['x'] = x
    coor['y'] = y
    coor['z'] = z
    
    return coor

def makesphere(r:float,step:int=100)->pd.DataFrame:
    """

    Parameters
    ----------
    r : float, optional
        Radius of the circles -> sphere. The default is 1.

    step : int, optional
        Number of circles in the sphere. The default is 100.

    Returns
    -------
    coor : pd.DataFrame
        DataFrame with XYZ coordinates of the sphere.

    """
    incre = (2*np.pi/step) #incrementation for desired steps
    coor = makecircle(r,rotate=np.pi,step=step)
    
    circle_size = len(coor)#reference size of circle
    angles =  np.arange(-np.pi,np.pi,incre)
    sphere_size = circle_size*len(angles)
    
    sphere = np.zeros([sphere_size,3])#sphere array allocated
    
    #for every angle a circle is created
    for i,a in enumerate(angles):
        values = makecircle(r,rotate=a,step=step)
        sphere[i*circle_size:(i+1)*circle_size] = np.array(values) 
    
    coor = pd.DataFrame(sphere,columns= ['x','y','z'])#DataFrame formed
    
    return coor


def marksurface(arr:np.ndarray)->pd.DataFrame:
    """
    Marks surface atoms of XYZ np.ndarray into a dataframe.

    Parameters
    ----------
    arr : np.ndarray
        XYZ coordinates of the particles.

    Returns
    -------
    blob : pd.DataFrame
        DataFrame with particles marked as 0 for core 1 for surface in as an
        additional column. XYZ coordinates are preserved.

    """
    
    r_max = np.ceil(np.sqrt(np.max(np.sum(arr**2,axis=1)))) # max of r_sqr is sqrt'ed
    size = len(arr)
    size_sqrt = np.sqrt(size)
    
    step = np.ceil(size_sqrt)*4
    if r_max>size_sqrt:
        step = np.ceil(np.sqrt(size))*5
        
    
    sphere = np.array(makesphere(r_max,step=step))
    # print(len(sphere))
    # print(f'{len(sphere):.1f}',end='\t')
    blob = pd.DataFrame(arr,columns= ['x','y','z'])
    blo_type = np.zeros(len(arr))#type is surface or not surface (1 or 0)
    closest = np.zeros(len(sphere))#to store closest partcl. for each sp. point 
    
    for i,s in enumerate(sphere):
        dist_sqr = np.sum((arr-s)**2,axis=1)#xyz(all)-xyz[i]
        # dist = np.sqrt(dist_sqr)
        index_min = np.argmin(dist_sqr)
        closest[i] = index_min
        
    unique_close = np.unique(closest)
    blo_type[unique_close.astype(int)] = 1
    blob['surface'] = blo_type
    
    return blob

def cluster_surface(frame:pd.DataFrame)-> pd.DataFrame:
    
    naps = frame[frame.type==5]
    surf_arr = np.zeros(len(naps))
    #naps['surface'] = np.zeros(len(naps))
    naps['id'] = np.arange(len(naps))
    maxC = naps.cluster.max()
    
    for i in range(int(maxC)):
        cl = naps[naps.cluster==i+1]
        ids = np.array(cl['id'])
        xyz = np.array(cl[['x','y','z']])
        
        avg_xyz = np.average(xyz,axis=0)
        xyz_norm = xyz-avg_xyz
        
        cl_surf = marksurface(xyz_norm)
        
        surf_arr[ids] = cl_surf.surface
        
    naps['surface'] = surf_arr

    return naps[naps.cluster>0]

def slice_clusters(fname:str='states.npz',
                   dump_name:str='dump.npy',
                   slice_num:int=5,
                   duration:int=100,
                   ignore:int=1001,
                   part:str='sp_state')->list:
    """
    Loads a state part of '.npz'file and slice it according to position of the
    particle with respect to the clusters acquired from dump file. 

    Parameters
    ----------
    fname : str, optional
        filename of the '.npz' file. The default is 'states.npz'.
    dump_name : str, optional
        filename of the dump file. The default is 'dump.npy'.
    slice_num : int, optional
        Number of time points to be randomly selected. The default is 5.
    duration : int, optional
        Time length of the slice. The default is 100.
    ignore : int, optional
        Steps to be ignored in the dump file. The default is 1001.
    part : str, optional
        Part or the state to be worked on. The default is 'sp_state'.

    Returns
    -------
    list
        List of numpy.ndarrays which are core,surface, and free slices 
        respectively.

    """
    #state_all is time consuming check if already exists first
    if not os.path.exists(fname):
        state_all(ignore=ignore)#creates states.npz

    state = np.load(fname)[part]#loading npz file and taking a part of it
     
    tmax=len(state[0])#just finding time interval
    
    time_points = np.sort(np.random.randint(1,[tmax-duration]*slice_num))

    stateL = state[0::2]
    stateR = state[1::2]
    state = stateR+stateL
    state[state==2]=1
    
    #initing list to store slices arrays
    core_slices = []
    surf_slices = []
    free_slices = []
    for t in time_points:
        #framing clusters
        frame = cluster_single_frame(np.load(dump_name)[t])
        #marking surface atoms
        clx = cluster_surface(frame)
        if len(clx) == 0:
            break
             
        #indexes of core,surface and free (non-cluster atoms)
        core_ind = np.array(clx[clx.surface==0].id)
        surf_ind = np.array(clx[clx.surface==1].id)
        clx_ind = np.concatenate((surf_ind,core_ind))
        free_ind = [i for i in range (len(state)) if i not in clx_ind]
        free_ind = np.array(free_ind)
    

        #slices of the time intervals for core,surface, and free atoms.
        core_slice = state[core_ind,t:t+duration].copy() 
        surf_slice = state[surf_ind,t:t+duration].copy()
        free_slice = state[free_ind,t:t+duration].copy()
        #appending slices to their respective lists
        core_slices.append(core_slice)
        surf_slices.append(surf_slice)
        free_slices.append(free_slice)
    return [core_slices,surf_slices,free_slices]
       

if __name__ == '__main__':
    pd.set_option('mode.chained_assignment', None) 

    parts = ['sp_state.npy','ns_state.npy']
    
    path = os.getcwd()
    sepr = '/'
    if 'win' in sys.platform:
        sepr = '\\'
        
    path_parts = path.split(sepr)
    slice_name = '_'.join(path_parts[-2:])
    save_loc = sepr.join(path_parts[:-2])+sepr

    for part in parts:
        
        type_slices = slice_clusters(fname='states.npz',
                                            ignore=1001,
                                            part=part)
        
        df = pd.DataFrame()
        for n,ts in zip(['core','surf','free'],type_slices):
            for j,s in enumerate(ts):
                name = f'{n}{j}'
                df[name] = read_residence(s)
        
        
        p_name = part.split('_')[0]
        save_name = f'{save_loc}{slice_name}_{p_name}_slice.csv'
        df.to_csv(save_name,index=False)
           
    # for i in ([0,5,10]):
    #     ax = plt.figure()
    #     matplotlib.use('agg')
    #     sns.set(rc = {'figure.figsize':(10,6)})
    #     sns.set_style(style = 'white')
    #     sns.violinplot(data=df[df.columns[i:i+5]],
    #                    palette='magma')
        
    #     plt.savefig(f'{i}.png')