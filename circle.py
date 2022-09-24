import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import math



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
    xco = np.zeros(len(angles))
    yco = np.zeros(len(angles))
    zco = np.zeros(len(angles))
    
    #for every angle a point (or particle) is created
    for i,a in enumerate(angles):
        xco[i] = r*np.cos(a)*np.sin(rotate)
        yco[i] = r*np.sin(a)*np.sin(rotate)
        zco[i] = r*np.cos(rotate)
        
    coor = pd.DataFrame()#empty DataFrame
    #coordinates are assigned as columns
    coor['x'] = xco
    coor['y'] = yco
    coor['z'] = zco
    
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

# coor = makesphere(r=10,step=200)
# sphere = np.array(coor)
# blo = np.random.rand(2000,3)*10-5
# blob = pd.DataFrame(blo,columns= ['x','y','z'])
# blo_type = np.zeros(len(blo))
# closest = np.zeros(len(sphere))

# for i,s in enumerate(sphere):
#     dist_sqr = np.sum((blo-s)**2,axis=1)#xyz(all)-xyz[i]
#     # dist = np.sqrt(dist_sqr)
#     index_min = np.argmin(dist_sqr)
#     closest[i] = index_min

# unique_close = np.unique(closest)
# blo_type[unique_close.astype(int)] = 1
# blob['surface'] = blo_type


def marksurface(arr:np.ndarray)->pd.DataFrame:
    
    r_max = np.sqrt(np.max(np.sum(arr**2,axis=1)))+1 # max of r_sqr is sqrt'ed 
    sphere = np.array(makesphere(r_max,step=int(r_max*10)))
    print(f'{len(sphere):.1f}',end='\t')
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

n = 10
vs = np.zeros([2,200])
for j in range(1,201):
    sizeb = 20*j
    blo = np.random.rand(sizeb,3)*n*2-n
    print(f'-----{j}-----')

    blob = marksurface(blo)
    surf = blob[blob.surface==1]
    # print(f'surface particles:{len(surf)} surface/total:{len(surf)/len(blob)}')
    core = blob[blob.surface==0]

    vs[0,j-1] = sizeb
    vs[1,j-1] = len(surf)/len(blob)
    
sns.scatterplot(x=vs[0], y=vs[1])
    # fig = plt.figure(figsize=(7,7))
    # ax = plt.axes(projection='3d')
    # ax.scatter(core.x,core.y,core.z,
    #            marker='o',alpha=1,
    #            color='grey')
    
    # ax.scatter(surf.x,surf.y,surf.z,
    #            marker='o',alpha=1,
    #            color='green')