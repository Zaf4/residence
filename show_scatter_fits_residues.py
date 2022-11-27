import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.optimize import curve_fit
import warnings
import os
from time import perf_counter as perf

warnings.filterwarnings('ignore')

def deleteNaN(y):
    """
    delete NaN parts of the input array
    and returns time and respective values.

    """
    
    t = np.arange(len(y))+1
    val = np.array(y)
    
    t = t[~np.isnan(val)]
    val = val[~np.isnan(val)]
    

    return t,val


def double_exp(x,a,b,c,d):
    return a*np.exp(-x*b) + (d)*np.exp(-x*c)

def powerlaw(x,a,b):
	return a*x**(-b)


def exp_decay(x,a,b):
    return a*np.exp(-x*b)

def tri_exp(x,a,b,c,d,e,f):
    return a*np.exp(-x*b) + c*np.exp(-x*d)+e*np.exp(-x*f)

def quad_exp(x,a,b,c,d,e,f,g,h):
    return a*np.exp(-x*b) + c*np.exp(-x*d)+e*np.exp(-x*f)+g*np.exp(-x*h)

def value_fit(val,eq):
    
    t_range = np.arange(len(val))+1
    
    residual_t = np.zeros([len(val),2])
    
    t,val = deleteNaN(val)
    
    popt, pcov= curve_fit(eq, t, val, maxfev=2000000)
    residuals = (val- eq(t, *popt))/np.max(val)*np.max(t_range)#norm. to max value
    res_norm = residuals/len(val)#*len(t_range)#norm. to size
    ss_res_norm = np.sum(res_norm**2)
    # ss_res_norm = ss_res/len(t)
    
    y_fit = eq(t_range, *popt)#full time length
    y_fit[y_fit<1] = np.nan#too small values to be removed
    y_fit[y_fit>np.max(val)*1.2] = np.nan#too big values removed
    
    return y_fit,ss_res_norm

def arr_minimize(arr,method='median'):
    """
    Minimized 1d array by removing repeats,
    according to the given method. 
    ---
    methods: 'median' or 'average'
    """

    search = np.unique(arr) #arr of unique elements
    search = search[search>0] #remove nans
    
    arr1 = arr.copy()
    
    for s in search:
        positions, = np.where(arr==s)
        if method == 'median':
            mid = int(np.median(positions))
    
        elif method == 'average': 
            mid = int(np.average(positions))
        
        arr1[positions] = np.nan
        arr1[mid] = s #mid value is kept
        
    return arr1

def df_minimize(df):
    """
    minimizes dataframe values by removing repeating values.
    """
    for i in range(len(df.columns)):
        df.iloc[:,i] = arr_minimize(df.iloc[:,i]) #values minimized and returned

    return df

def hex2rgb(hex)->tuple:
  return tuple(int(hex[i:i+2], 16) for i in (1, 3, 5))

def hexlist2rgb(hex_colors)->list:
    return [hex2rgb(hex_color) for hex_color in hex_colors ]


def generate_palette(colors:list,num_colors:int)->np.ndarray:
    """
    

    Parameters
    ----------
    colors : list
        List of colors either by size 2 or 3.
    num_colors : int
        Number of the colors to proce.

    Returns
    -------
    palette
        Color palette of given colors with given number of colors.

    """
    
    rgb_colors = np.array(hexlist2rgb(colors))
    r,g,b = rgb_colors[:,0],rgb_colors[:,1],rgb_colors[:,2]
    #initing palette array --RGBA
    palette = np.zeros([num_colors,3])
    if len(colors)==2:
        #reds
        palette[:,0] = np.linspace(r[0],r[-1],num_colors)

        #greens
        palette[:,1] = np.linspace(g[0],g[-1],num_colors)

        #blues
        palette[:,2] = np.linspace(b[0],b[-1],num_colors)

    
    elif len(colors)==3:
        mid = int(num_colors/2)#middle of the array
        #reds
        palette[:mid,0] = np.linspace(r[0],r[1],mid)
        palette[mid:,0] = np.linspace(r[1],r[2],num_colors-mid)
        #greens
        palette[:mid,1] = np.linspace(g[0],g[1],mid)
        palette[mid:,1] = np.linspace(g[1],g[2],num_colors-mid)
        #blues
        palette[:mid,2] = np.linspace(b[0],b[1],mid)
        palette[mid:,2] = np.linspace(b[1],b[2],num_colors-mid)
        
    return palette/256

def graphit(df:pd.DataFrame,fits:pd.DataFrame,ind:int=0)->plt.Axes:
    """
    

    Parameters
    ----------
    df : pd.DataFrame
        Main data.
    fits : pd.DataFrame
        Equation fit of the main data.
    ind : int, optional
        index to start. The default is 0.

    Returns
    -------
    ax : plot
        seaborn scatter graph.

    """
    
    df.index = np.arange(len(df))+1
    fits.index = np.arange(len(df))+1
    ax = plt.figure()
    font = {'family': 'Arial',
            'weight': 'light',
            'size': 14,
            }
    
    sns.set(style='ticks',
        rc = {'figure.figsize':(5.6,4.2),
                  'font.weight':'light',
                  'font.family':'sans-serif',
                  'axes.spines.top':'False',
                  'axes.spines.right':'False',
                  'ytick.minor.size':'0',
                  'xtick.minor.size':'0',
                  'ytick.major.size':'10',
                  'xtick.major.size':'10'
                  
                  }
            )
    

    cols = len(df.columns)
    jump = int(cols/5)
    
    vlo = df.columns[ind]
    low = df.columns[ind+jump]
    mid = df.columns[ind+2*jump]
    hi =  df.columns[ind+3*jump]
    ns =  df.columns[ind+4*jump]
    
    deep = generate_palette(['#2456E5','#C70039','#581845'], 5)

    deeper = sns.color_palette(deep)
    paor = ['#2456E5','#E57E24','#454649','#C70039','#581845']
    

    # paor = ['#FF5733','#C70039','#581845', '#FF7955','#C7225B','#583967']
    
    #scatter plots
    sns.scatterplot(data=df[vlo],color=paor[0],s=50,edgecolor='black',linewidth=0.25)
    sns.scatterplot(data=df[low],color=paor[1],s=50,edgecolor='black',linewidth=0.25)
    sns.scatterplot(data=df[mid],color=paor[2],s=50,edgecolor='black',linewidth=0.25)
    sns.scatterplot(data=df[hi], color=paor[3],s=50,edgecolor='black',linewidth=0.25)
    sns.scatterplot(data=df[ns], color=paor[4],s=50,edgecolor='black',linewidth=0.25)
    #lineplot
    sns.lineplot(data=fits[vlo],color=paor[0],linestyle='dashed',alpha=0.8)
    sns.lineplot(data=fits[low],color=paor[1],linestyle='dashed',alpha=0.8)
    sns.lineplot(data=fits[mid],color=paor[2],linestyle='dashed',alpha=0.8)
    sns.lineplot(data=fits[hi], color=paor[3],linestyle='dashed',alpha=0.8)
    sns.lineplot(data=fits[ns], color=paor[4],linestyle='dashed',alpha=0.8)
    
    
    plt.legend(['very low kT','low kT','mid kT','high kT','Nonspecific'])
    plt.yscale('log')
    plt.xscale('log')
    plt.xlim([0.7,3e3])
    plt.ylim([0.5,1e6])
    plt.xlabel('Duration (a.u.)',fontdict=font)
    plt.ylabel('Occurence',fontdict=font)
    

    return ax

def graphit_better(df:pd.DataFrame,fits:pd.DataFrame)->plt.Axes:
    """
    

    Parameters
    ----------
    df : pd.DataFrame
        Main data.
    fits : pd.DataFrame
        Equation fit of the main data.

    Returns
    -------
    ax : plot
        seaborn scatter graph.

    """
    
    df.index = np.arange(len(df))+1
    fits.index = np.arange(len(df))+1
    ax = plt.figure()
    font = {'family': 'Arial',
            'weight': 'light',
            'size': 14,
            }
    
    sns.set(style='ticks',
        rc = {'figure.figsize':(5.6,4.2),
                  'font.weight':'light',
                  'font.family':'sans-serif',
                  'axes.spines.top':'False',
                  'axes.spines.right':'False',
                  'ytick.minor.size':'0',
                  'xtick.minor.size':'0',
                  'ytick.major.size':'10',
                  'xtick.major.size':'10'
                  
                  }
            )
    
    #creating color palette
    # colors = ['#2456E5','#C70039','#581845']
    colors = ['#2456E5','#C70039']
    deep = generate_palette(colors, df.shape[1])
    light = generate_palette(colors, df.shape[1])*1.7
    light[light>1] = 1
    deeper = sns.color_palette(deep)
    lighter = sns.color_palette(light)    
    
    for i,col in enumerate(df):
        #scatterplot
        sns.scatterplot(data=df[col],color=deeper[i],s=50,
                        edgecolor='black',linewidth=0.25)
    for i,col in enumerate(fits):
        #lineplot
        sns.lineplot(data=fits[col], color=lighter[i],
                     linestyle='dashed',alpha=0.95)

    name_list = [f'{x[:4]} kT' for x in df]
    
    plt.legend(name_list)
    plt.yscale('log')
    plt.xscale('log')
    plt.xlim([0.7,3e3])
    plt.ylim([0.5,1e6])
    plt.xlabel('Duration (a.u.)',fontdict=font)
    plt.ylabel('Occurence',fontdict=font)
    

    return ax



def graphit_df(df:pd.DataFrame,fits:pd.DataFrame,pl:str='viridis')->plt.Axes:
    
    df.index = np.arange(len(df))+1

    ax = plt.figure()
    font = {'family': 'Arial',
            'weight': 'light',
            'size': 14,
            }
    
    sns.set(style='ticks',
        rc = {'figure.figsize':(5.6,4.2),
                  'font.weight':'light',
                  'font.family':'sans-serif',
                  'axes.spines.top':'False',
                  'axes.spines.right':'False',
                  'ytick.minor.size':'0',
                  'xtick.minor.size':'0',
                  'ytick.major.size':'10',
                  'xtick.major.size':'10'
                  
                  }
            )
    
    #scatter plots
    sns.scatterplot(data=df,palette=pl,markers=(['o']*df.shape[1]),
                    s=50,edgecolor='black',linewidth=0.25)
    #lineplots
    sns.lineplot(data=fits,palette=pl,linestyle='dashed',alpha=0.8)


    # plt.legend(['low kT','mid kT','high kT'])
    plt.yscale('log')
    plt.xscale('log')
    plt.xlim([0.7,3e3])
    plt.ylim([0.5,1e6])
    plt.xlabel('Duration (a.u.)',fontdict=font)
    plt.ylabel('Occurence',fontdict=font)
    

    return ax

def residuals2boxplot(residues:pd.DataFrame,cpalette='flare'):
    #residues of equations
    ax = plt.figure()
    
    font = {'family': 'Arial',
            'weight': 'light',
            'size': 14,
            }
    sns.set(style='ticks',
        rc = {'figure.figsize':(5.6,4.2),
              'font.weight':'light',
              'font.family':'Arial',
              'axes.spines.top':'False',
              'axes.spines.right':'False',
              'ytick.minor.size':'0',
              'ytick.major.size':'10',
              'xtick.major.size':'10'
              
              }
        )
    
    #kwargs for the boxplot
    props = {
    'boxprops':{'edgecolor':'black'},
    'medianprops':{'color':'black'},
    'whiskerprops':{'color':'black'},
    'capprops':{'color':'black'}
    
    }
    
    #boxplot for the residual comparison
    sns.boxplot(data=residues,palette=cpalette,
                saturation=1,linewidth=0.7,showfliers=False,
                **props)

    plt.ylabel('Σ(Residuals Normalized)²',fontdict=font)
    plt.yscale('log')
    plt.ylim([1e-1,1e5])
    plt.savefig('./scatters/eqVSresidues.png', dpi=400,
                transparent=True,
                bbox_inches='tight')

    # residues.to_csv('residual.csv',index=False)

    
    return ax

def equation_fit(datafile:str)->None:
    """
    From the durations-occurences dataframe produces fits for given equations
    and saves the files.
    

    Parameters
    ----------
    datafile : str
        Path/to/datafile.

    Returns
    -------
    None
 
    """
    
    
    #durations file is read as the main raw data -->preprocessing
    durations = pd.read_csv(datafile,index_col=None)
    durations[durations == 0] = np.nan 
    durations = df_minimize(durations)
    durations.to_csv('./data/durations_minimized.csv',index=False)
    
    #names and functions are listed for the for loop
    eqnames = ['ED','DED','TED','QED','Powerlaw']
    equations = [exp_decay,double_exp,tri_exp,quad_exp,powerlaw]

    #initing residuals dataframe
    residues = pd.DataFrame()
    for name,equation in zip(eqnames,equations):

        fits = pd.DataFrame() #a fit data frame is opened for each equation
        ress = np.zeros([20]) #array to store residual values
        
        #for each column of durations file
        for i,d in enumerate(durations):
            #fitting to equation is done
            fits[d],resid_sum_sqr=value_fit(np.array(durations[d]),eq=equation)
            ress[i] = resid_sum_sqr #collecting residuals
    
        residues[name] = ress #residuals appended
        
        #fits saved to file with their shorthened names + .csv
        fits.to_csv(f'./data/{name}.csv',index=False)
        #residuals saved as .csv 
        residues.to_csv('./data/residuals.csv',index=False)        
                    
    return

def generate_fit_graph(datafile:str,fitfile:str,keyword:str)->None:
    """
    

    Parameters
    ----------
    datafile : str
        path/to/datafile.
    fitfile : str
        path/to/fitfile.

    Returns
    -------
    None

    """
    
    durations = pd.read_csv(datafile,index_col=None)
    fits = pd.read_csv(fitfile,index_col=None)
    fname = fitfile.split('/')[-1]
    fname = fname[:-4]
    
    cols = [x for x in durations if keyword in x] #columns wtih the keyword
    partial_data = durations[cols]
    partial_fits = fits[cols]
    ax = graphit_better(partial_data,partial_fits)
    co = partial_data.columns[1][-2:] #concentration
    plt.text(x=1,y=1,
             color='grey',
             s=f'Fit equation:{fname}\nConcentration:{co} µM')
    plt.savefig(f'./scatters/{fname}_{co}.png',dpi=400,
                transparent=False,bbox_inches='tight')
                    
    return 


def partial_dataframe(df:pd.DataFrame,keyword:str):
    
    sub_df_cols = [x for x in df if keyword in x]
    
    return df[sub_df_cols]

    
if __name__ == '__main__':

    # residues = equation_fit('./data/durations.csv')
    # residues = pd.read_csv('./data/residuals.csv',index_col=None)
    # ax = residuals2boxplot(residues,cpalette='husl')
    generate_fit_graph(datafile='./data/durations_minimized.csv',
                       fitfile='./data/QED.csv',
                       keyword='3.50')
    

            
            
            
            
            
            
            
            
            
            
            
            