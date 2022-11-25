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


def graphit(df:pd.DataFrame,fits:pd.DataFrame,ind:int=0)->plt.Axes:
    
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
    jump = int(cols/3)
    low,mid,hi = df.columns[ind],df.columns[ind+jump],df.columns[ind+2*jump]
    
    paor = ['#2456E5','#E57E24','#454649','#3567F6','#F68F35','#56575A']
    paor = ['#FF5733','#C70039','#581845', '#FF7955','#C7225B','#583967']
    
    #scatter plots
    sns.scatterplot(data=df[low],color=paor[0],s=50,edgecolor='black',linewidth=0.25)
    sns.scatterplot(data=df[mid],color=paor[1],s=50,edgecolor='black',linewidth=0.25)
    sns.scatterplot(data=df[hi], color=paor[2],s=50,edgecolor='black',linewidth=0.25)
    #lineplot
    sns.lineplot(data=fits[low],color=paor[3],linestyle='dashed',alpha=0.9)
    sns.lineplot(data=fits[mid],color=paor[4],linestyle='dashed',alpha=0.9)
    sns.lineplot(data=fits[hi], color=paor[5],linestyle='dashed',alpha=0.9)
    
    plt.legend(['low kT','mid kT','high kT'])
    plt.yscale('log')
    plt.xscale('log')
    plt.xlim([0.7,3e3])
    plt.ylim([0.5,1e6])
    plt.xlabel('Duration (a.u.)',fontdict=font)
    plt.ylabel('Occurence',fontdict=font)
    

    return ax

def eq_residuals(residues:pd.DataFrame,cpalette='flare'):
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

def equation_fit():
    durations = pd.read_csv('./data/durations.csv',index_col=None)
    durations[durations == 0] = np.nan
    durations = df_minimize(durations)
    
    eqnames = ['ED','DED','TED','QED','Powerlaw']
    equations = [exp_decay,double_exp,tri_exp,quad_exp,powerlaw]
    
    rexist = os.path.exists('./data/residuals.csv')
    if rexist:
        residues=pd.read_csv('./data/residuals.csv',index_col=None)
    else:
        residues = pd.DataFrame()
        
    for name,equation in zip(eqnames,equations):
        
        fexist = os.path.exists(f'./data/{name}.csv') 

        
        if fexist:
            fits = pd.read_csv(f'./data/{name}.csv',index_col=None)
        else:
            fits = pd.DataFrame()
        
            ress = np.zeros([12])
            for i,d in enumerate(durations):
                
                fits[d],resid_sum_sqr=value_fit(np.array(durations[d]),eq=equation)
                
                if not rexist:  
                    ress[i] = resid_sum_sqr
        
            if not rexist:
                residues[name] = ress
            
            fits.to_csv(f'./data/{name}.csv',index=False)
        
        #creating scatter and fits plots..
        for i,co in enumerate(['10','20','40','60']):
            plot = graphit(durations,fits,ind=i)
            plt.text(x=1,y=1,s=f'Fit equation:{name}\nConcentration:{co}',color='grey')
            plt.savefig(f'./scatters/{name}_{co}.png',dpi=400,
                        transparent=True,bbox_inches='tight')
        
        if not rexist:
            residues.to_csv('./data/residuals.csv',index=False)
                    
    return residues

if __name__ == '__main__':
    residues = equation_fit()
    # residues = pd.read_csv('./data/residuals.csv',index_col=None)
    ax = eq_residuals(residues,cpalette='husl')
    

            
            
            

            
            
    
            
            
            
            
            
            
            
            
            
            
            
            