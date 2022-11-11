import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.optimize import curve_fit
import warnings

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
            'weight': 'bold',
            'size': 14,
            }
    sns.set(rc = {'figure.figsize':(8,6),
                  'font.weight':'bold',
                  'font.family':'Arial',
                  'font.size':12}
            )
    
    sns.set_style("white")
    cols = len(df.columns)
    jump = int(cols/3)
    low,mid,hi = df.columns[ind],df.columns[ind+jump],df.columns[ind+2*jump]
    
    sns.scatterplot(data=df[low],color='#2456E5',s=50,edgecolor='black')
    sns.scatterplot(data=df[mid],color='#E57E24',s=50,edgecolor='black')
    sns.scatterplot(data=df[hi],color='#454649',s=50,edgecolor='black')
    #lineplot
    sns.lineplot(data=fits[low],color='#2456E5',linestyle='dashed',alpha=1)
    sns.lineplot(data=fits[mid],color='#E57E24',linestyle='dashed',alpha=1)
    sns.lineplot(data=fits[hi],color='#454649',linestyle='dashed',alpha=1)
    
    plt.legend([low,mid,hi])
    plt.yscale('log')
    plt.xscale('log')
    plt.xlim([0.5,25e3])
    plt.xlabel('Duration (a.u.)',fontdict=font)
    plt.ylabel('Occurence',fontdict=font)
    

    return ax

if __name__ == '__main__':
    durations = pd.read_csv('durations.csv')
    durations[durations == 0] = np.nan
    durations = df_minimize(durations)
    
    eqnames = ['exp_decay','double_exp','tri_exp','quad_exp','powerlaw']
    equations = [exp_decay,double_exp,tri_exp,quad_exp,powerlaw]
    
    residues = pd.DataFrame()
    for name,equation in zip(eqnames,equations):
        
        fits = pd.DataFrame()
        
        ress = np.zeros([12])
        for i,d in enumerate(durations):
            fits[d],resid_sum_sqr=value_fit(np.array(durations[d]),eq=equation)
            ress[i] = resid_sum_sqr
            
        residues[name] = ress
        
        for i,co in enumerate(['10','20','40','60']):
            plot = graphit(durations,fits,ind=i)
            plt.title(f'{name}_{co}')
            plt.savefig(f'./scatters/{name}_{co}.jpeg',dpi=400)
            
            
            
    #residues of equations
    plt.figure()

    sns.set(style='ticks',
        rc = {'figure.figsize':(9,6),
              'font.weight':'bold',
              'font.family':'Arial',
              'font.size':12
              }
        )

    sns.boxplot(data=residues,palette='viridis')
    plt.yscale('log')
    plt.xlabel('Equations')
    plt.ylabel('Sum of Squares of Residuals Normalized')
    plt.savefig('/scatters/eqVSresidues.png', dpi=200)
    
        
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            