import numpy as np
import glob
import os
import sys
import time
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.optimize import curve_fit
import warnings

warnings.filterwarnings('ignore')



def deleteNaN(y):
    """
    delete NaN parts of the arrays
    and returns them.

    """
    
    t = np.arange(len(y))+1
    val = np.array(y)
    
    t = t[~np.isnan(val)]
    val = val[~np.isnan(val)]
    

    return t,val

def double_exp_decay(x,a,b,c,d):
    return a*np.exp(-x*b) + (d)*np.exp(-x*c)

def powerlaw(x,a,b):
	return a*x**(-b)


def exp_decay(x,a,b):
    return a*np.exp(-x*b)

def tri_exp(x,a,b,c,d,e,f):
    return a*np.exp(-x*b) + c*np.exp(-x*d)+e*np.exp(-x*f)

def value_fit(val,eq):
    
    t_range = np.arange(len(val))+1
    
    t,val = deleteNaN(val)
    
    popt, pcov= curve_fit(eq, t, val, maxfev=200000) 
    y_fit = eq(t_range, *popt)#full time length
    y_fit[y_fit<1] = np.nan#too small values to be removed
    y_fit[y_fit>np.max(val)*1.2] = np.nan#too big values removed
    
    return y_fit

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


def df_sum_types(df:pd.DataFrame)->pd.DataFrame:
    #finding column names with types
    corenames = [x for x in df.columns if 'core' in x]
    surfnames = [x for x in df.columns if 'surf' in x]
    freenames = [x for x in df.columns if 'free' in x]

    #summing up types in another dataframe
    sums = pd.DataFrame()
    sums['core'] = df[corenames].sum(axis=1)
    sums['surf'] = df[surfnames].sum(axis=1)
    sums['free'] = df[freenames].sum(axis=1)

    return sums

def graphit(df:pd.DataFrame)->plt.Axes:
    
    ax = plt.figure()
    font = {'family': 'Arial',
            'weight': 'bold',
            'size': 14,
            }
    sns.set(rc = {'figure.figsize':(6.6,4.4),
                  'font.weight':'bold',
                  'font.family':'Arial',
                  'font.size':12}
            )
    
    sns.set_style("white")
    sns.scatterplot(data=df.core,color='#2456E5',s=50,edgecolor='black')
    sns.scatterplot(data=df.surf,color='#E57E24',s=50,edgecolor='black')
    sns.scatterplot(data=df.free,color='#454649',s=50,edgecolor='black')
    #lineplot
    sns.lineplot(data=df.corefit,color='#2456E5',linestyle='dashed',alpha=.6)
    sns.lineplot(data=df.surffit,color='#E57E24',linestyle='dashed',alpha=.6)
    sns.lineplot(data=df.freefit,color='#454649',linestyle='dashed',alpha=.6)
    
    plt.yscale('log')
    plt.xscale('log')
    types  = ['Core','Surf','Free']
    plt.legend(types)
    plt.xlabel('Duration (a.u.)',fontdict=font)
    plt.ylabel('Occurence',fontdict=font)
    
    return ax

    

csv_files = glob.glob('./csv/*.csv')#finding all csv files
slices60 = [csv for csv in csv_files if '60' in csv]#60um cases
slices40 = [csv for csv in csv_files if '40' in csv]#40um cases

font = {'family': 'Arial',
        'weight': 'bold',
        'size': 13,
        }
font2 = {'family': 'Calibri',
        'size': 13,
        'color':'#FD4610'
        }
for slice60 in slices40:
    df = pd.read_csv(slice60)
    
    sums = df_sum_types(df)
    sums.iloc[-1] = np.nan####if full just delete -- data outlier
    sums = df_minimize(sums)
    eq = tri_exp
    sums['corefit'] = value_fit(sums.core,eq)
    sums['surffit'] = value_fit(sums.surf,eq)
    sums['freefit'] = value_fit(sums.free,eq)
    sums.index+=1
    graph = graphit(sums)
    name = slice60[6:-10]
    fname = name +'.png'
    descr = f'U : {name[:4]}kT\nC : {name[5:7]}µM\nS : {name[-2:].upper()}\nF : TED'
    plt.text(1, 1, descr, fontdict=font2,style='italic')
    plt.savefig('ted'+fname,dpi=200)
    df.to_csv(name+'.csv',index=False)
    

# for i in range(len(slices60)):
#     df = slices60[i]
#     df = pd.read_csv(df)

#     df[df==0] = np.nan
#     df = df_minimize(df)

#     df.to_csv(f'clean{i}.csv',index=False)













    
