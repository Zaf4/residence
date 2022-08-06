import pandas as pd
import glob
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.optimize import curve_fit
import warnings
from zk import zk
import re

warnings.filterwarnings('ignore')


def deleteNaN(x,y,xmin=0):
    """
    delete NaN parts of the arrays
    and returns them.

    """
    a = np.array(x,dtype ='float32')[xmin:]
    b = np.array(y,dtype ='float32')[xmin:]
    
    c = a[~np.isnan(b)]
    d = b[~np.isnan(b)]
    

    return c,d
    
def powerlaw(x,a,b):
	return a*x**(-b)


def exp_decay(x,a,b):
    return a*np.exp(-x*b)


def double_exp_decay(x,a,b,c,d):
    return a*np.exp(-x*b) + (a-d)*np.exp(-x*c)

def name2case(cases):
    #conc
    rows = len(cases)
    kt,conc = np.zeros([rows]),np.zeros([rows])
    for i,c in enumerate(cases):
        r1,r2 = re.findall(r'[\d\.]+',c)
        r1,r2 = float(r1),float(r2)
        
        for val in [r1,r2]:
            if val<10:
                kt[i] = val
            else:
                conc[i] = val
                
    return kt,conc

def type2describe(event):
    evendict = {'4X'    : 'Potentials:SP>NS ',
                'NSI'   : 'Potentials:SP=NS ',
                'sp'    : 'Binding Type:Specific Sites ',
                'ns'    : 'Binding Type:Nonspecific Sites ',
                'feet'  : 'Legs:Legs Total ',
                'fl'    : 'Binding Type:NS and SP (Either Site) ',
                'lfoot' : 'Left Leg ',
                'rfoot' : 'Right Leg ',
                'nap'   : 'NAP ',
                'napstate': '(Either Site) '}
    
    
    descr = ''
    for key in evendict:
        if key in event:
            descr += evendict[key] + '\n'     
    
    return descr
    
    

def fitNshow(fname,savefits=True,savefigs=False,showfits=True,keyword=60):
    df = pd.read_csv(fname) 
    df = zk.df_minimize(df)#minimzing array to the median
    keyword = str(keyword)
    fname=fname.replace('_df.csv','')
    ###fonkisyon
    t = df['Time']
    values = df[df.columns[1:]]
    rows = len(values.columns)
    fits = pd.DataFrame()
    fits['Time'] = t

    #array to store exponential rates
    rates = np.zeros([rows,3],dtype='float32')
    rate_names = [' ']*rows
    case_names = [' ']*rows

    eqs = [powerlaw,exp_decay,double_exp_decay]
    eqnames = ['PL','ED','DED']

    for i,v in enumerate(values):
        val = values[v]
        
        rate_names[i] = fname
        case_names[i] = v

        for j,(eq,name) in enumerate(zip(eqs,eqnames)):
            x,y = deleteNaN(t, val)
            if len(y)<10:

                y_nan = np.zeros([len(t)])
                y_nan[:] = np.nan
                fits[f'{name}({v})'] = y_nan
                continue

            if eq == powerlaw:
                xmax = np.max(x)
                xmin = int(np.sqrt(xmax)*1.2)
                if xmin<len(x)+10:
                    x,y = x[xmin:],y[xmin:]

                if len(y)<10:

                    y_nan = np.zeros([len(t)])
                    y_nan[:] = np.nan
                    fits[f'{name}({v})'] = y_nan
                    continue

            #print(name,v,fname,y[:5],x[:5])

            popt,_ = curve_fit(eq, x, y, maxfev=200000)
            rates[i,j] = popt[1]

            y_fit = eq(t, *popt)
            y_fit[y_fit<1] = np.nan 
            y_fit[y_fit>np.max(val)*1.2] = np.nan
            fits[f'{name}({v})'] = y_fit

    if savefits == True:
        fits.to_csv(f'{fname}_fits.csv',index=False,encoding = 'UTF_8')
    
    #fits figure
    if showfits == True:
        fname +='with_fits'
        cols = fits.columns
        x = fits['Time']
        for col in cols:
            if keyword in col:
                y = fits[col]
                plt.plot(x,y,'k-')
                


    #figure part
    if savefigs == True:
        
        sns.set(rc = {'figure.figsize':(8,12)})
        sns.set_style("white")

        col2draw = [c for c in values.columns if keyword in c] 
        markers = ['o']*(len(values[col2draw].columns))
        sns.scatterplot(data=values[col2draw], markers= markers,palette='inferno')
        
        plt.yscale('log')
        plt.xscale('log')
        plt.xlim((0.8,25000))
        plt.ylim((0.8,10**6))
        plt.xlabel('Time (A.U)')
        plt.ylabel('Occurence')
        plt.title(fname)
        plt.savefig(f'{fname}_fig_.png', dpi=300)
        plt.clf()
    
    return rates,np.array([rate_names,case_names])

if __name__ == '__main__':
    
    fnames = glob.glob('*df.csv')
    
    store = np.zeros([3])
    store_names = np.array([['hahahaha'],['sbhfjdsfdsj']])
    
    for fname in fnames:
       # fitNshow(fname,savefits=False,savefigs=True,showfits=False,keyword=60)
       # fitNshow(fname,savefits=False,savefigs=True,showfits=True,keyword=60)
        rates,the_names= fitNshow(fname,savefits=False,savefigs=False,showfits=False)
        store = np.vstack((store,rates))
        store_names = np.hstack((store_names,the_names))
    
    namex = store_names[:,1:].copy()
    ratex = store[1:].copy()
    ratex = ratex.T
    ratex[ratex==0] = np.nan
    
    kt,conc = name2case(namex[1])
    
    
    #assinging the types of the cases
    rows = len(namex[0])
    Utype = ['']*rows
    for i,n in enumerate(namex[0]):
        if 'NSI' in n:
            Utype[i] = 'ns'
        else:
            Utype[i] = 'sp'
    
    #Dataframing part
    rates_df = pd.DataFrame()
    rates_df['type'] = namex[0]
    rates_df['kt'] = kt
    rates_df['concentration'] = conc
    rates_df['PL'] = ratex[0]
    rates_df['ED'] = ratex[1]
    rates_df['DED'] = ratex[2]
    rates_df['utype'] = Utype
    

    
    #for each event type
    types = np.array(namex[0])
    types = np.unique(types)
    
    for t in types:
        
        td = rates_df[rates_df['type']==t]
        sp = td[td['utype']=='sp']
        ns = td[td['utype']=='ns']
        
        if '4X' in t:
            #draw sp by conc
            sns.lineplot(data=sp,x='concentration',y='PL')
            sns.lineplot(data=sp,x='concentration',y='ED')
            sns.lineplot(data=sp,x='concentration',y='DED')
            plt.title(type2describe(t)+'\n SP by Concentration')
            plt.savefig(t+'_c.png',dpi=300)
            plt.show()
    
                    
            #draw sp by kt
            sns.lineplot(data=sp,x='kt',y='PL')
            sns.lineplot(data=sp,x='kt',y='ED')
            sns.lineplot(data=sp,x='kt',y='DED')
            plt.title(type2describe(t)+'\n SP by kT')
            plt.savefig(t+'_kt.png',dpi=300)
            plt.show()

        else:
            #draw ns by kt single conc 60um
            sns.lineplot(data=ns,x='kt',y='PL')
            sns.lineplot(data=ns,x='kt',y='ED')
            sns.lineplot(data=ns,x='kt',y='DED')
            plt.title(type2describe(t) + '\n NS by kT')
            plt.savefig(t+'_kt.png',dpi=300)
            plt.show()
            
        
        
    







