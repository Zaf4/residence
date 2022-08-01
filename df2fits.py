import pandas as pd
import glob
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.optimize import curve_fit
import warnings
from zk import zk

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


def fitNshow(fname,savefits=True,savefigs=False,showfits=True,keyword=60):
    df = pd.read_csv(fname)
    df = zk.df_minimize(df)
    keyword = str(keyword)

    fname=fname.replace('_df.csv','')
    ###fonkisyon
    t = df['Time']
    values = df[df.columns[1:]]
    fits = pd.DataFrame()

    fits['Time'] = t

    eqs = [powerlaw,exp_decay,double_exp_decay]
    eqnames = ['PL','ED','DED']

    for v in values:
        val = values[v]

        for eq,name in zip(eqs,eqnames):
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

            print(name,v,fname,y[:5],x[:5])

            popt,_ = curve_fit(eq, x, y, maxfev=200000)
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
    
    return

if __name__ == '__main__':
    
    fnames = glob.glob('*df.csv')
    
    for fname in fnames:
        fitNshow(fname,savefits=False,savefigs=True,showfits=False,keyword=60)
        fitNshow(fname,savefits=False,savefigs=True,showfits=True,keyword=60)
