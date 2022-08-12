import pandas as pd
import glob
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from zk import zk




fnames = glob.glob('*df.csv')

print(fnames)

np = pd.read_csv(fnames[0]) 
ns = pd.read_csv(fnames[1])
sp = pd.read_csv(fnames[2])

np = zk.df_minimize(np)
ns = zk.df_minimize(ns)
sp = zk.df_minimize(sp)


sns.scatterplot(data=np[np.columns[14:16]],palette='Greys',markers =['o']*2)
sns.scatterplot(data=ns[np.columns[14:16]],palette='Blues',markers =['o']*2)
sns.scatterplot(data=sp[np.columns[14:16]],palette='Reds', markers =['o']*2)
sns.set(rc = {'figure.figsize':(12,12)})
sns.set_style("white")
plt.yscale('log')
plt.xscale('log')
plt.title('GREYS: All sites,\n BLUES: ns sites ----- REDS: sp sites')
plt.savefig('A350.png',dpi =600)








