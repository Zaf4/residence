import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

df = pd.read_csv('rates_sp.csv')
ns = df[df['type']=='ns']
sp = df[df['type']=='sp']

eq = ['powerlaw','exponential','d_exp_1']

#ns plot
plt.figure()
ax = sns.lineplot(data=ns,x='kt',y='powerlaw')
ax = sns.lineplot(data=ns,x='kt',y='exponential')
ax = sns.lineplot(data=ns,x='kt',y='d_exp_1')
plt.legend(eq)
plt.title('SP=NS')
plt.ylim((10**-2,10**1))
plt.yscale('log')
plt.savefig('SPeqNS(sp).png',dpi=300)

#sp plot
plt.figure()
ax = sns.lineplot(data=sp,x='kt',y='powerlaw')
ax = sns.lineplot(data=sp,x='kt',y='exponential')
ax = sns.lineplot(data=sp,x='kt',y='d_exp_1')
plt.legend(eq)
plt.title('SP>NS')
plt.ylim((10**-2,10**1))
plt.yscale('log')
plt.savefig('SPmtNS(sp).png',dpi=300)