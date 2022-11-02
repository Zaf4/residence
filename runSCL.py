import state_clust_slice as scs
import os

locs = []
cwd = os.getcwd()
init = cwd+'/5X10X'
cases = ['/2.80','/3.00','/3.50']
concs = ['/10','/20','/40','/60']


for ca in cases:
    for co in concs:
        name =init+ca+co
        locs.append(name)

os.mkdir('./csv')
save_loc = cwd+'/csv'
for loc in locs:
    os.chdir(loc)
    scs.main_work(20,save_loc=save_loc)
    