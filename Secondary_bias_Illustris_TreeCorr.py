import treecorr
import numpy as np
import matplotlib.pyplot as plt




file_name='/home/ignacio/Documents/HAB/Catalogs/dfG_halos.fits'
cat = treecorr.Catalog(file_name, x_col='xHalo',y_col='yHalo',z_col='zHalo')
dd=treecorr.NNCorrelation(bin_size=0.1,min_sep=1e-6,max_sep=15)
dd.process(cat)

#Randoms
n = len(cat.x)
boxsize = max(cat.x)
rand_N = 3*n
seed = 42
np.random.seed(seed)
rand_X = np.random.uniform(0, boxsize, rand_N)
rand_Y = np.random.uniform(0, boxsize, rand_N)
rand_Z = np.random.uniform(0, boxsize, rand_N)
rand = treecorr.Catalog(x=rand_X,y=rand_Y,z=rand_Z)
rr = treecorr.NNCorrelation(min_sep=1e-6, max_sep=15, bin_size=0.1)
rr.process(rand)

#data-random cross correlation, DR estimators

dr=treecorr.NNCorrelation(min_sep=1e-6, max_sep=15, bin_size=0.1)
dr.process(cat,rand)
print(dr)

xi, varxi = dd.calculateXi(rr)
r = np.exp(dd.meanlogr)
sig = np.sqrt(varxi)


