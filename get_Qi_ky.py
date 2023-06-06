#!/usr/bin/env python

import numpy as np
import numpy.ma as ma
from glob import glob
import os

from netCDF4 import Dataset

def get_Qi_ky_nc(dirname,ncfile='stella.out.nc', tfit = 85):
    ncfile = dirname + "/" + ncfile
    if os.path.isfile(ncfile):
        with Dataset(ncfile,'r') as f:
            #f.replace('--', np.nan)
            #
            Qikxky = f.variables["qflx_kxky"][()][:,0,0,:,:,:] # t, zed, kx, ky
            jacob = f.variables['jacob'][()][:,0] # nzed, nalpha
            Qikxky = np.sum(Qikxky * jacob[None,:,None,None],axis=1)/np.sum(jacob) # t, kx, ky
            t = f.variables["t"][()]
            ky = f.variables["ky"][()]
            kx = f.variables["kx"][()]
            I = np.argsort(kx)
            kx = kx[I]
            Qikxky = Qikxky[:,I,:]
            Qiky = np.trapz(Qikxky,kx,axis = 1) #t, ky
        # remove invalid entries
        # at the end of vector
        # netcdf4 returns a masked array, not normal array. See
        # https://numpy.org/doc/stable/reference/maskedarray.generic.html
        #print(t.mask)
        if ma.is_masked(t):
            while t.mask[-1] == True:
                t = t[:-1]
                Qiky = Qi[:-1]

    istart = np.argmin(np.fabs(t - tfit))
    T = t[-1] - t[istart]
    Qiky = np.trapz(Qiky[istart:],t[istart:], axis=0)/T
    
    return Qiky, ky

if __name__ == "__main__":
    import sys
    import matplotlib.pyplot as plt
    if len(sys.argv) > 1:
        ds = sys.argv[1:]
    else:
        ds = ['.']

    fig, ax = plt.subplots()

    Ndirs = len(ds)
    cmap = plt.get_cmap("brg",Ndirs)
    legend = []
    
    for i,d in enumerate(ds):    
        Qi,ky = get_Qi_ky_nc(d)
        if not np.isnan(Qi).all():
            ax.plot(ky,Qi,marker='.',color=cmap(i))
            legend.append(d)
    ax.set_xlabel(r'$k_y \rho$')
    ax.set_ylabel(r'$Q_{\rm{i}}/[Q_{\rm{gB}}]$')
    ax.legend(legend)
    plt.show()
