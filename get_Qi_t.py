#!/usr/bin/env python

import numpy as np
import numpy.ma as ma
from glob import glob
import os

from netCDF4 import Dataset, default_fillvals

def get_Qi_t_stdout(dirname,stdout_file='out.*'):
    stdout_file = sorted(glob(dirname + "/" + stdout_file))[-1]
    
    with open(stdout_file,'r') as f:
        lines = f.readlines()

    t = []
    Qi = []
    for l in lines:
        if "gx: Step" in l:
            ls = l.rsplit(':',1)[-1].split('=')
            for i,entry in enumerate(ls):
                if 'Time' in entry:
                    t.append(float(ls[i+1].split(',',1)[0]))
                if 'Q_i' in entry:
                    Qi.append(float(ls[i+1].split(None,1)[0]))
    return(np.array(Qi),np.array(t))


def get_Qi_t_nc(dirname,ncfile='stella.out.nc'):
    ncfile = dirname + "/" + ncfile
    if os.path.isfile(ncfile):
        with Dataset(ncfile,'r') as f:
            #f.replace('--', np.nan)
            Qi = f["qflx_kxky"][()][:,0,0] #t, (species), (tube), zed, kx, ky
            Qi = np.sum(Qi,axis = (2,3)) # sum over kx and ky
            jacob = f.variables['jacob'][()][:,0] # nzed, nalpha
            Qi = np.sum(Qi * jacob,axis=1)/np.sum(jacob)
            t = f.variables["t"][()]

        # remove invalid entries
        # at the end of vector
        # netcdf4 returns a masked array, not normal array. See
        # https://numpy.org/doc/stable/reference/maskedarray.generic.html
        #print(t.mask)
        if ma.is_masked(t):
            while t.mask[-1] == True:
                t = t[:-1]
                Qi = Qi[:-1]
        ret = (np.array(Qi),np.array(t))
    else:
        print("File does not exist '" + ncfile + "'")
        ret = (np.nan,np.nan)
    return ret



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
        Qi,t = get_Qi_t_nc(d)
        if not np.isnan(Qi).all():
            ax.plot(t,Qi,marker='.',color=cmap(i))
            legend.append(d)
    ax.set_xlabel(r'$t/[a/v_{\rm{Ti}}]$')
    ax.set_ylabel(r'$Q_{\rm{i}}/[Q_{\rm{gB}}]$')
    ax.legend(legend)
    plt.show()
