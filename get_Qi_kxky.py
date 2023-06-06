#!/usr/bin/env python

import numpy as np
import numpy.ma as ma
from glob import glob
import os



from netCDF4 import Dataset

def get_Qi_kxky_nc(dirname,ncfile='stella.out.nc', tfit = 300):
    ncfile = dirname + "/" + ncfile
    if os.path.isfile(ncfile):
        with Dataset(ncfile,'r') as f:
            #f.replace('--', np.nan)
            # qflx_kxky(t, species, tube, zed, kx, ky)
            Qikxky = f.variables["qflx_kxky"][()][:,0,0,:,:,:] # t, zed, kx, ky
            jacob = f.variables['jacob'][()][:,0] # nzed, nalpha
            Qikxky = np.sum(Qikxky * jacob[None,:,None,None],axis=1)/np.sum(jacob) # t, kx, ky
            t = f.variables["t"][()]
            ky = f.variables["ky"][()]
            kx = f.variables["kx"][()]
            I = np.argsort(kx)
            kx = kx[I]
            Qikxky = Qikxky[:,I,:]
            
        # remove invalid entries
        # at the end of vector
        # netcdf4 returns a masked array, not normal array. See
        # https://numpy.org/doc/stable/reference/maskedarray.generic.html
        #print(t.mask)
        if ma.is_masked(t):
            while t.mask[-1] == True:
                t = t[:-1]
                Qikxky = Qikxky[:-1]

        istart = np.argmin(np.fabs(t - tfit))
        T = t[-1] - t[istart]
        Qikxky = np.trapz(Qikxky[istart:],t[istart:], axis=0)/T
    else:
        Qikxky = np.nan
        kx = np.nan
        ky = np.nan
    
    return Qikxky, kx, ky

if __name__ == "__main__":
    import sys
    import matplotlib as mpl
    from mpl_toolkits.axes_grid1 import make_axes_locatable
    import matplotlib.pyplot as plt
    if len(sys.argv) > 1:
        ds = sys.argv[1:]
    else:
        ds = ['.']

    
    Ndirs = len(ds)
    cmap = plt.get_cmap("brg",Ndirs)
    legend = []
    nrows = 1+ Ndirs//4
    ncols = min((4,Ndirs))
    
    fig, axes = plt.subplots(squeeze=False,nrows=nrows, ncols = ncols)
    axes = axes.flatten()
    vmax = 0.0
    vmin = 0.0

    Qis = []
    kxs = []
    kys = []
    n = 0
    for i,d in enumerate(ds):    
        Qi,kx,ky = get_Qi_kxky_nc(d)
        if not np.isnan(Qi).all():  
            Qis.append(Qi)
            kxs.append(kx)
            kys.append(ky)
            _vmax = np.nanmax(Qi)
            _vmin = np.nanmin(Qi)
            n = n + 1
            if _vmax > vmax:
                vmax = _vmax
            if _vmin < vmin:
                vmin = _vmin
                
    for i in range(n):
        Qi = Qis[i]
        kx = kxs[i]
        ky = kys[i]
        axes[i].pcolormesh(ky,kx,Qi,vmax=vmax,vmin=vmin)
        divider = make_axes_locatable(axes[i])
        cax = divider.append_axes('right', size='5%', pad=0.05)
        norm = mpl.colors.Normalize(vmin=vmin, vmax=vmax)
        cb1 = mpl.colorbar.ColorbarBase(cax, cmap='viridis', norm=norm,
                                        orientation='vertical')
        cb1.mappable.set_clim(vmin=vmin, vmax=vmax)
        cb1.draw_all()
                    
        axes[i].set_title(ds[i])
        axes[i].set_xlabel(r'$k_y \rho$')
        axes[i].set_ylabel(r'$k_x \rho$')

    plt.tight_layout(h_pad=-2,w_pad=-2.5)
    fig.suptitle(r"$Q_{kx,ky}$")
    plt.show()
