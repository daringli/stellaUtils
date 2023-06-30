#!/usr/bin/env python

import numpy as np
from scipy.io import netcdf
import matplotlib.pyplot as plt
import os

from debug_ql_fluxes import  ql_flux_debug_from_dir


if __name__ == "__main__":
    import sys
    argv = sys.argv
    argc = len(argv)
    if argc > 1:
        dirs = argv[1:]
    else:
        dirs = ['.']


    fig, axes =plt.subplots(2,sharex=True)
    for d in dirs:
        ky, impurity_flux, energy_flux,Imomega_avg, phi2_vs_kxky, Ky2, zed, jacob = ql_flux_debug_from_dir(d)
        nzed = len(zed)
        delzed = np.diff(zed)
        dl_over_b = jacob[:,0] * delzed[0]
        print(np.sum(dl_over_b)) # should be roughly resolution independent

        Qql = energy_flux[1:,0]/phi2_vs_kxky[:,0]/ky**2
        Gammazql = impurity_flux[1:,0]/phi2_vs_kxky[:,0]/nzed/ky**2
        NGammaz = np.count_nonzero(~np.isnan(Gammazql))
        NQ = np.count_nonzero(~np.isnan(Qql))
        Qqlsum = np.sum(np.nan_to_num(Qql))/NQ
        Gammazqlsum = np.sum(np.nan_to_num(Gammazql))/NGammaz
        axes[0].plot(ky,Gammazql)
        axes[1].plot(ky, Qql)
        
        print("sum Q, sum Gammaz, sumGammaz/sumQ ratio")
        print(Qqlsum) # should be roughly resolution independent
        print(Gammazqlsum) # should be roughly resolution independent
        print(Gammazqlsum/Qqlsum) # should be roughly resolution independent
        


    axes[-1].set_xlabel(r"$k_y$")
    axes[0].set_ylabel(r"Quasilinear $\Gamma_z$")
    axes[1].set_ylabel(r"Quasilinear $Q$")
        
    axes[0].axhline(0,color='silver',linestyle='dotted')
    axes[1].axhline(0,color='silver',linestyle='dotted')
        
    axes[0].legend(dirs)    
    plt.show()
