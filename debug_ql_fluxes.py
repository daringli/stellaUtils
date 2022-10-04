#!/usr/bin/env python

import numpy as np
from scipy.io import netcdf
import matplotlib.pyplot as plt
import os

def read_omega(filename = "stella.omega"):
    with open(filename,'r') as f:
        lines = f.readlines()

    t = []
    ky = []
    kx = []
    Reomega = []
    Imomega = []
    Reomega_avg = []
    Imomega_avg = []
    
    for l in lines:
        l = l.split("#")[0].strip()
        if len(l) > 0:
            ls = l.split()
            t.append(float(ls[0]))
            ky.append(float(ls[1]))
            kx.append(float(ls[2]))
            Reomega.append(float(ls[3]))
            Imomega.append(float(ls[4]))
            Reomega_avg.append(float(ls[5]))
            Imomega_avg.append(float(ls[6]))
    t = np.array(t)
    ky = np.array(ky)
    kx = np.array(kx)
    Reomega = np.array(Reomega)
    Imomega = np.array(Imomega)
    Reomega_avg = np.array(Reomega_avg)
    Imomega_avg = np.array(Imomega_avg)
    return (t,ky,kx,Reomega,Imomega,Reomega_avg,Imomega_avg)

def get_latest_omega(omega_tuple):
    t,ky,kx,Reomega,Imomega,Reomega_avg,Imomega_avg = omega_tuple
    last_time = t[-1]
    last_ky = ky[-1]
    i = np.where(t == last_time)[0][0]
    j = np.where(ky[i:] == last_ky)[0][0]
    kx = kx[i:][j:]
    Nkx = len(kx)
    if Nkx > 1:
        print("Warning! We have never ran this script with Nkx > 1. It will just plot for the first kx value...")
    ky = ky[i::Nkx]
    Nky = len(ky)
    Reomega = Reomega[i:].reshape((Nky,Nkx), order='C') 
    Imomega = Imomega[i:].reshape((Nky,Nkx), order='C')
    Reomega_avg = Reomega_avg[i:].reshape((Nky,Nkx), order='C') 
    Imomega_avg = Imomega_avg[i:].reshape((Nky,Nkx), order='C')
    
    return (last_time,ky,kx,Reomega,Imomega,Reomega_avg,Imomega_avg)


def ql_flux_debug_from_dir(dirname,average_from = None):
    tmp = os.getcwd()
    os.chdir(dirname)

    if os.path.exists("my_growthrate.npy"):
        Imomega_avg = np.load("my_growthrate.npy")
        ky = np.load("my_ky.npy")
        kx = np.array([0])
    else:
        omega_tuple = read_omega()
        t, ky, kx, Imomega, Reomega, Reomega_avg, Imomega_avg = get_latest_omega(omega_tuple)


    Ky,Kx = np.meshgrid(ky,kx,indexing='ij')

    Nky = len(ky)
    Nkx = len(kx)
    with netcdf.netcdf_file('stella.out.nc','r',mmap=False) as f:
        if average_from is not None:
            t = f.variables['t'][()]
            i = np.where(t>=average_from)[0][0]
            print(t[i])
            Nt = len(t[i:])
            phi2_vs_kxky  = np.transpose(np.sum(f.variables['phi2_vs_kxky'][()][i:],axis=0)/Nt)
        
        else:
            phi2_vs_kxky  = np.transpose(f.variables['phi2_vs_kxky'][()][-1])
        nz = f.variables['dens'][()][-1]
        zed = f.variables["zed"][()]
        # pflx_kxky(t, species, tube, zed, kx, ky)
        if average_from is not None:
            _pfluxes_vs_kxky = np.sum(f.variables['pflx_kxky'][()][i:],axis = 0)/Nt # time average
        else:
            _pfluxes_vs_kxky = f.variables['pflx_kxky'][()][-1] # last time step
        Nspecies = len(_pfluxes_vs_kxky)
        pfluxes_vs_kxky = np.zeros((Nspecies,Nkx,Nky))
        for i in range(Nspecies):
            pfluxes_vs_kxky[i] = np.sum(_pfluxes_vs_kxky[i],axis=(0,1))
        print(pfluxes_vs_kxky.shape)
        impurity_flux = pfluxes_vs_kxky[-1]/nz

        if average_from is not None:
            _qfluxes_vs_kxky = np.sum(f.variables['qflx_kxky'][()][i:],axis=0)/Nt
        else:
            _qfluxes_vs_kxky = f.variables['qflx_kxky'][()][-1]  # last time step
        qfluxes_vs_kxky = np.zeros((Nspecies,Nkx,Nky))
        for i in range(Nspecies):
            qfluxes_vs_kxky[i] = np.sum(_qfluxes_vs_kxky[i],axis=(0,1))
        print(qfluxes_vs_kxky.shape)
        print("-----------------")
        energy_flux = np.sum(qfluxes_vs_kxky,axis=0)

        energy_flux = np.transpose(energy_flux)
        impurity_flux = np.transpose(impurity_flux)
        jacob = f.variables["jacob"][()]

    
    normalization = Imomega_avg[1:]/(phi2_vs_kxky[1:] * Ky[1:]**2) # units??
    ql_impurity_flux = normalization[:,0] * impurity_flux[1:,0]
    ql_energy_flux = normalization[:,0] * energy_flux[1:,0]
    os.chdir(tmp)
    return ky[1:], impurity_flux, energy_flux, Imomega_avg[1:], phi2_vs_kxky[1:], Ky[1:]**2, zed, jacob


if __name__ == "__main__":
    import sys
    argv = sys.argv
    argc = len(argv)
    if argc > 1:
        dirs = argv[1:]
    else:
        dirs = ['.']


    fig, axes =plt.subplots(5,sharex=True)
    for d in dirs:
        ky, impurity_flux, energy_flux,Imomega_avg, phi2_vs_kxky, Ky2, zed, jacob = ql_flux_debug_from_dir(d)
        nzed = len(zed)
        delzed = np.diff(zed)
        dl_over_b = jacob[:,0] * delzed[0]
        print(np.sum(dl_over_b)) # should be roughly resolution independent

        Qql = energy_flux[1:,0]/phi2_vs_kxky[:,0]/ky**2
        Gammazql = impurity_flux[1:,0]/phi2_vs_kxky[:,0]/nzed/ky**2
        Qqlsum = np.sum(np.nan_to_num(Qql))
        Gammazqlsum = np.sum(np.nan_to_num(Gammazql))
        axes[0].plot(ky, impurity_flux[1:,0])
        axes[1].plot(ky, Qql)
        axes[2].plot(ky,phi2_vs_kxky[:,0])
        axes[3].plot(ky,Gammazql)
        axes[4].plot(ky,impurity_flux[1:,0]/phi2_vs_kxky[:,0]/nzed)

        print("sum Q, sum Gammaz, ratio")
        print(Qqlsum) # should be roughly resolution independent
        print(Gammazqlsum) # should be roughly resolution independent
        print(Gammazqlsum/Qqlsum) # should be roughly resolution independent
        


        axes[-1].set_xlabel(r"$\gamma$")
        axes[0].set_ylabel(r"raw $\Gamma_z$")
        axes[1].set_ylabel(r"raw $Q$ over $<|\varphi^2|>$ over $k_y$")
        axes[2].set_ylabel(r"$<|\varphi^2|>$")
        axes[3].set_ylabel(r"raw $\Gamma_z$ over $<|\varphi^2|>$ over $k_y$")
        axes[4].set_ylabel(r"raw $\Gamma_z$ over $<|\varphi^2|>$")
        
        axes[1].axhline(0,color='silver',linestyle='dotted')
        axes[3].axhline(0,color='silver',linestyle='dotted')
        axes[4].axhline(0,color='silver',linestyle='dotted')
        

        #axes[0].set_ylim([-0.01,1e8])
        #axes[2].set_ylim([-0.01,1e8])        
        #axes[4].set_ylim([-0.01,0.04])


    axes[0].legend(dirs)    
    plt.show()
