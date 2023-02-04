#!/usr/bin/env python

import numpy as np
from netcdf_util import netcdf_file
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


def get_best_omega(dirname, use_my = True):
    tmp = os.getcwd()
    os.chdir(dirname)

    if os.path.exists("my_growthrate.npy") and use_my:
        Imomega_avg = np.load("my_growthrate.npy")
        ky = np.load("my_ky.npy")
        kx = np.array([0])
        Imomega_avg = Imomega_avg[:,np.newaxis]
    else:
        omega_tuple = read_omega()
        t, ky, kx, Imomega, Reomega, Reomega_avg, Imomega_avg = get_latest_omega(omega_tuple)
    os.chdir(tmp)
    Imomega_avg = np.fmax(Imomega_avg,0.0)
    return Imomega_avg

def ql_flux_from_dir(dirname, average_from = None):
    tmp = os.getcwd()
    os.chdir(dirname)

    if os.path.exists("my_growthrate.npy"):
        Imomega_avg = np.load("my_growthrate.npy")
        ky = np.load("my_ky.npy")
        kx = np.array([0])
        Imomega_avg = Imomega_avg[:,np.newaxis]
    else:
        omega_tuple = read_omega()
        t, ky, kx, Imomega, Reomega, Reomega_avg, Imomega_avg = get_latest_omega(omega_tuple)


    Imomega_avg = np.fmax(Imomega_avg,0.0)
    Ky,Kx = np.meshgrid(ky,kx,indexing='ij')

    Nky = len(ky)
    Nkx = len(kx)
    with netcdf_file('stella.out.nc','r',mmap=False) as f:
        if average_from is not None:
            t = f.variables['t'][()]
            i = np.where(t>=average_from)[0][0]
            print(t[i])
            Nt = len(t[i:])
            phi2_vs_kxky  = np.transpose(np.sum(f.variables['phi2_vs_kxky'][()][i:],axis=0)/Nt)
        else:
            phi2_vs_kxky  = np.transpose(f.variables['phi2_vs_kxky'][()][-1])
        nz = f.variables['dens'][()][-1]
        # pflx_kxky(t, species, tube, zed, kx, ky)
        if average_from is not None:
            _pfluxes_vs_kxky = np.sum(f.variables['pflx_kxky'][()][i:],axis = 0)/Nt # time average
        else:
            _pfluxes_vs_kxky = f.variables['pflx_kxky'][()][-1] # last time step
        Nspecies = len(_pfluxes_vs_kxky)
        pfluxes_vs_kxky = np.zeros((Nspecies,Nkx,Nky))
        for i in range(Nspecies):
            pfluxes_vs_kxky[i] = np.sum(_pfluxes_vs_kxky[i],axis=(0,1))
        impurity_flux = pfluxes_vs_kxky[-1]/nz

        if average_from is not None:
            _qfluxes_vs_kxky = np.sum(f.variables['qflx_kxky'][()][i:],axis=0)/Nt
        else:
            _qfluxes_vs_kxky = f.variables['qflx_kxky'][()][-1]  # last time step
        qfluxes_vs_kxky = np.zeros((Nspecies,Nkx,Nky))
        for i in range(Nspecies):
            qfluxes_vs_kxky[i] = np.sum(_qfluxes_vs_kxky[i],axis=(0,1))
        energy_flux = np.sum(qfluxes_vs_kxky,axis=0)

        impurity_flux = np.transpose(impurity_flux)
        energy_flux = np.transpose(energy_flux)
    
    normalization = Imomega_avg/(phi2_vs_kxky * Ky**2) # units??
    ql_impurity_flux = normalization[:,0] * impurity_flux[:,0]
    ql_energy_flux = normalization[:,0] * energy_flux[:,0]
    os.chdir(tmp)
    return ky, ql_impurity_flux, ql_energy_flux


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

        if os.path.isfile(d + "/my_tfit.npy"):
            tfit = np.load(d + "/my_tfit.npy")
        else:
            tfit = None
        ky, ql_impurity_flux, ql_energy_flux = ql_flux_from_dir(d,tfit)
        axes[0].plot(ky, ql_impurity_flux)
        axes[1].plot(ky, ql_energy_flux)
        axes[1].set_xlabel(r"$k_y \rho$")
        axes[0].set_ylabel(r"$\Gamma_z$")
        axes[1].set_ylabel(r"$\sum_a Q_a$")

    axes[0].legend(dirs)    
    plt.show()
