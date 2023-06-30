#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt

from simsopt.mhd.vmec_diagnostics import vmec_fieldlines
from simsopt.mhd import Vmec

import sys
from stella_input import Stella_input

from scipy.integrate import cumtrapz


plot_phi = False

if len(sys.argv) > 1:
    dirnames = sys.argv[1:]
else:
    dirnames = ['.']

variables = ['modB', 'B_sup_theta_pest', 'B_sup_phi', 'B_cross_grad_B_dot_grad_alpha', 'B_cross_grad_B_dot_grad_psi',
                     'B_cross_kappa_dot_grad_alpha', 'B_cross_kappa_dot_grad_psi',
                     'grad_alpha_dot_grad_alpha', 'grad_alpha_dot_grad_psi', 'grad_psi_dot_grad_psi',
                     'bmag', 'gradpar_theta_pest', 'gradpar_phi', 'gbdrift', 'gbdrift0', 'cvdrift', 'cvdrift0', 'gds2', 'gds21', 'gds22']

variables = ['B_cross_kappa_dot_grad_alpha', 'grad_alpha_dot_grad_alpha', 'modB']
#variables = ['grad_alpha_dot_grad_alpha']
#variables = ['gradpar_phi']
                     
nrows = np.min((len(variables),4))
ncols = len(variables)//nrows

fig,axes = plt.subplots(nrows=nrows,ncols=ncols, squeeze=False,sharex=True)
axes = axes.flatten()

for d in dirnames:
    si = Stella_input(d)
    wout_filename = d +"/" + si.vmec_filename
    
    vmec = Vmec(wout_filename, mpi=None)
    s = si.torflux
    alpha0 = si.alpha0

    
    npol = 6
    ntheta = 999
    theta1d =  np.linspace(-npol * np.pi, npol * np.pi, num=ntheta)

    f = vmec_fieldlines(vmec, s, alpha0, theta1d=theta1d, phi_center=0, plot=False, show=False)

    phi = f.phi[0,0]
    Nphi = len(phi)
    if Nphi % 2 == 0:
        pass
    else:
        i0 = Nphi//2 + 1
        Isort = np.argsort(phi)
        phi = phi[Isort]
        gradpar_phi = f.gradpar_phi[0,0,Isort]
        L = np.zeros(Nphi)
        dphi = phi[1] - phi[0]
        L[(i0-1)::-1] = cumtrapz(dphi/gradpar_phi[i0::-1],phi[i0::-1])
        L[(i0+1):] = cumtrapz(dphi/gradpar_phi[i0:],phi[i0:])

        if f.phi[0,0,0] > 0:
            L = L[::-1]

        #print(L)
        #print(L[Nphi//2])

        istart = np.argmin(np.fabs(L + 1.25))
        istop = np.argmin(np.fabs(L - 1.25))

        #np.sum(f.B_cross_kappa_dot_grad_alpha[istart:istop])/L[istart:istop]
        if istart > istop:
            istart,istop = istop,istart
        avg_curv = np.mean(f.B_cross_kappa_dot_grad_alpha[0,0,istart:istop])
        print(d)
        print(-avg_curv)
        
        
    if plot_phi:
        for j, variable in enumerate(variables):
            if variable == 'B_cross_kappa_dot_grad_alpha':
                setattr(f,variable, -eval("f." + variable))
            axes[j].plot(f.phi[0, 0, :], eval("f." + variable + '[0, 0, :]'))
            axes[j].set_xlabel('Standard toroidal angle $\phi$')
            axes[j].set_title(variable)
    else:
        for j, variable in enumerate(variables):
            if variable == 'B_cross_kappa_dot_grad_alpha':
                setattr(f,variable, -eval("f." + variable))
            axes[j].plot(L, eval("f." + variable + '[0, 0, :]'))
            axes[j].set_xlabel('Arclenght $l(\phi)$')
            axes[j].set_title(variable)
        
        #plt.figtext(0.5, 0.995, f's={f.s[0]}, alpha={f.alpha[0]}', ha='center', va='top')
        
plt.tight_layout()
axes[-1].legend(dirnames)   
plt.show()
