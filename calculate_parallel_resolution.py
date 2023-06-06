#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import argrelmax, argrelmin

from simsopt.mhd.vmec_diagnostics import vmec_fieldlines
from simsopt.mhd import Vmec

import pdb

def remove_duplicate_extrema(y, imax, imin):
    imax = np.array(imax)
    imin = np.array(imin)
    
    # delete potentially spurious minima at the edge
    # they might be relevant, but since we only want to
    # estimate a needed resolution, they can probably be excluded.
    while imin[0] < imax[0]:
        imin = imin[1:]

    while imin[-1] > imax[-1]:
        imin = imin[:-1]

    Nmax = len(imax)
    Nmin = len(imin)
    if Nmax  == Nmin + 1:     
        mask1 = imax[:-1] < imin
        mask2 = imax[1:] > imin
        mask = np.logical_and(mask1,mask2)
        if not(mask.all()):
            do_careful = True
            print(mask1)
            print(mask2)
            print("Masks don't agree")
        else:
            filtered_imax = imax
            filtered_imin = imin
            do_careful = False
    else:
        do_careful = True

    #do_careful = True
    if do_careful:
        # this should always work, but is slower
        print("Warning! Something odd with the maxes and mins. Starting careful well construction")

        filtered_imax = list(imax)
        filtered_imin = []
        i = 0
        while i < len(filtered_imax)-1:
            i0 = filtered_imax[i]
            i1 = filtered_imax[i+1]
            #print(i1)
            #print(imin)
            indices_between = np.logical_and(imin<i1, imin>i0)
            imins_between = imin[indices_between]
            Nimins_between = len(imins_between)
            #pdb.set_trace()
            if Nimins_between > 1:
                # take lowest minimum
                # and proceed to next maximum
                _tmp = np.argmin(y[imins_between])
                filtered_imin.append(imins_between[_tmp])
                i = i + 1
            elif Nimins_between == 0:
                # there is no minimum between the maximum
                # remove the smallest of the two maximum
                # and repeat
                if y[filtered_imax[i]] > y[filtered_imax[i+1]]:
                    del filtered_imax[i+1]
                else:
                    del filtered_imax[i]
            else:
                # just 1 minimum between
                filtered_imin.append(imins_between[0])
                i = i + 1
                
        filtered_imax = np.array(filtered_imax)
        filtered_imin = np.array(filtered_imin)
            
    return filtered_imax, filtered_imin


def remove_shallow_wells(y, imax, imin, well_threshold=0.05):
    filtered_imax = list(imax)
    filtered_imin = list(imin)
    # You need to be this deep to be considered a minimum
    yrange = np.max(y) - np.min(y)
    i = 0 
    while i < len(filtered_imax) - 1:
        Dy1 = y[filtered_imax[i]] - y[filtered_imin[i]]
        Dy2 = y[filtered_imax[i+1]] - y[filtered_imin[i]]
        # smallest Dy sets the depth of well
        if Dy1 > Dy2:
            Dy = Dy2
            to_del = i + 1
        else:
            Dy = Dy1
            to_del = i

        #DBs.append(DB)
        if Dy/yrange > well_threshold:
            #filtered_imax.append(maxmin[i][0])
            i = i + 1
        else:
            del filtered_imax[to_del]
            del filtered_imin[i]
            
    imin = np.array(filtered_imin)
    imax = np.array(filtered_imax)
    imax, imin = remove_duplicate_extrema(y,imax,filtered_imin)
    return imax, imin
    
def get_Ntheta(dirname, well_threshold = 0.05, Ntheta_per_well=3, extrema_threshold = 0.05, plot=False, axis=None, show=False, wout_filename = "wout_vmec.nc"):
    wout_filename  = dirname + "/" + wout_filename
    vmec = Vmec(wout_filename, mpi=None)


    s = 0.25
    alpha = 0.0
    npol = 6
    ntheta = 2**12 + 1 # should be large and odd to have a point at theta1d = 0.0

    if ntheta % 2 == 0:
        print("Warning: ntheta should be odd to have a point at 0!")

    theta1d =  np.linspace(-npol * np.pi, npol * np.pi, num=ntheta,endpoint=True)
    Ltheta = theta1d[-1] - theta1d[0]


    f = vmec_fieldlines(vmec, s, alpha, theta1d=theta1d, phi_center=0, plot=False, show=False)
    modB = f.modB[0][0]
    Brange = np.max(modB) - np.min(modB)

    imax = np.sort(argrelmax(modB, order=5,mode='wrap')[0])
    imin = np.sort(argrelmin(modB, order=5,mode='wrap')[0])
    imax, imin = remove_duplicate_extrema(modB, imax, imin)
    imax, imin = remove_shallow_wells(modB, imax, imin, well_threshold=well_threshold)
    Nmax = len(imax)
    Nmin = len(imin)


    maxmin = list(zip(imax[:-1],imin))
    thetamaxmin = list(zip(theta1d[imax[:-1]],theta1d[imin]))
    Bmaxmin = list(zip(modB[imax[:-1]],modB[imin]))
    
    # needed later
    B_extrema = []
    for Bmm in Bmaxmin:
        B_extrema.append(Bmm[0])
        B_extrema.append(Bmm[1])
                         
    
    Dthetas = []
    theta_extrema = [] # needed later
    for i,tmm in enumerate(thetamaxmin):
        theta_extrema.append(tmm[0])
        theta_extrema.append(tmm[1])
        Dtheta = tmm[1] - tmm[0]
        if i < (Nmax-2):
            # if well is not symmetric, take the smallest radius
            Dtheta2 = thetamaxmin[i+1][0] - tmm[1]
            Dtheta = np.min([Dtheta,Dtheta2])
        Dthetas.append(Dtheta)

    Dthetas = np.array(Dthetas)

    # look at width of wells with depth over a certain percentage of the maximum well depth
    # set Ntheta resolution such that the smallest width of these well should have at least 5 grid points.
    minDtheta = np.min(Dthetas)

    # since Dtheta is the smallest "radius" of the well
    # 2*minDtheta is a lower (=stricter) bound on the well diameter
    dtheta = 2*minDtheta/Ntheta_per_well

    ntheta_lowres = int(np.ceil(Ltheta/dtheta))
    if ntheta_lowres % 2 == 0:
        ntheta_lowres = ntheta_lowres + 1

    # In principle, we should have enough theta resolution to resolve the narrowest well with 5 points.
    # However, the wells may not be well aligned with our theta grid if the extrema are between grid points
    # We increase the resolution until we capture the values of the optima within a set threshold.
    # Up to a factor of 2 at worst, since the worst case scenario had a extrema exactly between 2 grid points.
    # (We probably don't have to care about excluding shallow wells,
    #  since those would naturally have low modB-variation and thus be easy to get within threshold at lower resolution.)

    Nthetas = np.linspace(ntheta_lowres - 1, 2*(ntheta_lowres - 1),dtype=int,num=8)
    Nthetas = Nthetas + (1 - Nthetas%2)
    print(Nthetas)
    for Ntheta in Nthetas:
        theta1d_lowres =  np.linspace(-npol * np.pi, npol * np.pi, num=Ntheta,endpoint=True)
        f_lowres = vmec_fieldlines(vmec, s, alpha, theta1d=theta1d_lowres, phi_center=0, plot=False, show=False)
        modB_lowres = f_lowres.modB[0][0]


        B_extrema_lowres = np.interp(theta_extrema,theta1d_lowres, modB_lowres)
        if np.max(np.fabs(B_extrema_lowres - B_extrema)/B_extrema) < extrema_threshold:
            break

    print(Ntheta)
    if plot:
        if axis is None:
            fig, ax = plt.subplots()
        ax.plot(theta1d,modB)
        ax.plot(theta1d_lowres,modB_lowres)



        ax.plot(thetamaxmin,Bmaxmin,'o')
        if show:
            plt.show()

    return Ntheta

if __name__ == "__main__":
    import sys

    if len(sys.argv) > 1:
        dirname = sys.argv[1]
    else:
        dirname = '.'
    get_Ntheta(dirname, plot=True, show=True, well_threshold = 0.1,Ntheta_per_well=3)

    
    
