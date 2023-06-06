#!/usr/bin/env python

import numpy as np

from get_Qi_t import get_Qi_t_nc

twindow = 50
marker = None

mthreshold = 0.1
sigmathreshold = 0.1
autocor_threshold = np.exp(-2)




def linear_regression_window(t,Q,twindow=50):
    Nt = len(t)
    nwindow = np.argmin(np.fabs(t-twindow))

    ms = np.zeros(Nt-nwindow)
    cs = np.zeros(Nt-nwindow)
    avgQ = np.zeros(Nt-nwindow)
    rmsd = np.zeros(Nt-nwindow)
    tstarts = t[:(Nt-nwindow)]
    
    for i in range(Nt-nwindow):
        _t = t[i:(nwindow+i)]
        _t = (_t - _t[0])/(_t[-1] - _t[0])
        _Q = Q[i:(nwindow+i)]
        avgQ[i] = np.trapz(_Q,_t)
        _Q = _Q/avgQ[i]
        # we fit a linear relation
        #     Q = mt + c
        # based on  Q = [{t}, {1}] * [m,c]
        A = np.vstack([_t, np.ones(nwindow)]).T
        ret = np.linalg.lstsq(A, _Q, rcond=None)
        m,c = ret[0]
        residuals = ret[1]

        #rmsd[i] = np.sqrt(np.sum(residuals)/(nwindow - 2))
        rmsd[i] = np.sqrt(np.sum((_Q - 1)**2)/(nwindow - 2))
        ms[i] = m
        cs[i] = c
    return tstarts, ms, cs, rmsd, avgQ

def Q_linear_regression_in_dir(dirname,ncfile='gx.nc', twindow=50):
    Qi, t = get_Qi_t_nc(dirname, ncfile)
    if not np.isnan(Qi).all():
        tstarts, ms, cs, rmsd, avgQ = linear_regression_window(t,Qi,twindow)
        ret = tstarts, ms, cs, rmsd, avgQ
    else:
        ret = (np.nan, np.nan, np.nan, np.nan, np.nan)
    return ret

def get_end_of_transient_time(dirname, ncfile='gx.nc', twindow=50, mthreshold=mthreshold, sigmathreshold=sigmathreshold, plot = False, axes=None, color='b'):
    tstarts, ms, cs, rmsd, avgQ = Q_linear_regression_in_dir(dirname,twindow=twindow)
    if not np.isnan(ms).all():
        i_m = np.argmax(ms < mthreshold)
        i_sigma = np.argmax(rmsd < sigmathreshold)
        t_m = tstarts[i_m]
        t_sigma = tstarts[i_sigma]

        t_f = np.max((t_m,t_sigma))

        if plot:
            if axes is None:
                fig, axes = plt.subplots(3,sharex=True)
            
            axes[0].plot(tstarts,ms,marker=marker,color=cmap(i))
            axes[1].plot(tstarts,cs,marker=marker,color=cmap(i))
            axes[2].plot(tstarts,rmsd,marker=marker,color=cmap(i))

            axes[0].axvline(t_m, color=cmap(i), label="_nolegend_", linestyle='dashed')
            axes[2].axvline(t_sigma, color=cmap(i), label="_nolegend_", linestyle='dashed')
    else:
        t_f = np.nan
    return t_f

def autocorrelation(t, Q, autocor_threshold=autocor_threshold, plot=False, ax = None, color='b'):
    avgQ = np.trapz(Q,t)/(t[-1] - t[0]) # \int_0^T dt Q/T
    q = Q - avgQ
    N = len(q)
    autocor = np.correlate(q,q, mode="full")[N-1:]
    autocor = autocor/autocor[0]
    #print(t)
    Nt = 100
    tau = 2*np.trapz(autocor[:Nt],t[:Nt]) #correlation time, assuming autocor ~ exp(-t/tau)
    #print(tau)
    if plot:
        if ax is None:
            fig, ax = plt.subplots()
        ax.plot(t,autocor,color=color)
        ax.axvline(t[0] + tau,color=color,linestyle='dashed')

    return avgQ


def batching_averages(t,Q, twindow):
    t = t - t[0]
    T = t[-1] - t[0]
    Nwindows = T/twindow
    window_remainder = Nwindows -  T//twindow

    # if Nwindows is not an integer, we shift the start time forward
    # since the start of the simulation is more dubious than the end
    if (window_remainder) < 0.9:
        t0 = window_remainder * twindow
        i0 = np.argmin(np.fabs(t0 -t))
        Nwindows = int(np.floor(Nwindows))
    else:
        # Live with the final window being slightly smaller
        t0 = 0
        i0 = 0
        Nwindows = int(np.ceil(Nwindows))
        
        
    T = t[-1] - t[i0]
    Qavg = np.trapz(Q[i0:],t[i0:])/T
    
    Qwindowavg = np.zeros(Nwindows)
    
    for i in range(Nwindows):
        tstart = t0 + i * twindow
        tstop = t0 + (i+1) * twindow
        istart = np.argmin(np.fabs(tstart -t))
        istop = np.argmin(np.fabs(tstop -t))

        # Twindow should be almost identical to twindow
        # but may have some slight rounding due to non-uniform and finite time steps.
        Twindow = t[istop] - t[istart]
        #print(i)
        print("Twindow = " + str(Twindow))
        Qwindowavg[i] = np.trapz(Q[istart:istop],t[istart:istop])/Twindow

    # Qavg is the same as Qwindowavg since the average is linear
    Qwindowvar = np.sum((Qwindowavg - Qavg)**2)/(Nwindows - 1)
    return Qwindowavg, Qwindowvar

    
if __name__ == "__main__":
    import sys
    import matplotlib.pyplot as plt
    if len(sys.argv) > 1:
        ds = sys.argv[1:]
    else:
        ds = ['.']

    fig, axes = plt.subplots(3,sharex=True)

    Ndirs = len(ds)
    cmap = plt.get_cmap("brg",Ndirs)
    legend = []
    tfs = []
    Qis = []
    ts = []

    for i,d in enumerate(ds):    
        t_f = get_end_of_transient_time(d, twindow=twindow, mthreshold=mthreshold, sigmathreshold=sigmathreshold, plot = True, color=cmap(i), axes=axes)
        if not np.isnan(t_f):
            #print(d)
            #print(t_f)
            legend.append(d)
            tfs.append(t_f)

            Qi, t = get_Qi_t_nc(d)
            ts.append(t)
            Qis.append(Qi)
        else:
            ts.append(np.nan)
            Qis.append(np.nan)
            tfs.append(np.nan)
    
    axes[0].axhline(mthreshold,color='silver',label="_nolegend_",linestyle='dotted')
    axes[1].axhline(1,color='silver',label="_nolegend_",linestyle='dotted')
    axes[2].axhline(sigmathreshold,color='silver',label="_nolegend_",linestyle='dotted')
    axes[-1].set_xlabel(r'start $t/[a/v_{\rm{Ti}}]$')
    axes[0].set_ylabel(r'$m$')
    axes[1].set_ylabel(r'$c$')
    axes[2].set_ylabel(r'$\sigma$')
    
    axes[0].legend(legend)
    axes[0].set_title(r"$\hat{Q} = m\hat{t} + c$, window size $50$")

    fig2, axes2 = plt.subplots(2,sharex=True)
    for i,d in enumerate(ds):
        if not np.isnan(tfs[i]):
            axes2[0].plot(ts[i],Qis[i],color=cmap(i))
            axes2[0].axvline(tfs[i],color=cmap(i), label="_nolegend_", linestyle='dashed')
            i_tf = np.argmax(ts[i]>tfs[i])
            #print(i_tf)
            #print(ts[i][i_tf:])
            autocorrelation(ts[i][i_tf:],Qis[i][i_tf:], ax=axes2[1],plot=True,color=cmap(i),autocor_threshold=autocor_threshold)
    
    axes2[0].set_ylabel(r'$Q_{\rm{i}}/[Q_{\rm{gB}}]$')
    axes2[1].set_ylabel(r'autocorrelation')

    axes2[1].axhline(autocor_threshold,color='silver', label="_nolegend_", linestyle='dotted')
    axes2[1].axhline(0.0,color='silver', label="_nolegend_", linestyle='dotted')
    
    axes2[-1].set_xlabel(r'$t/[a/v_{\rm{Ti}}]$')
    
    
    axes2[0].legend(legend)

    fig3, axes3 = plt.subplots(2,sharex=True)
    Ntwindows = 150
    ns = np.arange(2,Ntwindows + 2)
    
    
    for i,d in enumerate(ds):
        Qwindowvar = np.zeros(Ntwindows)
        twindows = np.zeros(Ntwindows)
        tau = np.zeros(Ntwindows)
        if not np.isnan(tfs[i]):
            i_tf = np.argmax(ts[i]>tfs[i])
            T = ts[i][-1] - ts[i][i_tf]

            Qavg = np.trapz(Qis[i][i_tf:], ts[i][i_tf:])/T
            
            Qvar = np.trapz((Qis[i][i_tf:] - Qavg)**2, ts[i][i_tf:])/T


            twindows = np.linspace(0.01*T,T/10,Ntwindows)
            for j in range(Ntwindows):
                print(j)
                Qwindowavg, Qwindowvar[j] = batching_averages(ts[i][i_tf:],Qis[i][i_tf:], twindows[j])
                tau[j] = twindows[j] * Qwindowvar[j]/Qvar
                
            axes3[0].plot(twindows/T, Qwindowvar)
            axes3[1].plot(twindows/T, tau)

        
    axes3[-1].set_xlabel(r'$t_{\rm{window}} / T$')
    axes3[0].set_ylabel(r'$Var[Q_i^{\rm{window}}]$')
    axes3[1].set_ylabel(r'$\tau$')
    
    
    axes3[0].legend(legend)

    plt.show()

    
   

    
