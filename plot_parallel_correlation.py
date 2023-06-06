import matplotlib.pyplot as plt
import numpy as np

from plot_phi2_vs_z import get_latest_phi2

def parallel_correlation(dirname):
    
     z , phi2 =get_latest_phi2(dirname)

     phi2avg = np.trapz(phi2,z)/(z[-1]-z[0])

     A = phi2 - phi2avg

     N = len(A)
     autocor = np.correlate(A,A, mode="full")[N-1:]
     autocor = autocor/autocor[0]
     return z - z[0], autocor
     
if __name__ == "__main__":
    import sys
    argv = sys.argv
    argc = len(argv)
    if argc > 1:
        dirs = argv[1:]
    else:
        dirs = ['.']

    fig, ax =plt.subplots(1,sharex=True)
    for d in dirs:
        z,autocorr = parallel_correlation(d)
        
        #ax.plot(z,np.fabs(autocorr))
        ax.plot(z,autocorr)
        ax.axhline(0,color='silver')
        ax.axhline(-np.exp(-2),color='silver',linestyle='dotted')
        ax.axhline(np.exp(-2),color='silver',linestyle='dotted')
        
        ax.set_xlabel(r"$z$ distance")
        ax.set_ylabel(r"$|\phi|^2$ autocorrelation")
    ax.legend(dirs)
    plt.show()
