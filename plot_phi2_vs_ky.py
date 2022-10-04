from scipy.io import netcdf

import matplotlib.pyplot as plt

with netcdf.netcdf_file('stella.out.nc','r',mmap=False) as f:
    phi2_vs_kxky  = f.variables['phi2_vs_kxky'][()]
    ky = f.variables['ky'][()]
    kx = f.variables['kx'][()]
    

if len(kx) == 1:
    # last time-step
    phi2_vs_ky  = phi2_vs_kxky[-1,0]
    plt.plot(ky,phi2_vs_ky)


    plt.xlabel(r"$\rho k_y$")
    plt.ylabel(r"$|\phi|^2$")
    
    plt.show()
