#!/usr/bin/env python

from netcdf_util import netcdf_file


import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib import cm
from mpl_toolkits.axes_grid1 import make_axes_locatable
import numpy as np


import sys

BIGNUM = 1e200

class ClickHandler(object):
    def __init__(self, stella_output = 'stella.out.nc'):
        self.imax_pressed = False
        self.imin_pressed = False
        self.maxlines = []
        self.minlines = []
        self.imin = 0
        self.imax = None
        with netcdf_file(stella_output,'r',mmap=False) as f:
            self.phi2_vs_kxky  = f.variables['phi2_vs_kxky'][()]
            self.ky = f.variables['ky'][()]
            self.kx = f.variables['kx'][()]
            self.t = f.variables['t'][()]
        self.gamma = np.zeros_like(self.ky)
        self.Nkx = len(self.kx)
        self.Nky = len(self.ky)
        self.phi2s = np.transpose(self.phi2_vs_kxky[:,0,:])
        self.iinfs = [None] * self.Nky
        
        for iky, phi2 in enumerate(self.phi2s):
            iinf = np.where(phi2 > BIGNUM)[0]
            if len(iinf) == 0:
                iinf = np.where(np.isnan(phi2))[0]
            if len(iinf) > 0:
                self.iinfs[iky] = iinf[0]
        
        if self.Nkx != 1:
            raise ValueError("Nkx != 1 is not supported.")

        
        
    def plot(self):
        self.fig, self.axes = plt.subplots(2)
        cid = self.fig.canvas.mpl_connect('button_press_event', self.onclick)
        cmap = cm.get_cmap("rainbow",self.Nky)
        for iky,aky in enumerate(self.ky):
            y  = self.phi2s[iky][:self.iinfs[iky]]
            x = self.t[:self.iinfs[iky]]
            self.axes[0].semilogy(x,y,color=cmap(iky))

        self.axes[0].set_xlabel(r"$t$")
        self.axes[0].set_ylabel(r"$|\phi|_k^2$")
        norm = mpl.colors.Normalize(vmin=self.ky[0], vmax=self.ky[-1])
        divider = make_axes_locatable(self.axes[0])
        cax = divider.append_axes('right', size='5%', pad=0.05)
        cb1 = mpl.colorbar.ColorbarBase(cax, cmap=cmap,
                                            norm=norm,
                                            orientation='vertical')

        #cb1.set_label(r'$k_y \rho$', labelpad=-40, y=1.05, rotation=0)
        cax.set_title(r'$k_y \rho$')
        self.gamma_line = self.axes[1].plot(self.ky,self.gamma)
        self.gamma_line = self.gamma_line[0]
        self.axes[1].set_xlabel(r"$k_y \rho$")
        self.axes[1].set_ylabel(r"$\gamma$")
        self.axes[1].set_ylim(bottom=0.0)

        #plt.show()
        plt.show()


    def onclick(self,event):
        if event.inaxes is not None:
            ax = event.inaxes
        else:
            print("Must press the topmost axis!")
            return
        if ax == self.axes[0]:
            self.gamma = []
            clicked_t = event.xdata
            if str(event.button) == "MouseButton.LEFT":
                minline = self.axes[0].axvline(clicked_t)
                self.imin  = (np.abs(self.t-clicked_t)).argmin()
                self.imin_pressed = True
                self.minlines.append(minline)
                if len(self.minlines) > 1:
                    self.minlines[0].remove()
                    del self.minlines[0]
                print("tmin: " + str(self.t[self.imin]))
                
            elif str(event.button) == "MouseButton.RIGHT":
                maxline = self.axes[0].axvline(clicked_t,color='red')
                self.imax  = (np.abs(self.t-clicked_t)).argmin()
                self.imax_pressed = True
                self.maxlines.append(maxline)
                if len(self.maxlines) > 1:
                    self.maxlines[0].remove()
                    del self.maxlines[0]
                print("tmax: " + str(self.t[self.imax]))

            # plotting
            for iky, phi2 in enumerate(self.phi2s):
                if self.imax_pressed:
                    if (self.iinfs[iky] is not None) and (self.imax > self.iinfs[iky]):
                        imax = self.iinfs[iky]
                    else:
                        imax = self.imax
                else:
                    imax = self.iinfs[iky]
                    
                y = np.log(phi2[self.imin:imax])
                x = self.t[self.imin:imax]
                # number of non-nan:
                N = np.count_nonzero(~np.isnan(y))
                if N >= 3:
                    coefs0 = np.polyfit(x,y,1,full=False)[0]
                else:
                    # no growth rate defined
                    print("Warning: no growth rate for iky = " + str(iky))
                    coefs0 = np.nan
                self.gamma.append(coefs0/2)
            self.gamma = np.array(self.gamma)

            self.gamma_line.set_ydata(self.gamma)
            #update limits must be done manually
            self.axes[1].relim()
            self.axes[1].autoscale()
            self.axes[1].set_ylim(bottom=0.0)
            plt.draw()
            np.save('my_growthrate.npy',self.gamma)
            np.save('my_ky.npy',self.ky)
            np.save("my_tfit.npy",self.t[self.imin])
            plt.savefig("phi2kykx_vs_t.png")


ch = ClickHandler('stella.out.nc')
ch.plot()
