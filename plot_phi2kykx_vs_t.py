#!/usr/bin/env python

from netcdf_util import netcdf_file


import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib import cm
from mpl_toolkits.axes_grid1 import make_axes_locatable
import numpy as np


import sys

BIGNUM = 1e200

class stellaClickHandler(object):
    def __init__(self, stella_output = 'stella.out.nc', full_kx = False, Ntheta0 = 0):
        self.imax_pressed = False
        self.imin_pressed = False
        self.maxlines = []
        self.minlines = []
        self.imin = 0
        self.imax = None
        self.full_kx = full_kx
        self.Ntheta0 = Ntheta0
        with netcdf_file(stella_output,'r',mmap=False) as f:
            self.phi2_vs_kxky  = f.variables['phi2_vs_kxky'][()] # t, kx, ky
            self.ky = f.variables['ky'][()]
            self.kx = f.variables['kx'][()]
            self.t = f.variables['t'][()]
            self.shat = f.variables['shat'][()]

        self.Nkx = len(self.kx)
        self.Nky = len(self.ky)
        print("t = " + str(self.t))

        self.phi2s = np.transpose(self.phi2_vs_kxky,axes=[1,2,0])  # kx, ky, t
        if not self.full_kx:
            # remove negative kx
            self.phi2s = self.phi2s[:(self.Nkx//2+1)]
            self.kx = self.kx[:(self.Nkx//2+1)]
            self.Nkx = self.Nkx//2 + 1
            self.kx0_index = 0
        else:
            I = np.argsort(self.kx)
            self.kx = self.kx[I]
            self.phi2s = self.phi2s[I]
            self.kx0_index = self.Nkx//2 + 1
            
        self.gamma = np.zeros((self.Nkx,self.Nky))
        self.iinfs = np.full((self.Nkx,self.Nky),None)
        
        for ikx, phi2_x in enumerate(self.phi2s):
            for iky, phi2 in enumerate(phi2_x):
                print("kx = " + str(self.kx[ikx]))
                print("ky = " + str(self.ky[iky]))
                print("phi_kxky(t) = " + str(phi2))
                iinf = np.where(phi2 > BIGNUM)[0]
                if len(iinf) == 0:
                    iinf = np.where(np.isnan(phi2))[0]
                if len(iinf) > 0:
                    self.iinfs[ikx,iky] = iinf[0]
        
        
    def plot(self):
        if self.Nkx == 1:
            self.fig, self.axes = plt.subplots(2)
        else:
            self.fig, self.axes = plt.subplots(3)
        cid = self.fig.canvas.mpl_connect('button_press_event', self.onclick)
        cmap = cm.get_cmap("rainbow",self.Nky)
        linestyles = ['solid'] + [(0, (i*0.5, i)) for i in range(1,self.Nkx)]
        for ikx,akx in enumerate(self.kx):
            for iky, aky in enumerate(self.ky):
                y  = self.phi2s[ikx,iky][:self.iinfs[ikx,iky]]
                x = self.t[:self.iinfs[ikx,iky]]
                self.axes[0].semilogy(x,y,color=cmap(iky),linestyle = linestyles[ikx])

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
        if self.Nkx == 1:
            self.gamma_plot = self.axes[1].plot(self.ky,self.gamma[0])
            self.gamma_plot = self.gamma_plot[0]
            
            self.axes[1].set_xlabel(r"$k_y \rho$")
            self.axes[1].set_ylabel(r"$\gamma$")
            self.axes[1].set_ylim(bottom=0.0)
        else:
            self.gamma_plot = self.axes[1].pcolormesh(self.ky,self.kx,self.gamma)
            divider = make_axes_locatable(self.axes[1])
            self.cax = divider.append_axes('right', size='5%', pad=0.05)
            #self.cbar = plt.colorbar(self.gamma_plot)
            #self.gamma_plot = self.gamma_plot[0]
            self.axes[1].set_xlabel(r"$k_y \rho$")
            self.axes[1].set_ylabel(r"$k_x \rho$")
            #self.axes[1].set_ylim(bottom=0.0)

            placehold = np.column_stack([self.gamma[0], self.gamma[0]] + self.Ntheta0 * [self.gamma[0]])
            self.gamma_plot2 = self.axes[2].plot(self.ky, placehold)
            self.axes[2].set_xlabel(r"$k_y \rho$")
            self.axes[2].set_ylabel(r"$\gamma$")
            self.axes[2].set_ylim(bottom=0.0)

        #plt.show()
        plt.show()


    def onclick(self,event):
        if event.inaxes is not None:
            ax = event.inaxes
        else:
            print("Must press the topmost axis!")
            return
        if ax == self.axes[0]:
            #self.gamma = []
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
            for ikx, phi2_x in enumerate(self.phi2s):
                for iky, phi2 in enumerate(phi2_x):
                    if self.imax_pressed:
                        if (self.iinfs[ikx,iky] is not None) and (self.imax > self.iinfs[ikx,iky]):
                            imax = self.iinfs[ikx,iky]
                        else:
                            imax = self.imax
                    else:
                        imax = self.iinfs[ikx,iky]
                    
                    y = np.log(phi2[self.imin:imax])
                    x = self.t[self.imin:imax]
                    # number of non-nan:
                    N = np.count_nonzero(~np.isnan(y))
                    if N >= 3:
                        coefs0 = np.polyfit(x,y,1,full=False)[0]
                    else:
                        # no growth rate defined
                        print("Warning: no growth rate for iky = " + str(iky) + ", ikx = " + str(ikx))
                        coefs0 = np.nan
                    self.gamma[ikx,iky] = coefs0/2

            if self.Nkx == 1:
                self.gamma_plot.set_ydata(self.gamma[0])
                self.axes[1].relim()
                self.axes[1].autoscale()
                self.axes[1].set_ylim(bottom=0.0)
                plt.draw()
            else:
                print(self.gamma)
                #self.gamma_plot.set_array(self.gamma.ravel())
                self.gamma_plot = self.axes[1].pcolormesh(self.ky,self.kx,self.gamma)
                vmax = np.nanmax(self.gamma)
                vmin = 0.0
                norm = mpl.colors.Normalize(vmin=vmin, vmax=vmax)
                cb1 = mpl.colorbar.ColorbarBase(self.cax, cmap='viridis', norm=norm,
                                            orientation='vertical')
                cb1.mappable.set_clim(vmin=vmin, vmax=vmax)
                cb1.draw_all()
                print(vmax)
                self.gamma_plot.set_clim(vmin=vmin,vmax=vmax)

                self.cax.set_title(r'$\gamma/(v_T/a)$')
                #self.cbar.remove()
                #self.axes[1].set_ylim(bottom=0.0)
                self.axes[1].set_xlim(left=0.0)
                #cbar_ticks = np.linspace(0., 1., num=6, endpoint=True)
                #self.cbar.ax.set_autoscale_on(True)
                #self.cbar.set_ticks(cbar_ticks)
               
                self.fig.colorbar(self.gamma_plot, cax=self.cax)

                # axes[2]
                ydata = np.nanmax(self.gamma,axis=0)
                ydata2 = self.gamma[self.kx0_index]
                self.gamma_plot2[0].set_ydata(ydata)
                self.gamma_plot2[1].set_ydata(ydata2)
                #self.axes[2].plot(self.ky,self.gamma[0])

                # create interpolation for each ky theta0
                print(self.Ntheta0)
                for j in range(1,self.Ntheta0+1):
                    _gammas = [0.0]
                    for i in range(1,self.Nky):
                        if not self.full_kx:
                            kx = j*np.fabs(2*np.pi * self.ky[i] * self.shat)
                        else:
                            kx = j*2*np.pi * self.ky[i] * self.shat
                        _gammas.append(np.interp(kx,self.kx,self.gamma[:,i],left=np.nan, right=np.nan))
                    print("HERE")
                    print(_gammas)
                    self.gamma_plot2[1+j].set_ydata(_gammas)
                    
                self.axes[2].relim()
                self.axes[2].autoscale()
                self.axes[2].set_ylim(bottom=0.0)
                self.axes[2].legend(['max','kx=0'] + ["theta0 = " + str(2*j) + "pi" for j in range(1,self.Ntheta0+1)])
                plt.draw()
                
            #update limits must be done manually
            
            
            np.save('my_growthrate.npy',self.gamma)
            np.save('my_ky.npy',self.ky)
            np.save('my_kx.npy',self.kx)
            np.save("my_tfit.npy",self.t[self.imin])
            plt.savefig("phi2kykx_vs_t.png")


if __name__=="__main__":
    import sys

    if len(sys.argv) > 1:
        dir = sys.argv[1]
    else:
        dir = '.'
    ch = stellaClickHandler(dir + '/stella.out.nc', full_kx=False)
    ch.plot()
