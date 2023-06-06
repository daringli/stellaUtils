#!/usr/bin/env python

import numpy as np

class Geometryfile(object):

    def __init__(self, geomfile = 'stella.geometry', geomtype = None):
        self.geomfile = geomfile
        gsplit = self.geomfile.rsplit('.',1)

        if geomtype is None:
            if len(gsplit) > 1:
                file_ending = gsplit[-1]
                if file_ending == 'nc':
                    self.geomtype = 2
                else:
                    self.geomtype = 1
            else:
                self.geomtype = 1
        else:
            self.geomtype = geomtype
            
        self.read_geometry()

        
    def read_geometry(self):

        if self.geomtype == 1:
            self.read_geometry_text()
        elif self.geomtype == 2:
            self.read_geometry_nc()


    def read_geometry_nc(self):
        from netcdf_util import netcdf_file
        with netcdf_file(self.geomfile) as f:
            self.zed =    f.variables['zed'][()]
            self.alpha =    f.variables['alpha'][()]
            self.bmag =     f.variables['bmag'][()]
            self.bdot_grad_z =  f.variables['gradpar'][()]
            self.grho  =    f.variables['grho'][()]
            self.gbdrift =  f.variables['gbdrift'][()]
            self.gbdrift0 = f.variables['gbdrift0'][()]
            self.cvdrift =  f.variables['cvdrift'][()]
            self.cvdrift0 = f.variables['cvdrift0'][()]
            self.gds2 =     f.variables['gds2'][()]
            self.gds21 =    f.variables['gds21'][()]
            self.gds22 =    f.variables['gds22'][()]
            # scalars
            #self.kxfac =    f.variables['kxfac'][()]
            self.q =        f.variables['q'][()]
            self.scale =    f.variables['scale'][()]
            #self.Rmaj =     f.variables['Rmaj'][()]
            self.shat =     f.variables['shat'][()]
            self.drhodpsi = f.variables['drhodpsi'][()]
            # metadata
            self.wout =     f.variables['wout'].filename
                                

        
        self.Nalpha = len(self.alpha)
        self.Nz = len(self.zed)
        self.zeta = np.zeros_like(self.zed) # not strictly needed
    
            
    def read_geometry_text(self):
        # states
        # 0: look for header for scalars
        # 1: look for scalars
        # 2: look for header for arrays
        # 3: look for arrays
        state = 0 
        with open(self.geomfile, 'r') as f:
            lines = f.readlines()
        for l in lines:
            L = l.strip()
            if len(L) == 0:
                continue

            if state == 0:
                if L[0] == '#':
                    if 'rhoc' in L:
                        header = L[1:].split()
                        state = 1
                        Nh = len(header)
                        
            elif state == 1:
                if L[0] == '#':
                    Ls = L[1:].split()
                    if len(Ls) == Nh:
                        for h,val in zip(header,Ls):
                            self.__dict__[h] = float(val)
                        state = 2
                        
            elif state == 2:
                if L[0] == '#':
                    Ls = L.split()
                    if 'alpha' in L:
                        header2 = L[1:].split()
                        Nh2 = len(header2)
                        state = 3
                        for h in header2:
                            self.__dict__[h] = []
                        
            elif state == 3:
                Ls = L.split()
                if len(Ls) == Nh2:
                    for h,val in zip(header2,Ls):
                        self.__dict__[h].append(float(val))

        if state == 3:
            for h in header2:
                self.__dict__[h] = np.array(self.__dict__[h])

    @property
    def gradpar(self):
        return self.bdot_grad_z
    
    def __str__(self):
        """placeholder-ish for debugging"""
        ret = str(self.zed) + "\n"
        ret = ret + str(self.rhotor)
        return ret


    def write_output(self, output_filename):
        # below is the corresponding fortran code we should mimic
        # (except we want more digits)
        # write (geometry_unit,*) is a newline:  https://pages.mtu.edu/~shene/COURSES/cs201/NOTES/chap02/write-1.html
        # (logic: each write is a newline, so an empty write is an empty newline)
#      do ia = 1, nalpha
#         do iz = -nzgrid, nzgrid
#            write (geometry_unit, '(15e12.4)') alpha(ia), zed(iz), zeta(ia, iz), bmag(ia, iz), b_dot_grad_z(ia, iz), &
#               gds2(ia, iz), gds21(ia, iz), gds22(ia, iz), gds23(ia, iz), &
#               gds24(ia, iz), gbdrift(ia, iz), cvdrift(ia, iz), gbdrift0(ia, iz), &
#               bmag_psi0(ia, iz), btor(iz)
#         end do
#         write (geometry_unit, *) 
#      end do


        
if __name__ == "__main__":
    import sys
    argc = len(sys.argv)
    if argc > 1:
        filename = sys.argv[1]
    
    g = Geometryfile(filename)
    g.bmag[:] = 0.0
    print(g)
    print(g.gradpar)
    
