#!/usr/bin/env python
# coding: utf-8

# # Script to take a bunch of SIR files for individual pixels and turn them into a atmosphere cube. 

# Make sure to run this outside the main pixels directory paying attention to the idexing of the pixels.




import numpy as np
import matplotlib.pyplot as plt
import glob
from astropy.io import fits




def read_pix(i,j):
    mod_names = np.sort(glob.glob('pixels/'+str(i)+'_'+str(j)+'/*.mod'))
    return np.loadtxt(mod_names[-1],skiprows=1)



class atmos:
    def __init__(self):
        self.tau  = np.array([])
        self.t    = np.array([])
        self.pe   = np.array([])
        self.vmic = np.array([])
        self.b    = np.array([])
        self.vlos = np.array([])
        self.inc  = np.array([])
        self.azi  = np.array([])
        self.z    = np.array([])
        self.p    = np.array([])
        self.rho  = np.array([])




class pixel:
    def __init__(self,i,j):
        self.tau  = read_pix(i,j)[:,0]
        self.t    = read_pix(i,j)[:,1]
        self.pe   = read_pix(i,j)[:,2]
        self.vmic = read_pix(i,j)[:,3]
        self.b    = read_pix(i,j)[:,4]
        self.vlos = read_pix(i,j)[:,5]
        self.inc  = read_pix(i,j)[:,6]
        self.azi  = read_pix(i,j)[:,7]
        self.z    = read_pix(i,j)[:,8]
        self.p    = read_pix(i,j)[:,9]
        self.rho  = read_pix(i,j)[:,10]
        




cube = atmos()
for i in np.arange(252):
    for j in np.arange(252):
        pix = pixel(i,j)
        cube.tau = np.concatenate((cube.tau,pix.tau)) 
        cube.t = np.concatenate((cube.t,pix.t)) 
        cube.pe = np.concatenate((cube.pe,pix.pe)) 
        cube.vmic = np.concatenate((cube.vmic,pix.vmic)) 
        cube.b = np.concatenate((cube.b,pix.b)) 
        cube.vlos = np.concatenate((cube.vlos,pix.vlos)) 
        cube.inc = np.concatenate((cube.inc,pix.inc)) 
        cube.azi = np.concatenate((cube.azi,pix.azi)) 
        cube.z = np.concatenate((cube.z,pix.z)) 
        cube.p = np.concatenate((cube.p,pix.p)) 
        cube.rho = np.concatenate((cube.rho,pix.rho)) 
        



cube.t = np.transpose(cube.t.reshape(252,252,64),axes = [1,0,2])
cube.pe = np.transpose(cube.pe.reshape(252,252,64),axes = [1,0,2])
cube.vmic = np.transpose(cube.vmic.reshape(252,252,64),axes = [1,0,2])
cube.vlos = np.transpose(cube.vlos.reshape(252,252,64),axes = [1,0,2])
cube.z = np.transpose(cube.z.reshape(252,252,64),axes = [1,0,2])
cube.p = np.transpose(cube.p.reshape(252,252,64),axes = [1,0,2])
cube.rho = np.transpose(cube.rho.reshape(252,252,64),axes = [1,0,2])




hdu = fits.PrimaryHDU(cube.t)
hdu.writeto('sir_t.fits')

hdu = fits.PrimaryHDU(cube.pe)
hdu.writeto('sir_pe.fits')

hdu = fits.PrimaryHDU(cube.vmic)
hdu.writeto('sir_vmic.fits')

hdu = fits.PrimaryHDU(cube.vlos)
hdu.writeto('sir_vlos.fits')

hdu = fits.PrimaryHDU(cube.z)
hdu.writeto('sir_z.fits')

hdu = fits.PrimaryHDU(cube.p)
hdu.writeto('sir_p.fits')

hdu = fits.PrimaryHDU(cube.rho)
hdu.writeto('sir_rho.fits')




