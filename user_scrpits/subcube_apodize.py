# Take big atmosphere and turn them into subcubes and also apodize the subcubes. 
import numpy as np
import muram as muram
import matplotlib.pyplot as plt
import sys
sys.path.append('/Users/kulkarniad/utils_rh/') # RH util functions
from rh import *
from astropy.io import fits
import pickle
from scipy.signal import general_gaussian
from skimage.util.shape import view_as_windows


with open('atmos.pkl', 'rb') as f:
    atmos = pickle.load(f)

NZ, NY, NX = atmos['T'].shape # dimensions of the cutout

def is3D(array):
    return len(array.shape) == 3

def mk_subcubes(arr_in,scube_shape, step):
    NZ,NY,NX = arr_in.shape
    nz,ny,nx = scube_shape
    dz,dy,dx = step
    
    if (dz!= NZ) or (nz!= NZ):
        raise Exception('NZ=nz=dz')
    if (NX-nx)%dx != 0:
        raise Exception('(NX - nx)%step_x == 0 ----- subcubes dont overlap correctly')
    if (NY-ny)%dy != 0:
        raise Exception('(NY - ny)%step_y == 0 ----- subcubes dont overlap correctly')
    subcubes = view_as_windows(arr_in, scube_shape, step=step)[0].copy()
    return subcubes
    

for param in atmos.keys():
    if(is3D(atmos[param])):
        
        # make subcubes:
        atmos[param] = mk_subcubes(atmos[param],(NZ,210,210), (NZ,190,190)) # hardcoded, change later

        # apodize edges
        exp = 6  # steepness of gaussian smoothing
        overlap = 20 # 210 - 190 hardocded
        Ny_scubes,Nx_scubes,NZ,ny,nx = atmos[param].shape
        mid = (nx - overlap)//2
        window = np.outer(general_gaussian(nx,exp,mid),general_gaussian(ny,exp,mid))
        for j in range(Ny_scubes):
            for i in range(Nx_scubes):
                for k in range(NZ):
                    atmos[param][j,i,k] =  window*atmos[param][j,i,k] + (1 - window)*np.mean(atmos[param][j,i,k])

NHydr = 1
Ny_scubes,Nx_scubes,nz,ny,nx = atmos['T'].shape

dx = atmos['DX']
dy = atmos['DY']

z_vals = atmos['z']

BcIrradiated = 0
BcZero = 1
BcThermalised = 2

upperBc = BcZero
lowerBc = BcThermalised

vmic = np.zeros((nz, ny, nx))
nelec = np.zeros((nz, ny, nx))

# save each of of the subcubes
for j in range(Ny_scubes):
    for i in range(Nx_scubes):
        p = xdrlib.Packer()
        p.pack_int(nx)
        p.pack_int(ny)
        p.pack_int(nz)
        p.pack_int(NHydr)
        p.pack_int(upperBc)
        p.pack_int(lowerBc)

        p.pack_double(dx)
        p.pack_double(dy)

        p.pack_farray(nz, z_vals, p.pack_double)
        p.pack_farray(nz*ny*nx, atmos['T'][j,i].flatten(), p.pack_double)
        p.pack_farray(nz*ny*nx, nelec.flatten().flatten(), p.pack_double)
        p.pack_farray(nz*ny*nx, vmic.flatten(), p.pack_double)
        p.pack_farray(nz*ny*nx, atmos['vx'][j,i].flatten(), p.pack_double)
        p.pack_farray(nz*ny*nx, atmos['vy'][j,i].flatten(), p.pack_double)
        p.pack_farray(nz*ny*nx, atmos['vz'][j,i].flatten(), p.pack_double)
        p.pack_farray(NHydr*nz*ny*nx,atmos['nh'][j,i].flatten(), p.pack_double)

        # NOTE(cmo): Dump the packer's data into a binary file
        with open(f"hires_{i}_{j}.atmos", 'wb') as f:
            f.write(p.get_buffer())
    