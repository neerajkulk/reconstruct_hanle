'''
Take a small cutout of a MURAM cube and save it as pickled atmosphere and an RH atmosphere.
'''

import xdrlib
import numpy as np
import muram as muram
import os.path
from os import path
import pickle

# Global variables

# dimensions of smaller cube
nx = 400
ny = 400

bottom_cut = 550  # cells to cut in bottom
top_cut = 64  # cells to cut in top
stride = 5

sim_dir = '../'  # path of MURaM directory
iter = 30000    # snaphot number

B_unit = np.sqrt(4*np.pi)
run = muram.MuramSnap(sim_dir, iter)
print('found', run.available)

NX, NZ, NY = run.T.shape  # dimensions of original MURaM cube
X, Z, Y = run.T.X/(1e5)  # in km
DX,DZ,DY = run.T.dX/1e5 # in km 

AMU = 1.6726219e-27  # in kg
mu = 1.3  # mean molecular weight


def create_z_grid(run):
    # construct z grid so z=0 <-> log tau=0
    tau_avg = np.mean(run.tau,axis=(0,2))
    z0_index = np.argmin(np.abs(np.log10(tau_avg)))
    NX, NZ, NY = run.T.shape  # dimensions of original MURaM cube

    start = -z0_index
    end = -z0_index + NZ
    z = (np.arange(start,end)*DZ)
    return z 
    

def create_atmos(run):
    atmos = {}
    params = ['rho', 'T', 'vx', 'vy', 'vz','Bx','By','Bz']
    z = create_z_grid(run)
    for param in params:
        cut_atmos = (getattr(run, param))[:, bottom_cut:-top_cut, :]
        atmos[param] = cut_atmos
        atmos[param] = atmos[param].transpose(1,2,0) #nz,ny,nx
        atmos[param] = np.flip(atmos[param],axis=0)

    z = z[bottom_cut:-top_cut] # cut z grid too 
    z = np.flip(z,axis=0)
    
    # downsample with stride and take horizontal cutout
    for param in params:
        atmos[param] = atmos[param][::stride,:ny,:nx]
    z = z[::stride]
 
    atmos['vx'] = atmos['vx']/1e5  # in km/s
    atmos['vy'] = atmos['vy']/1e5  # in km/s
    atmos['vz'] = atmos['vz']/1e5  # in km/s
    atmos['rho'] = atmos['rho']*1000  # in kg/m^3
    atmos['Bx'] = atmos['Bx']*B_unit
    atmos['By'] = atmos['By']*B_unit
    atmos['Bz'] = atmos['Bz']*B_unit

    #orient RH atmos
            
    atmos['z'] = z
    atmos['nh'] = atmos['rho']/(AMU*mu)
    atmos['DX'] = DX
    atmos['DY'] = DY

    f = open("atmos.pkl","wb")
    pickle.dump(atmos,f)
    f.close()
    
    return atmos


atmos = create_atmos(run)
# Also to write to RH atmosphere
NHydr = 1
nz,ny,nx = atmos['T'].shape

dx = DX
dy = DY

z_vals = atmos['z']

BcIrradiated = 0
BcZero = 1
BcThermalised = 2

upperBc = BcZero
lowerBc = BcThermalised

vmic = np.zeros((nz, ny, nx))
nelec = np.zeros((nz, ny, nx))


# NOTE(cmo): Pack the data into an xdrlib Packer
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
p.pack_farray(nz*ny*nx, atmos['T'].flatten(), p.pack_double)
p.pack_farray(nz*ny*nx, nelec.flatten().flatten(), p.pack_double)
p.pack_farray(nz*ny*nx, vmic.flatten(), p.pack_double)
p.pack_farray(nz*ny*nx, atmos['vx'].flatten(), p.pack_double)
p.pack_farray(nz*ny*nx, atmos['vy'].flatten(), p.pack_double)
p.pack_farray(nz*ny*nx, atmos['vz'].flatten(), p.pack_double)
p.pack_farray(NHydr*nz*ny*nx,atmos['nh'].flatten(), p.pack_double)

os.system('rm *npy')
# NOTE(cmo): Dump the packer's data into a binary file
with open('hires_deep.atmos', 'wb') as f:
    f.write(p.get_buffer())


