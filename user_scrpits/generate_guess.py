import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d

# define an atmosphere class to read in values. 
class atmos:
    def __init__(self,file):
        raw = np.loadtxt(file,skiprows=1)
        self.tau  = raw[:,0]
        self.t    = raw[:,1]
        self.pe   = raw[:,2]
        self.vmic = raw[:,3]
        self.b    = raw[:,4]
        self.vlos = raw[:,5]
        self.inc  = raw[:,6]
        self.azi  = raw[:,7]
        self.z    = raw[:,8]
        self.p    = raw[:,9]
        self.rho  = raw[:,10]



def generate_guess(ngrid, tau_low, tau_up):
    
    # generate guess model based on FALC and extrapolating to used defined grid.
    
    falc = readmodel('/Users/neku5162/sir/models/FALC11.mod')     

    # Manually define tau grid 
    tau = np.linspace(tau_low,tau_up,ngrid,endpoint=False)

    # create an interpolating function and interpolate FALC temperatures over user defined tau grid.
    temp_interp = interp1d(falc.tau,falc.t, kind = 'linear', fill_value='extrapolate')
    temp = temp_interp(tau)
        
    pe_interp = interp1d(falc.tau,np.log(falc.pe), kind = 'linear', fill_value='extrapolate')
    pe = np.e**(pe_interp(tau))
    
    vmicro = np.zeros(ngrid)

    b_strength = np.zeros(ngrid) 

    vlos = np.zeros(ngrid) + 10.0**(5.0)

    b_inc = np.zeros(ngrid)

    b_azi = np.zeros(ngrid)
    
    # Interpolate z just as temperature is interpolated
    
    z_interp = interp1d(falc.tau,falc.z, kind = 'linear', fill_value='extrapolate')
    z = temp_interp(tau)

    p_interp = interp1d(falc.tau,np.log(falc.p), kind = 'linear', fill_value='extrapolate')
    p = np.e**(p_interp(tau))

    rho_interp = interp1d(falc.tau,np.log(falc.rho), kind = 'linear', fill_value='extrapolate')
    rho = np.e**(rho_interp(tau))
    
    guess_model = np.column_stack((tau,temp,pe,vmicro,b_strength,vlos,b_inc,b_azi,z,p,rho))
    
    # write to file
    with open('guess.mod', 'wb') as f:
        f.write(b'0.00000      0.00000      0.00000\n')
        np.savetxt(f, guess_model, delimiter=' ')
        f.close()

    
    
