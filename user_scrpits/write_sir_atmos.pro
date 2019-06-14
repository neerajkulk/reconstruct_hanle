; IDL script to write reconstructed sir files into RH
; usage: @write_sir_atmos.pro

a = read3datmos('snapshot_B_0d0_252x126x252_le_00100.atmos')

sir_temp = readfits('recons_t.fits')

sir_vx = dblarr(a.Nx,a.Ny,a.Nz)
sir_vy = dblarr(a.Nx,a.Ny,a.Nz)
sir_vz = readfits('recons_vz.fits')
sir_vturb = readfits('recons_vturb.fits')

sir_nh = readfits('recons_nh.fits')
sir_nelec = readfits('recons_nelec.fits')


openw, unit, /GET_LUN, 'sir_reconstruct.atmos', /XDR
writeu, unit, long([a.Nx,a.Ny,a.Nz,a.Nhydr]) 
writeu, unit, long([1, 2])
writeu, unit, double(a.dx), double(a.dy)
writeu, unit, double(a.z)               
writeu, unit, double(sir_temp)
writeu, unit, double(sir_nelec)
writeu, unit, double(sir_vturb) 
writeu, unit, double(sir_vx)    
writeu, unit, double(sir_vy)
writeu, unit, double(sir_vz)
writeu, unit, double(sir_nh)
free_lun, unit
