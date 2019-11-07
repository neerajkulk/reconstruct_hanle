; basic IDL routine to read line and continuum opacity of an
; atmosphere. Results are written into a .fits file. 

; note: need to run @initrh before this


@files.common
@geometry.common
@atmos.common
@opacity.common
@spectrum.common

metalFile      = "metals.out"
moleculeFile   = "molecules.out"
opacFile       = 'opacity.out'
backgroundFile = 'background.dat'

result = readInput('input.out')
result = readgeometry('geometry.out')
result = readatmos('atmos.out')
result = readspectrum('spectrum.out')
result = openOpacity(opacFile)

opac_l = fltarr(geometry.Nx, geometry.Ny, geometry.Nz, spectrum.nspect)
opac_c = fltarr(geometry.Nx, geometry.Ny, geometry.Nz, spectrum.nspect)
emis_l = fltarr(geometry.Nx, geometry.Ny, geometry.Nz, spectrum.nspect)
emis_c = fltarr(geometry.Nx, geometry.Ny, geometry.Nz, spectrum.nspect)

FOR l = 0, spectrum.nspect-1 DO BEGIN
   readOpacity, l, 0
   opac_l[*, *, *, l] = chi_as
   opac_c[*, *, *, l] = chi_c
   emis_l[*, *, *, l] = eta_as
   emis_c[*, *, *, l] = eta_c   
ENDFOR

mwrfits, opac_l, 'opac_l.fits'
mwrfits, opac_c, 'opac_c.fits'
mwrfits, emis_l, 'emis_l.fits'
mwrfits, emis_c, 'emis_c.fits'

end
