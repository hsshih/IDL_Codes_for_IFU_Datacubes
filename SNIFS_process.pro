PRO parinfostruct__define

tmp = {parinfostruct, $
		value:0.D, $	
		fixed:0, $	
		limited:[0,0], $	
		limits:[0.D,0.D] $		
			}
END





PRO objinfostruct__define

tmp = {objinfostruct, $
		Name:'tmp',$
		Objtype:'RG or QSO',$		;Is object RG or QSO?
		Filenum: 0,$
		RA:0., $		;RA in degrees
		DEC:0., $		;DEC in degrees
		SpecIndex:0.,$
		Radio1400:0., $	
		Radio365:0., $	
		Redshift:0., $	
		Gooddata:'yes?', $
		Outflow:'no?', $
		Calc:'no', $		
		O3nuclear:0.D, $	
		O3extended:0.D, $	
		PA_cone:0., $		;Position angle from fitting bicone model in degrees
		inc_cone:0., $		;Inclination angle from fitting bicone model in degrees
		PA_disk:0., $		;Position angle from fitting disk model in degrees
		inc_disk:0., $		;Inclination angle from fitting disk model in degrees
		Guess_PA:'tmp', $
		Radio_PA:'tmp' $
			}
END





PRO objinfo

;Read in and save objinfo data structure

;READCOL, 'snifs_objinfo', objname, objtype, filenum, redshift, radio1400, radio365, specindex, ra, dec, gooddata, outflow, FORMAT='A,I,A,F,F,F,F,F,F,A,A'
READCOL, 'tmpsnifsinfo.txt', objname, filenum, objtype, redshift, radio1400, radio365, specindex, ra, dec, outflow, gooddata, calc, radio_pa, guess_pa, FORMAT='A,I,A,F,F,F,F,F,F,A,A,A,A,A'

nstruct = N_ELEMENTS(objname)
tmpobjinfo = {objinfostruct}
objinfo = REPLICATE(tmpobjinfo, nstruct)

FOR i=0,nstruct-1 DO BEGIN

	objinfo[i].Name = objname[i]
	objinfo[i].Filenum = filenum[i]
	objinfo[i].Objtype = objtype[i]
	objinfo[i].RA = ra[i]
	objinfo[i].DEC = dec[i]
	objinfo[i].SpecIndex = specindex[i]
	objinfo[i].Radio1400 = radio1400[i]
	objinfo[i].Radio365 = radio365[i]
	objinfo[i].Redshift = redshift[i]
	objinfo[i].Gooddata = gooddata[i]
	objinfo[i].Outflow = outflow[i]
	objinfo[i].Calc = calc[i]
	objinfo[i].Radio_PA = radio_pa[i]
	objinfo[i].Guess_PA = guess_pa[i]

ENDFOR

SAVE, objinfo, filename='snifs_objinfo.save'

END







FUNCTION name_extract, filename

tmpdat = READFITS(filename, filehead, /SILENT)
objname = SXPAR(filehead, 'OBJECT')

RETURN, objname

END



FUNCTION datenum_extract, cubename

extract = STRSPLIT(cubename, '_', /EXTRACT)
datenum = STRJOIN(extract[1:3], '_')
RETURN, datenum


END




PRO cube_extract, cubename, objname, xcen, ycen, radius, filenum

;Given a center and radius, extract a 'cylinder' from a data cube

datacube = READFITS(cubename, cubehead, /SILENT)
xdim = N_ELEMENTS(datacube[*,1,1])
ydim = N_ELEMENTS(datacube[1,*,1])
wavdim = N_ELEMENTS(datacube[1,1,*])

crval3 = SXPAR(cubehead, 'CRVAL3')
cd3_3 = SXPAR(cubehead, 'CD3_3')
object = SXPAR(cubehead, 'OBJECT')

collapse_cube = MEDIAN(datacube, DIMENSION=3)
CNTRD, collapse_cube, xcen, ycen, xcentr, ycentr, 2
IF (xcentr lt 0.0) THEN BEGIN
	xcentr = xcen
	ycentr = ycen
ENDIF


DIST_CIRCLE, circles, [xdim,ydim], xcentr-1, ycentr-1
goodidx = WHERE(circles LT radius)

extspec = FLTARR(wavdim)

FOR i=0,wavdim-1 DO BEGIN

tmpslice = datacube[*,*,i]
;extspec[i] = TOTAL(tmpslice[goodidx])
extspec[i] = MEDIAN(tmpslice[goodidx])

ENDFOR

new_crval1 = crval3
new_cd1_1 = cd3_3

FXADDPAR, cubehead, 'SIMPLE', 'T'
FXADDPAR, cubehead, 'BITPIX', -32
FXADDPAR, cubehead, 'NAXIS', 1
FXADDPAR, cubehead, 'OBJECT', object
FXADDPAR, cubehead, 'DISPAXIS', 1
FXADDPAR, cubehead, 'CRVAL1', new_crval1
FXADDPAR, cubehead, 'CD1_1', new_cd1_1

FXADDPAR, cubehead, 'CRVAL2', 1
FXADDPAR, cubehead, 'CD2_2', 1
FXADDPAR, cubehead, 'CRVAL3', 1
FXADDPAR, cubehead, 'CD3_3', 1

outspec = objname + "_extspec" + STRTRIM(filenum,2) + ".fits"

WRITEFITS, outspec, extspec, cubehead
;WRITEFITS, 'circle.fits', circles

END





PRO snifs_subsky, cubename, subskyname, xcen, ycen
;Use skyflat from 4C+14.11 sky (not sky subtracted) to flat field each slice 
;Use mastercorr.fits to apply correction to the flux calibration
;Subtract from each slice the median of the slice

	cube = READFITS(cubename, header, /SILENT)
	dim = SIZE(cube)
	nslice = dim[3]
	subcube = FLTARR(dim[1],dim[2],dim[3])
	skyflat = READFITS('snifs_skyflat.fits', /SILENT)

	mastercorr = READFITS('mastercorr.fits', /SILENT)
	cubesize = SIZE(cube)
	corrcube = FLTARR(cubesize[1], cubesize[2], cubesize[3])

	FOR i=0,cubesize[1]-1 DO BEGIN
	FOR j=0,cubesize[1]-1 DO BEGIN
		corrcube[i,j,*] = mastercorr
	ENDFOR
	ENDFOR

	cube = cube * corrcube


	FOR j=0,nslice-1 DO BEGIN

		slice = cube[*,*,j]
		slice = slice / skyflat
		sky = MEDIAN(slice)
		subcube(*,*,j) = slice - sky

	ENDFOR

	collapse_cube = MEDIAN(subcube, DIMENSION=3)

	;CNTRD, collapse_cube, xcen, ycen, xcentr, ycentr, 2
	;xshf = xcen - xcentr
	;yshf = ycen - ycentr

	;IF (xcentr lt 0.0) THEN BEGIN
		;xshf = 0.0
		;yshf = 0.0
	;ENDIF

	;shiftcube = FLTARR(dim[1],dim[2],dim[3])
	
	;FOR k=0,nslice-1 DO BEGIN
		;shiftcube(*,*,k) = FSHIFT(subcube[*,*,k], xshf, yshf)
	;ENDFOR

	WRITEFITS, subskyname, subcube, header
	;WRITEFITS, clspname[i], collapse_cube

END








FUNCTION line_slice, cube, line_coord

line = cube[*,*,line_coord[0]:line_coord[1]]
linesum = TOTAL(line, 3)

RETURN, linesum

END




FUNCTION flux_convert_slice, image_slice, redshift

;Original flux calibration from
;Convert flux scale on a slice to absolute surface brightness erg/s
;Require redshift input (calls, NWCC to calculate luminosity distance using Ned Wright's Javascript Cosmology Calculator in Mpc)

Ho=71
Omega_m = 0.27
Omega_vac = 0.73

result = NWCC(redshift, Ho=Ho, Omega_m=Omega_m, Omega_vac=Omega_vac)
dl = result[0] * 3.08567758 * 1.D24

app = 2.93	;Angstrom per pixel
;;arcpp = 0.4^2.0	;Arcsec^2 per pixel, if we want to convert to erg/s/arcsec^2
pi = 3.14159265359

conslice = image_slice * 4. * pi * dl^2.0 * app 	;/ arcpp

RETURN, conslice

END








PRO snifs_mko3sub, filename, subo3cubename, objname, redshift, filenum
;Get Median continuum
;; ***Make sure CRVAL3 and CD3_3 are correct in the header, if 'emptycube.fits' is copied, header need to be edited
;; 01/25/13 Modified the way veldisp is calculated such that instead of being divided by redshifted 5006.9, it is divided by the fit centroid of that line, but difference is too small to tell

data = READFITS(filename, head, /SILENT)

crval3 = SXPAR(head, 'CRVAL3')
cd3_3 = SXPAR(head, 'CD3_3')

o3_zwav = 5006.9 * (1. + redshift)
Hbeta_zwav = 4861. * (1. + redshift)

o3pix = FIX((o3_zwav - crval3) / cd3_3)
Hbetapix = FIX((Hbeta_zwav - crval3) / cd3_3)

o3_begin = o3pix - 15
o3_end = o3pix + 20

r_begin = o3_end + 10
r_end = r_begin + 100

b_end = Hbetapix - 30
b_begin = b_end - 100

bluecont  = data[*,*, b_begin:b_end]
redcont = data[*,*, r_begin:r_end]


nb = N_ELEMENTS(bluecont[1,1,*])
nr = N_ELEMENTS(redcont[1,1,*])

cont = FLTARR(15,15,nb+nr)
xc = 7
yc = 7

cont(*,*,0:nb-1) = bluecont
cont(*,*,nb:nb+nr-1) = redcont

contmed = MEDIAN(cont, DIMENSION=3)


;Subtract Continuum

new_crval3 = crval3 + (o3_begin * cd3_3)
FXADDPAR, head, 'CRVAL3', new_crval3

;radius = 2.0

cube = data[*,*, o3_begin:o3_end]

nx = N_ELEMENTS(cube[*,1,1])
ny = N_ELEMENTS(cube[1,*,1])
nslice = N_ELEMENTS(cube[1,1,*])

subcube = FLTARR(nx,ny,nslice)
scaled_cont_cube = FLTARR(nx,ny,nslice)

rmid = FIX(nr / 2) + r_begin
bmid = FIX(nb / 2) + b_begin
rmid_val = redcont[xc,yc,FIX(nr / 2)]
bmid_val = bluecont[xc,yc,FIX(nb / 2)]
o3mid = FIX(nslice / 2) + o3_begin

slope = (rmid_val - bmid_val) / (rmid - bmid)

startpix = o3_begin	;starting pixel of input cube
;midpoint = o3_end	;midpoint of the two regions used for continuum?????
midpoint = o3mid	;midpoint of the two regions used for continuum
value_mp = bmid_val + (slope * (o3mid - bmid)) ;interpolated continuum flux at midpoint

print, bmid_val, bmid, o3mid, slope, value_mp

;maxcont = MAX(contmed[20:40,40:60])

FOR i = 0,nslice-1 DO BEGIN

	tmpslice = cube[*,*,i]

	;Scale and subtract continuum

	value = value_mp + (((startpix + i) - midpoint) * slope)   ;interpolated continuum flux at given slice
	ratio = value / value_mp
	scaled_cont = contmed * ratio

	subtmp = tmpslice - scaled_cont
	subcube(*,*,i) = subtmp	
	scaled_cont_cube(*,*,i) = scaled_cont

ENDFOR

;Write Sub_O3cube file
WRITEFITS, subo3cubename, subcube, head

;Make continuum image and absolute flux-converted-combined O3 image
line_coord = [o3_begin,o3_end]
line_slice = LINE_SLICE(data, line_coord)
convert_slice = FLUX_CONVERT_SLICE(line_slice, redshift)

contname = objname + 'contimg' + STRTRIM(filenum,2) + '.fits'
line_slice_name = objname + 'o3slice' + STRTRIM(filenum,2) + '.fits'
WRITEFITS, contname, contmed
WRITEFITS, line_slice_name, convert_slice


END




PRO snifs_mko2sub, filename, subo2cubename, objname, redshift, filenum
;Get Median continuum
;; ***Make sure CRVAL3 and CD3_3 are correct in the header, if 'emptycube.fits' is copied, header need to be edited
;; 01/25/13 Modified the way veldisp is calculated such that instead of being divided by redshifted 5006.9, it is divided by the fit centroid of that line, but difference is too small to tell

data = READFITS(filename, head, /SILENT)

crval3 = SXPAR(head, 'CRVAL3')
cd3_3 = SXPAR(head, 'CD3_3')

o2_zwav = 3727. * (1. + redshift)

o2pix = FIX((o2_zwav - crval3) / cd3_3)

o2_begin = o2pix - 20
o2_end = o2pix + 20

r_begin = o3_end + 10
r_end = r_begin + 100

b_begin = o2_begin - 10
b_end = b_begin - 100

bluecont  = data[*,*, b_begin:b_end]
redcont = data[*,*, r_begin:r_end]


nb = N_ELEMENTS(bluecont[1,1,*])
nr = N_ELEMENTS(redcont[1,1,*])

cont = FLTARR(15,15,nb+nr)
xc = 7
yc = 7

cont(*,*,0:nb-1) = bluecont
cont(*,*,nb:nb+nr-1) = redcont

contmed = MEDIAN(cont, DIMENSION=3)


;Subtract Continuum

new_crval3 = crval3 + (o2_begin * cd3_3)
FXADDPAR, head, 'CRVAL3', new_crval3

;radius = 2.0

cube = data[*,*, o3_begin:o3_end]

nx = N_ELEMENTS(cube[*,1,1])
ny = N_ELEMENTS(cube[1,*,1])
nslice = N_ELEMENTS(cube[1,1,*])

subcube = FLTARR(nx,ny,nslice)
scaled_cont_cube = FLTARR(nx,ny,nslice)

rmid = FIX(nr / 2) + r_begin
bmid = FIX(nb / 2) + b_begin
rmid_val = redcont[xc,yc,FIX(nr / 2)]
bmid_val = bluecont[xc,yc,FIX(nb / 2)]
o2mid = FIX(nslice / 2) + o2_begin

slope = (rmid_val - bmid_val) / (rmid - bmid)

startpix = o2_begin	;starting pixel of input cube
midpoint = o2mid	;midpoint of the two regions used for continuum
value_mp = bmid_val + (slope * (o3mid - bmid)) ;interpolated continuum flux at midpoint

print, bmid_val, bmid, o3mid, slope, value_mp

;maxcont = MAX(contmed[20:40,40:60])

FOR i = 0,nslice-1 DO BEGIN

	tmpslice = cube[*,*,i]

	;Scale and subtract continuum

	value = value_mp + (((startpix + i) - midpoint) * slope)   ;interpolated continuum flux at given slice
	ratio = value / value_mp
	scaled_cont = contmed * ratio

	subtmp = tmpslice - scaled_cont
	subcube(*,*,i) = subtmp	
	scaled_cont_cube(*,*,i) = scaled_cont

ENDFOR

;Write Sub_O2cube file
WRITEFITS, subo2cubename, subcube, head

;Make continuum image and absolute flux-converted-combined O2 image
line_coord = [o2_begin,o2_end]
line_slice = LINE_SLICE(data, line_coord)
convert_slice = FLUX_CONVERT_SLICE(line_slice, redshift)

contname = objname + 'contimg' + STRTRIM(filenum,2) + '.fits'
line_slice_name = objname + 'o3slice' + STRTRIM(filenum,2) + '.fits'
WRITEFITS, contname, contmed
WRITEFITS, line_slice_name, convert_slice


END













PRO snifs_cube_fitting, fitsname, objname, parinfo, redshift, filenum, compnum=compnum, fix_ctr_coord=fix_ctr_coord

;compnum = 1 or 2

IF (not KEYWORD_SET(compnum)) THEN compnum = 1

npar = compnum * 3 + 1

cube = READFITS(fitsname, cubehead, /SILENT)

IF (not KEYWORD_SET(fix_ctr_coord)) THEN BEGIN
	collapse_cube = MEDIAN(cube, DIMENSION=3)
	CNTRD, collapse_cube, 7, 7, xcentr, ycentr, 2
	IF ((xcentr GT 0.0) AND (ycentr GT 0.0)) THEN BEGIN
		fix_ctr_coord=[xcentr, ycentr]
	ENDIF ELSE BEGIN
		fix_ctr_coord=[7,7]
	ENDELSE
ENDIF


READCOL, 'CoordOrder', x, y, FORMAT='I,I'

npix = N_ELEMENTS(x)

cubex = N_ELEMENTS(cube[*,1,1])
cubey = N_ELEMENTS(cube[1,*,1])
cubez = N_ELEMENTS(cube[1,1,*])

wavelength = FINDGEN(cubez)

fit_result = FLTARR(npar,npix)
fit_error = FLTARR(npar,npix)

vel_slice1 = FLTARR(cubex,cubey) - 10000.0
veldisp_slice1 = FLTARR(cubex,cubey)
flux_slice1 = FLTARR(cubex,cubey)

IF (compnum EQ 2) THEN BEGIN

vel_slice2 = FLTARR(cubex,cubey) - 10000.0
veldisp_slice2 = FLTARR(cubex,cubey)
flux_slice2 = FLTARR(cubex,cubey)

ENDIF



crval3 = SXPAR(cubehead, 'CRVAL3')
cd3_3 = SXPAR(cubehead, 'CD3_3')
o3_redshifted = 5006.9 * (1 + redshift)

txtname = objname + 'result.txt'
OPENW, 1, txtname

;Scale everything by a large constant, MPFIT can't deal with small values for some reason...
constscale = 1.e20
parinfo[3].value = parinfo[3].value * constscale
IF (compnum EQ 2) THEN parinfo[6].value = parinfo[6].value * constscale

;Extract the initial original initial fit parameter
iniparm = parinfo[*].value


FOR i=0,npix-1 DO BEGIN

	spectrum = cube[x[i],y[i],*]
	spectrum_scaled = spectrum * constscale

	modelname = objname + 'modspec' + STRTRIM(FIX(x[i]),2) +'_' + STRTRIM(FIX(y[i]),2)
	modelresult = MYFIT(wavelength, spectrum_scaled, 1, modelname, parinfo, compnum)
	IF (compnum EQ 1) THEN model = modelresult(0)+gauss1(wavelength, modelresult(1:3))
	IF (compnum EQ 2) THEN model = modelresult(0)+gauss1(wavelength, modelresult(1:3))+gauss1(wavelength, modelresult(4:6))

	;model_err = MYMC(wavelength, model, spectrum_err, parinfo, compnum)

	fit_result[*,i] = modelresult
	;fit_error[*,i] = model_err

	;Using sigma for velocity dispersion sigma = FWHM/2.355
	;Instrumentation Broadening from arc images (There were slight deviations between two gratings 3.74 vs 3.96, but close enough....)
	;sig_inst = 1.64

	vel1 = ((((cd3_3 * modelresult[1]) + crval3) - o3_redshifted) / o3_redshifted) * 3.E5
	;;;;sigma1 = SQRT((modelresult[2]^2.) - (sig_inst^2.))
	sigma1 = modelresult[2]
	IF (NOT (sigma1 GE 0.0)) THEN sigma1 = 0.0
	;;;;veldisp1 = ((sigma1 * cd3_3) / o3_redshifted) * 3.E5
	veldisp1 = ((sigma1 * cd3_3) / ((cd3_3 * modelresult[1]) + crval3)) * 3.E5

	vel_slice1[x[i],y[i]] = vel1
	veldisp_slice1[x[i],y[i]] = veldisp1	;Subtract instrumentation broadening in quadrature
	flux_slice1[x[i],y[i]] = modelresult[3] / constscale

	IF (compnum EQ 1) THEN PRINTF, 1, x[i], y[i], fit_result[*,i], FORMAT='(F9.2,F9.2,F9.2,F9.2,F9.2,F11.1)' 

	IF (compnum EQ 2) THEN BEGIN

		vel2 = ((((cd3_3 * modelresult[4]) + crval3) - o3_redshifted) / o3_redshifted) * 3.E5
		;;;;sigma2 = SQRT((modelresult[5]^2.) - (sig_inst^2.))
		sigma2 = modelresult[5]
		IF (NOT (sigma2 GE 0.0)) THEN sigma2 = 0.0
		;;;;veldisp2 = ((sigma2 * cd3_3) / o3_redshifted) * 3.E5
		veldisp2 = ((sigma2 * cd3_3) / ((cd3_3 * modelresult[4]) + crval3)) * 3.E5
		vel_slice2[x[i],y[i]] = vel2
		veldisp_slice2[x[i],y[i]] = veldisp2
		flux_slice2[x[i],y[i]] = modelresult[6] / constscale

		PRINTF, 1, x[i], y[i], fit_result[*,i], FORMAT='(F9.2,F9.2,F9.2,F9.2,F9.2,F13.1,F9.2,F9.2,F13.1)' 

	ENDIF 
	
;Set next startparm to this model
;startparm = modelresult
IF (i EQ 0) THEN firstvel = vel1

ENDFOR

CLOSE, 1


fitresultname = objname + 'fitresult_' + STRTRIM(filenum,2) + '.txt'
vel_slicename1 = objname + 'vel_slice1_' + STRTRIM(filenum,2) + '.fits'
veldisp_slicename1 = objname + 'veldisp_slice1' + STRTRIM(filenum,2) + '.fits'
flux_slicename1 = objname + 'flux_slice1' + STRTRIM(filenum,2) + '.fits'
vel_slicename2 = objname + 'vel_slice2' + STRTRIM(filenum,2) + '.fits'
veldisp_slicename2 = objname + 'veldisp_slice2' + STRTRIM(filenum,2) + '.fits'
flux_slicename2 = objname + 'flux_slice2' + STRTRIM(filenum,2) + '.fits'

zeroidx = WHERE(flux_slice1 EQ 0.0)
vel_slice1[zeroidx] = -999999.0
veldisp_slice1[zeroidx] = 0.0

;WRITECOL, fitresultname, xnode, ynode, fit_result[0,*], fit_result[1,*], fit_result[2,*], fit_result[3,*]

WRITEFITS, vel_slicename1, vel_slice1
WRITEFITS, veldisp_slicename1, veldisp_slice1
WRITEFITS, flux_slicename1, flux_slice1

IF (compnum EQ 2) THEN BEGIN

WRITEFITS, vel_slicename2, vel_slice2
WRITEFITS, veldisp_slicename2, veldisp_slice2
WRITEFITS, flux_slicename2, flux_slice2

ENDIF

END



FUNCTION duplicate, list, element

;Find out if 'element' already exits in 'list'

idx = WHERE(list EQ element)


IF (idx[0] LT 0) THEN RETURN, 0
IF (idx[0] GE 0) THEN RETURN, 1

END















PRO get_cube_fit

READCOL, 'namefile_3', objnames, redshift, cubenames, flux, filenum, FORMAT='A,F,A,F,I'	;All obs number *003*

xcen=8.0
ycen=8.0
radius=2.0

parinfo_struct = {parinfostruct}
parinfo = REPLICATE(parinfo_struct, 4)
;parinfo = REPLICATE(parinfo_struct, 7)

;;parinfo[*].value = [0.0, 12.0, 2.0, 4.e-17]	;HB0110+297
;;parinfo[*].value = [0.0, 14.0, 4.0, 2.0e-16]	;3C099

parinfo[1].limited=[1,1]
parinfo[1].limits=[-500.0, 500.0]
parinfo[1].limited=[1,1]
parinfo[1].limits=[2.0,20.0]
parinfo[2].limited=[1,1]
parinfo[2].limits=[0.01,20.0]
parinfo[3].limited=[1,0]
parinfo[3].limits=[0.0,10000000000000.0]

;;parinfo[*].value = [0.0, 60.0, 10.0, 1.e-17, 70.0, 10.0, 1.e-17]
;;parinfo[4].limited=[1,1]
;;parinfo[4].limits=[65.0,90.0]
;;parinfo[5].limited=[1,1]
;;parinfo[5].limits=[5.0,22.0]
;;parinfo[6].limited=[1,0]
;;parinfo[6].limits=[0.0,100000000.0]

ncubes = N_ELEMENTS(cubenames)

FOR i=0,ncubes-1 DO BEGIN
;FOR i=ncubes-1,ncubes-1 DO BEGIN
;FOR i=0,0 DO BEGIN
;FOR i=96,96 DO BEGIN

subskyname = objnames[i] + "_subskycube_" + STRTRIM(filenum[i],2) + ".fits"
subo3cubename = objnames[i] + "subo3cube_" + STRTRIM(filenum[i],2) + ".fits"

;;parinfo[*].value = [par1[i], par2[i], par4[i], par4[i]]
parinfo[*].value = [0.0, 16.0, 2.0, flux[i]]

xcen=8
ycen=8
radius=3

;SNIFS_SUBSKY, cubenames[i], subskyname, 8, 8
;CUBE_EXTRACT, subskyname, objnames[i], xcen, ycen, radius, filenum[i]
;SNIFS_MKO3SUB, subskyname, subo3cubename, objnames[i], redshift[i], filenum[i]
;CUBE_EXTRACT, subo3cubename, objnames[i] + '_o3', xcen, ycen, radius, filenum[i]
;SNIFS_CUBE_FITTING, subo3cubename, objnames[i], parinfo, redshift[i], filenum[i]

SNIFS_MKEELRIMG, objnames[i], filenum[i], redshift[i], cubenames[i], mincluster=12
;SNIFS_JPG, objnames[i], filenum[i], redshift[i]


ENDFOR

END




PRO get_eelrimg

;Only get those with confirmed EELR
READCOL, 'EELR_namefile', objnames, filenum, redshift, linebeg, lineend, nocontsub, FORMAT='A,I,F,F,F,I'
ncubes = N_ELEMENTS(objnames)

FOR i=0,ncubes-1 DO BEGIN
;FOR i=26,26 DO BEGIN
	SNIFS_MKEELRIMG, objnames[i], filenum[i], redshift[i], o3sumbeg=linebeg[i], o3sumend=lineend[i], nocontsub=nocontsub[i], mincluster=5
	;SNIFS_JPG, objnames[i], filenum[i], redshift[i]
ENDFOR

END









PRO snifs_mkeelrimg, objname, filenum, redshift, cubename, o3sumbeg=o3sumbeg, o3sumend=o3sumend, nocontsub=nocontsub, mincluster = mincluster

;mincluster = Mininum number of points in a cluster to count as significant

RESTORE, 'snifs_objinfo.save'
structnum = WHERE((objinfo.Name EQ objname) AND (objinfo.Filenum EQ filenum))

;Scale and subtract continuum contribution to EELR

contname = objname + 'contimg' + STRTRIM(filenum,2) + '.fits'
subo3cubename = objname + "subo3cube_" + STRTRIM(filenum,2) + ".fits"
eelrimg_name = objname + 'eelrimg' + STRTRIM(filenum,2) + '.fits'
coneelrimg_name = objname + 'con_eelrimg' + STRTRIM(filenum,2) + '.fits'

cont = READFITS(contname, /SILENT)
subo3cube = READFITS(subo3cubename, slicehead, /SILENT)



;Find the slices containing the bulk of the OIII emission for summing

IF ((KEYWORD_SET(o3sumbeg)) AND (KEYWORD_SET(o3sumend))) THEN BEGIN
nslice = N_ELEMENTS(subo3cube[0,0,*])
crval3 = SXPAR(slicehead, 'CRVAL3')
cd3_3 = SXPAR(slicehead, 'CD3_3')
tmp = FINDGEN(nslice)
tmp1 = crval3 + (tmp*cd3_3) - o3sumbeg
tmp2 = crval3 + (tmp*cd3_3) - o3sumend
tmpmin1 = MIN(ABS(tmp1), o3sumbegidx)
tmpmin2 = MIN(ABS(tmp2), o3sumendidx)
slice = TOTAL(subo3cube[*,*,o3sumbegidx: o3sumendidx], 3)
ENDIF ELSE BEGIN
slice = TOTAL(subo3cube, 3)
ENDELSE


;Find centroid of continuum image
CNTRD, cont, 7, 7, xcentr, ycentr, 2
IF ((xcentr GT 0.0) AND (ycentr GT 0.0)) THEN BEGIN
	ctr_coord=[xcentr, ycentr]
ENDIF ELSE BEGIN
	ctr_coord=[7,7]
ENDELSE

;Fit and subtract continuum such that the residual in 'circle' is minimized
xdim = N_ELEMENTS(cont[*,1])
ydim = N_ELEMENTS(cont[1,*])

DIST_CIRCLE, circles, [xdim,ydim], ctr_coord[0], ctr_coord[1]
circleidx = WHERE(circles LE 2.0)
maskidx = WHERE(circles LE 3.0)
mask = FLTARR(xdim,ydim)
mask[maskidx] = 1.0


guessratio = median(slice[circleidx]) / median(cont[circleidx])
IF (guessratio LT 0.0) THEN guessratio = 1.0

functargs = {CIRCLEIDX:circleidx, EELR:slice, CONT:cont}
bestratio = MPFIT('snifs_scale_cont', guessratio, FUNCTARGS=functargs, /QUIET)
IF (bestratio LT 0.0) THEN bestratio = 1.0

;Mask continuum image so that only central circle of 4 pixel radius is used for subtraction
cont = cont * mask
IF (NOT KEYWORD_SET(nocontsub)) THEN eelrimg = slice - (cont * bestratio[0])
IF (KEYWORD_SET(nocontsub)) THEN eelrimg = slice
;WRITEFITS, eelrimg_name, eelrimg, slicehead



;Convert EELR image to absolute luminosity using luminosity distance calculated from redshift
convert_eelr = FLUX_CONVERT_SLICE(eelrimg, redshift)
convert_slice = FLUX_CONVERT_SLICE(slice, redshift)
;WRITEFITS, coneelrimg_name, convert_eelr, slicehead

;Get nuclear OIII from the inner 2 pixel radius
nuclearo3 = TOTAL(convert_slice[circleidx])

goodidx = WHERE((circles GT 2.0) AND (convert_eelr GT 0.0) AND (convert_eelr LT 5.e+41))

;Make a cut and only keep pixels that have values above the mean, but below (9*sigma + mean)
good_eelr = convert_eelr[goodidx]
MEANCLIP, good_eelr, eelrmean, eelrsig, CLIPSIG=6, MAXITER=20, CONVERGE_NUM=0.2, SUBS=meanidx
eelridx = WHERE((convert_eelr GT (eelrmean + 0.*eelrsig)) AND (convert_eelr LT (eelrmean + 9*eelrsig)))


;Use hierarchical clustering to find the clumps of the high value data points
nrows = 15
eelrx = (eelridx mod nrows) 
eelry = FIX(eelridx / nrows)

nidx = N_ELEMENTS(eelridx)
carray = FLTARR(2,nidx)
carray[0,*] = eelrx
carray[1,*] = eelry

distance = DISTANCE_MEASURE(carray)
clusters = CLUSTER_TREE(distance, linkdistance)
;PRINT, [clusters, TRANSPOSE(linkdistance)], FORMAT='(I3, I7, F10.2)'
;WINDOW, XSIZE=1200, YSIZE=900
;!p.multi=[0,1,2]
;PLOT, eelrx, eelry, psym=1, xrange=[0,15], yrange=[0,15], POSITION=[0.30, 0.55, 0.70, 0.95]
;DENDRO_PLOT, clusters, linkdistance
;imgname =  objname + 'cluster'+ STRTRIM(filenum,2)
;imgval = TVREAD(filename=imgname,/jpeg,/nodialog,quality=100)

;;;Get date and number of rawfile
;;;datenum = DATENUM_EXTRACT(cubename)

clusterlevels = DIFFERENT(linkdistance)
IF (N_ELEMENTS(clusterlevels) LT 2) THEN BEGIN

	totalo3 = TOTAL(convert_eelr)

ENDIF ELSE BEGIN

	result = CUT_TREE(clusters, linkdistance, 1)
	cluster = N_ELEMENTS(result.group)
	clustersizes = result.length
	goodcluster = WHERE(clustersizes GE mincluster)
	ngoodcluster = N_ELEMENTS(goodcluster)

	IF (goodcluster[0] GE 0) THEN BEGIN
		FOR i=0,ngoodcluster-1 DO BEGIN
		tmpelements = result[goodcluster[i]].elements
		goodidx = WHERE(tmpelements GE 0)
		tmpelements = tmpelements[goodidx]

		IF (i EQ 0) THEN goodelement = tmpelements
		IF (i GT 0) THEN goodelement = [goodelement ,tmpelements]

		ENDFOR
;WINDOW, XSIZE=600, YSIZE=600
;!p.multi=[0,1,1]
;plot, eelrx[goodelement], eelry[goodelement], psym=1, xrange=[0,15], yrange=[0,15], POSITION=[0.10, 0.15, 0.90, 0.95]
;imgname =  objname + 'summed_points'+ STRTRIM(filenum,2)
;imgval = TVREAD(filename=imgname,/jpeg,/nodialog,quality=100)
		goodx = eelrx[goodelement]
		goody = eelry[goodelement]
		totalo3 = TOTAL(convert_eelr[goodx, goody])

	ENDIF ELSE BEGIN
		totalo3 = 0.0
	ENDELSE
ENDELSE

IF (structnum GE 0) THEN BEGIN
	objinfo[structnum].O3nuclear = nuclearo3
	objinfo[structnum].O3extended = totalo3
ENDIF

SAVE, objinfo, filename='snifs_objinfo.save'

END









PRO print_struct

RESTORE, 'snifs_objinfo.save'
nstruct = N_ELEMENTS(objinfo.name)


FOR i=0,nstruct-1 DO BEGIN

	PRINT, objinfo[i].name, '	', objinfo[i].O3nuclear, '	', objinfo[i].O3extended, '	', objinfo[i].gooddata

ENDFOR


END









FUNCTION snifs_scale_cont, ratio, CIRCLEIDX=circleidx, EELR=eelr, CONT=cont

scalefactor = 1.e20
deviates = (eelr[circleidx] - (cont[circleidx] * ratio[0])) * scalefactor
RETURN, deviates

END






PRO snifs_jpg, objname, filenum, redshift


Ho=71
Omega_m = 0.27
Omega_vac = 0.73

result = NWCC(redshift, Ho=Ho, Omega_m=Omega_m, Omega_vac=Omega_vac)
jpgname = objname + 'display'+ STRTRIM(filenum,2)
kpc_per_arcsec = result[1]

!p.multi=[0,3,1]
window, RETAIN=2, XSIZE=900, YSIZE=300		;Set RETAIN=2 prevents error message "% X windows protocol error: BadMatch (invalid parameter attributes)."

;Make jpg image of O III slice

o3slice_name = objname + 'o3slice' + STRTRIM(filenum,2) + '.fits'
o3slice = READFITS(o3slice_name, /SILENT)

jpgimage = o3slice/1.d40
scmax = MAX(jpgimage) * 0.3
scmin = MIN(jpgimage) * 0.2
bar_interval = FIX(scmax - scmin)/2

MKPLOT2, label=objname + '[O III] image', image= jpgimage, newfilename='tmpjpg', rotnum=0, rotang=0, scmin=scmin, scmax=scmax, center=[8,8], bar_interval= bar_interval, unit='E40', /maxred, /kpc_scale, kpc_per_arcsec=kpc_per_arcsec, /indcolorbar



;Make jpg image of EELR image

eelr_name = objname + 'eelrimg' + STRTRIM(filenum,2) + '.fits'
eelr = READFITS(eelr_name, /SILENT)

jpgimage = eelr
scmax = MAX(jpgimage) * 0.6
scmin = MIN(jpgimage) * 0.8
bar_interval = FIX(scmax - scmin)/2

MKPLOT2, label=objname + 'EELR image', image= jpgimage, newfilename='tmpjpg', rotnum=0, rotang=0, scmin=scmin, scmax=scmax, center=[8,8], bar_interval= bar_interval, unit=' ', /maxred, /kpc_scale, kpc_per_arcsec=kpc_per_arcsec, /indcolorbar




;Make jpg image of velocity slice

vel_slicename1 = objname + 'vel_slice1_' + STRTRIM(filenum,2) + '.fits'
velslice = READFITS(vel_slicename1, /SILENT)

jpgimage = velslice
medvel = MEDIAN(jpgimage[4:10,4:10])
stddev_vel = STDDEV(jpgimage[4:10,4:10])
scmax = medvel + (2.0 * stddev_vel)
scmin = medvel - (2.0 * stddev_vel)
IF (scmax GT 500.0) THEN scmax = 500.0
IF (scmin LT -500.0) THEN scmin = -500.0
IF (scmax LT scmin) THEN scmax = scmin+1000.0
bar_interval = FIX(scmax - scmin)/2
;print, bar_interval , scmin, scmax, medvel, stddev_vel

MKPLOT2, label=objname + 'Velocity image', image= jpgimage, newfilename=jpgname, rotnum=0, rotang=0, scmin=scmin, scmax=scmax, center=[8,8], bar_interval= bar_interval, unit='km/s', /maxred, /kpc_scale, kpc_per_arcsec=kpc_per_arcsec, /indcolorbar



END









PRO get_namefile

;Temporary fix.....

READCOL, 'cubenames', cubenames, FORMAT='A'
READCOL, 'rawfilenames', rawfilenames, FORMAT='A'
READCOL, 'redshift_flux', objname1, redshift, flux, FORMAT='A,F,F'

ncubes = N_ELEMENTS(cubenames)
objlist = STRARR(ncubes)

k=0 ;Counter for objname comparison

OPENW, 1, 'namefile_tmp', WIDTH=1000

FOR i=0,ncubes-1 DO BEGIN

objname = STRTRIM(NAME_EXTRACT(rawfilenames[i]),2)
obj=objname


FOR j=1,10 DO BEGIN
	IF (DUPLICATE(objlist, obj) EQ 0) THEN BREAK
	IF (DUPLICATE(objlist, obj) EQ 1) THEN obj = objname + STRTRIM(j,2)
ENDFOR

objlist[i] = obj
filenum=j-1

IF (objname1[k] EQ objname) THEN BEGIN
	PRINTF, 1, objname, '	', redshift[k], '	', cubenames[i], '	', flux[k], '	', filenum
	k += 1
ENDIF

ENDFOR

CLOSE, 1

END






FUNCTION cut_tree, clusters, linkdistance, cutlevel

nc = N_ELEMENTS(linkdistance)
mcnum = INDGEN(nc)

rcluster1 = clusters[0,*]
rcluster2 = clusters[1,*]
mergenums = mcnum + nc + 1
levels = linkdistance

counted = FLTARR(nc)
clusterlevels = DIFFERENT(linkdistance)
totalmerges = N_ELEMENTS(WHERE(linkdistance LE clusterlevels[cutlevel]))

WHILE (TOTAL(counted) LT totalmerges) DO BEGIN
	
	notcounted = WHERE(counted EQ 0)
	startidx = notcounted[0]
	group = [rcluster1[startidx], rcluster2[startidx]]
	level = levels[startidx]
	mergenum = mergenums[startidx]
	counted[startidx] = 1.

	WHILE (level LT clusterlevels[cutlevel]) DO BEGIN

		result = FINDNEXTCLUS(rcluster1, rcluster2, mergenum)
		nidx = result[0]
		nclus = result[1]
		level = levels[nidx]
		mergenum = mergenums[nidx]
		IF((nclus LE nc) AND (level LT clusterlevels[cutlevel]) AND (counted[nidx] NE 1)) THEN group = [group, nclus]
		
		;If we encounter a cluster number larger than the number of data points we have, back track recursively to find all the elements merged into that cluster
		IF ((nclus GT nc) AND (level LT clusterlevels[cutlevel])) THEN BEGIN
			group = BACKTRACK(rcluster1, rcluster2, mergenums, nclus, group, counted, countedidx=countedidx)
			counted[countedidx]=1.
		ENDIF
	counted[nidx] = 1.
	ENDWHILE

	IF(startidx EQ 0) THEN tmpleaves = [-3, mergenum, group]
	IF(startidx GT 0) THEN tmpleaves = [tmpleaves,-3,mergenum,group]	
ENDWHILE


;Return an structure with merged group number, number of elements, and group member elements
tmpleaves = [tmpleaves,-3]
delim = WHERE(tmpleaves EQ -3)
nleaves = N_ELEMENTS(delim) - 1
lengths = delim[1:nleaves] - delim[0:nleaves-1] - 2
maxlength = MAX(lengths)

strleaf = {group:0., length:0., elements:fltarr(maxlength)-1.}
leaves = REPLICATE(strleaf, nleaves)

FOR i = 0, nleaves-1 DO BEGIN

groupname = tmpleaves[delim[i]+1]
groupmembers = tmpleaves[delim[i]+2:delim[i+1]-1]

leaves[i].group = groupname
leaves[i].length = lengths[i]
leaves[i].elements[0:lengths[i]-1] = groupmembers

ENDFOR

RETURN, leaves

END



FUNCTION findrawclus, rclus1, rlcus2, mergenums, mergenum
;Return the two pre-merged clusters given a merged cluster number
	mergeidx = WHERE(mergenums EQ mergenum)
	rawclus1 = rclus1[mergeidx]
	rawclus2 = rlcus2[mergeidx]

RETURN, [rawclus1, rawclus2, mergeidx]
END



FUNCTION findnextclus, rclus1, rlcus2, mergenum

	next1 = WHERE(rclus1 EQ mergenum)
	next2 = WHERE(rlcus2 EQ mergenum)

	IF (next1 GE 0) THEN BEGIN
		nidx = next1[0]
		nclus = rlcus2[nidx]
	ENDIF

	IF (next2 GE 0) THEN BEGIN
		nidx = next2[0]
		nclus = rclus1[nidx]
	ENDIF	
RETURN, [nidx, nclus]
END



FUNCTION backtrack, rcluster1, rcluster2, mergenums, nclus, group, counted, countedidx=countedidx

FORWARD_FUNCTION backtrack
countedidx = [0]
a = 0
nc = N_ELEMENTS(rcluster1)
WHILE (a eq 0) DO BEGIN
	result = FINDRAWCLUS(rcluster1, rcluster2, mergenums, nclus)
	rawclusters = result[0:1]
	clusteridx = result[2] 
	countedidx = [countedidx,clusteridx]
								
	IF ((rawclusters[0] LE nc) AND (counted[clusteridx] NE 1)) THEN group = [group, rawclusters[0]]
	IF ((rawclusters[1] LE nc) AND (counted[clusteridx] NE 1)) THEN group = [group, rawclusters[1]]
	IF (rawclusters[0] GT nc) THEN nclus = rawclusters[0]
	IF (rawclusters[1] GT nc) THEN nclus = rawclusters[1]
	IF ((rawclusters[1] LE nc) AND (rawclusters[0] LE nc)) THEN a = 1
	IF ((rawclusters[1] GT nc) AND (rawclusters[0] GT nc)) THEN group = BACKTRACK(rcluster1, rcluster2, mergenums, rawclusters[0], group, counted, countedidx=countedidx)
	counted[clusteridx] = 1.	
ENDWHILE


RETURN, group


END




FUNCTION different, x

;Extract the distinct elements from input

diff = x[0]
comp = diff
j = 0
n = N_ELEMENTS(x)

	FOR i=0,n-1 DO BEGIN
		IF (x[i] NE comp) THEN BEGIN
			diff = [diff,x[i]]
			comp = x[i]
		ENDIF
	ENDFOR

RETURN, diff

END



















PRO get_line_slice

;READCOL, 'rawfilenames', rawfilenames, FORMAT='A'

;ncubes = N_ELEMENTS(rawfilenames)
;objnames=STRARR(ncubes)

;OPENW, 1, 'objnames'

;FOR i=0,ncubes-1 DO BEGIN

;objname = NAME_EXTRACT(rawfilenames[i])
;objname = STRTRIM(objname,2)
;PRINTF, 1, objname

;subskyname = objname + '_subskycube.fits'

;IF (FILE_TEST(subskyname) EQ 1) THEN BEGIN
;
;	subskyname = objname + "_subskycube_1.fits"
;	IF (FILE_TEST(subskyname) EQ 1) THEN subskyname = objname + "_subskycube_2.fits"
;
;ENDIF


;crval3 = 5101.13
;cd3_3 = 2.93

;cubename = 'HB0110+297_subskycube_1.fits'
;line_coord = [584,594]
;outfits = 'HB0110+297_O3_1.fits'

;cubename = '4C_subskycube_1.fits'
;line_coord = [770,790]
;outfits = '4C+38.35_O3_1.fits'

;cubename = '3C284_subskycube_2.fits'
;line_coord = [370,385]
;outfits = '3C284_O3_2.fits'

;cubename = '4C+34.42_subskycube.fits'
;line_coord = [650,660]
;outfits = '4C+34.42_O3.fits'

;cubename = '4C+47.37_subskycube.fits'
;line_coord = [610,620]
;outfits = '4C+47.37_O3.fits'

;cubename = '4C+00.46_subskycube.fits'
;line_coord = [678,693]
;outfits = '4C+00.46_O3.fits'

cubename = '3C166_subskycube.fits'
line_coord = [370,390]
outfits = '3C166_O3.fits'

cube = READFITS(cubename, cubehead)
line_slice = LINE_SLICE(cube, line_coord)
convert_slice = FLUX_CONVERT_SLICE(line_slice, 1925.7)
WRITEFITS, outfits, convert_slice, cubehead


;ENDFOR

CLOSE, 1
END







