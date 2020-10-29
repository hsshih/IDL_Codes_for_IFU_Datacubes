PRO parinfostruct__define

tmp = {parinfostruct, $
		value:0.D, $	
		fixed:0, $	
		limited:[0,0], $	
		limits:[0.D,0.D] $		
			}
END




PRO mko3sub, filename, objname, redshift, b_begin, b_end, r_begin, r_end, o3_begin, o3_end
;Get Median continuum
;; ***Make sure CRVAL3 and CD3_3 are correct in the header, if 'emptycube.fits' is copied, header need to be edited
;; 01/25/13 Modified the way veldisp is calculated such that instead of being divided by redshifted 5006.9, it is divided by the fit centroid of that line, but difference is too small to tell

data = READFITS(filename, head, /SILENT)
bluecont  = data[*,*, b_begin:b_end]
redcont = data[*,*, r_begin:r_end]

nb = N_ELEMENTS(bluecont[1,1,*])
nr = N_ELEMENTS(redcont[1,1,*])

;For small cubes with 0.2" sampling
cont = FLTARR(32,50,nb+nr)
xc = 16
yc = 25

cont(*,*,0:nb-1) = bluecont
cont(*,*,nb:nb+nr-1) = redcont

contmed = MEDIAN(cont, DIMENSION=3)


;Subtract Continuum

crval3 = SXPAR(head, 'CRVAL3')
cd3_3 = SXPAR(head, 'CD3_3')
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

;Trim off useless parts of the combined cube

;For small cubes with 0.2" sampling
xstart = 6
xend = 25
ystart = 13
yend = 37


;For big cubes with 0.1" sampling
;xstart = 12
;xend = 55
;ystart = 25
;yend = 75

zstart = 0
zend = N_ELEMENTS(subcube[1,1,*]) - 1

trimcube = subcube[xstart:xend, ystart:yend, zstart:zend]

fitsname = objname + 'sub_o3cube.fits'
trimcubename = objname + 'trimcube.fits'

;WRITEFITS, 'tmpcontmed.fits', contmed
WRITEFITS, fitsname, subcube, head
WRITEFITS, trimcubename, trimcube, head


END





PRO makecubeslice, fitsname, objname, minsn

;Read in the continuum subtracted [O III] cube to determine S/N
cube = READFITS(fitsname, cubehead, /SILENT)

cubex = N_ELEMENTS(cube[*,1,1])
cubey = N_ELEMENTS(cube[1,*,1])
cubez = N_ELEMENTS(cube[1,1,*])

medcube = FLTARR(cubex,cubey,cubez)

xarr = FLTARR(cubex, cubey)
yarr = FLTARR(cubex, cubey)

FOR i=0,cubex-1 DO BEGIN
FOR j=0,cubey-1 DO BEGIN

	spectrum = cube[i,j,*]
	stddev = STDDEV(spectrum)

FOR k=0,cubez-1 DO BEGIN

	val = cube[i,j,k]

	nmed = 20
	nmed_half = nmed / 2
	std = 50

	IF ((k GT nmed_half) AND (k LT (cubez-nmed_half))) THEN medval = MEDIAN(spectrum[k-nmed_half:k+nmed_half])
	IF (k LT nmed_half) THEN medval = MEDIAN(spectrum[k:k+nmed])
	IF (k GT (cubez-nmed_half)) THEN medval = MEDIAN(spectrum[k-nmed:k])

	;IF ((val GT (medval + std * stddev)) OR (val LT (medval + std * stddev))) THEN val = medval

	
	medcube[i,j,k] = medval
	cube[i,j,k] = val

ENDFOR
	xarr[i,j]=i
	yarr[i,j]=j
ENDFOR
ENDFOR


cubeslice = TOTAL(cube, 3)
negidx = WHERE(cubeslice LT 1) 
cubeslice[negidx] = 1

noiseslice = SQRT(cubeslice)

;Eliminate any pixels that does not reach the minimum signal to noise
sigtonoise = cubeslice / noiseslice
goodsigidx = WHERE(sigtonoise GT minsn)	

cubeslicename = 'cubeslice' + objname + '.fits'
noiseslicename = 'noiseslice' + objname + '.fits'
sigtonoisename = 'StoNslice' + objname + '.fits'
medcubename = 'medcube' + objname + '.fits'
sub_o3txtname = objname + 'sub_o3.txt'


WRITEFITS, cubeslicename, cubeslice[xarr[goodsigidx], yarr[goodsigidx]], cubehead
WRITEFITS, noiseslicename, noiseslice[xarr[goodsigidx], yarr[goodsigidx]], cubehead
WRITEFITS, sigtonoisename, sigtonoise, cubehead

xarr = FIX(xarr + 1)
yarr = FIX(yarr + 1)

WRITEFITS, medcubename, medcube, cubehead
WRITECOL, sub_o3txtname, xarr[goodsigidx], yarr[goodsigidx], cubeslice[goodsigidx], noiseslice[goodsigidx]


END





pro dovoronoi, objname, targetSN

sub_o3txtname = objname + 'sub_o3.txt'
rdfloat, sub_o3txtname, x, y, signal, noise

; Load a colortable and open a graphic window
;
loadct, 5
r = GET_SCREEN_SIZE()
window, xsize=r[0]*0.4, ysize=r[1]*0.8

; Perform the actual computation. The vectors
; (binNum, xNode, yNode, xBar, yBar, sn, nPixels, scale)
; are all generated in *output*
;
voronoi_2d_binning, x, y, signal, noise, targetSN, $
    binNum, xNode, yNode, xBar, yBar, sn, nPixels, scale, /PLOT, /QUIET

; Save to a text file the initial coordinates of each pixel together
; with the corresponding bin number computed by this procedure
;

outname1 = objname + 'sub_o3_output.txt'
outname2 = objname + 'sub_o3_output2.txt'

nbins = N_ELEMENTS(xNode)
;bins_ordered = INDGEN(nbins)

astrolib
forprint, x, y, binNum, TEXTOUT=outname1, $
    COMMENT='         X"          Y"      BIN_NUM'
forprint, xNode, yNode, sn, nPixels, TEXTOUT=outname2, $
    COMMENT='         XNODE          YNODE 	SN		NPIX'

END





PRO cube_fitting, fitsname, errcubename, objname, parinfo, redshift, compnum, fix_ctr_coord=fix_ctr_coord

IF (not KEYWORD_SET(fix_ctr_coord)) THEN fix_ctr_coord=[10,12]

npar = compnum * 3 + 1

cube = READFITS(fitsname, cubehead, /SILENT)
errcube = READFITS(errcubename, /SILENT)

outname1 = objname + 'sub_o3_output.txt'
outname2 = objname + 'sub_o3_output2.txt'

READCOL, outname1, x, y, binnum, FORMAT='F,F,I'
READCOL, outname2, xnode, ynode, sn, npix, FORMAT='F,F,F,I'

nbin = N_ELEMENTS(xnode)

cubex = N_ELEMENTS(cube[*,1,1])
cubey = N_ELEMENTS(cube[1,*,1])
cubez = N_ELEMENTS(cube[1,1,*])

wavelength = FINDGEN(cubez)

fit_result = FLTARR(npar,nbin)
fit_error = FLTARR(npar,nbin)

vel_slice1 = FLTARR(cubex,cubey) - 10000.0
veldisp_slice1 = FLTARR(cubex,cubey)
flux_slice1 = FLTARR(cubex,cubey)

IF ((compnum EQ 2) OR (compnum EQ 3) OR (compnum EQ 4)) THEN BEGIN

vel_slice2 = FLTARR(cubex,cubey) - 10000.0
veldisp_slice2 = FLTARR(cubex,cubey)
flux_slice2 = FLTARR(cubex,cubey)

ENDIF

IF ((compnum EQ 3) OR (compnum EQ 4)) THEN BEGIN

vel_slice3 = FLTARR(cubex,cubey) - 10000.0
veldisp_slice3 = FLTARR(cubex,cubey)
flux_slice3 = FLTARR(cubex,cubey)

ENDIF

IF (compnum EQ 4) THEN BEGIN

vel_slice4 = FLTARR(cubex,cubey) - 10000.0
veldisp_slice4 = FLTARR(cubex,cubey)
flux_slice4 = FLTARR(cubex,cubey)

ENDIF



crval3 = SXPAR(cubehead, 'CRVAL3')
cd3_3 = SXPAR(cubehead, 'CD3_3')
o3_redshifted = 5006.9 * (1 + redshift)

IF (objname EQ '3C124') THEN o3_redshifted = 3727.0 * (1 + redshift)

txtname = objname + 'result.txt'
OPENW, 1, txtname

;Scale everything by a large constant, MPFIT can't deal with small values for some reason...
constscale = 1.e20
parinfo[3].value = parinfo[3].value * constscale
IF ((compnum EQ 2) OR (compnum EQ 3) OR (compnum EQ 4)) THEN parinfo[6].value = parinfo[6].value * constscale
IF ((compnum EQ 3) OR (compnum EQ 4)) THEN parinfo[9].value = parinfo[9].value * constscale
IF (compnum EQ 4) THEN parinfo[12].value = parinfo[12].value * constscale

;Extract the initial original initial fit parameter
iniparm = parinfo[*].value

FOR i=0,nbin-1 DO BEGIN

	binidx = WHERE(binnum EQ i)
	subcube = cube[x[binidx]-1,y[binidx]-1,*]
	tmpplane = TOTAL(subcube, 1)
	spectrum = TOTAL(tmpplane,1)

	;Find neighbor fit result and make it the initial fit parameter
	;For pixels closer to the central pixel, use original initial fit parameter
	ctrxs = x[binidx]
	ctrys = y[binidx]
	fix_xctr = fix_ctr_coord[0]
	fix_yctr = fix_ctr_coord[1]
	dist_from_center = SQRT((ctrxs[0]-fix_xctr)^2. + (ctrys[0]-fix_yctr)^2.)
	
	
	neighboridx=WHERE((ABS(xnode-ctrxs[0]) LE 1.0) AND (ABS(ynode-ctrys[0]) LE 1.0))
;	dist_from_pixel = SQRT((xnode-ctrxs[0])^2. + (ynode-ctrys[0])^2.)
;	neighboridx=WHERE(dist_from_pixel LE 2.0)
	IF ((neighboridx[0] NE -1) AND (dist_from_center GE 2.0)) THEN BEGIN

		neighboridx_nonz= WHERE(fit_result[1,neighboridx] NE 0.0)

		IF (neighboridx_nonz[0] NE -1) THEN BEGIN

			neighboridx = neighboridx[neighboridx_nonz]
			IF (N_ELEMENTS(neighboridx) GT 1) THEN neighbor_parm = TOTAL(fit_result[*, neighboridx], 2) / N_ELEMENTS(neighboridx)
			IF (N_ELEMENTS(neighboridx) EQ 1) THEN neighbor_parm = fit_result[*, neighboridx]
			;PRINT, ctrxs[0], ctrys[0], neighbor_parm
			parinfo[*].value=neighbor_parm

		ENDIF

	ENDIF ELSE BEGIN
		parinfo[*].value=iniparm
		;PRINT, "No neighbor parm"
	ENDELSE

	subcube_err = errcube[x[binidx]-1,y[binidx]-1,*]
	tmpplane_err = TOTAL(subcube_err, 1)
	spectrum_err = TOTAL(tmpplane_err,1)

	spectrum_scaled = spectrum * constscale
	spectrum_err_scaled = spectrum_err * constscale
	

	modelname = objname + 'modspec' + STRTRIM(FIX(xnode[i]),2) +'_' + STRTRIM(FIX(ynode[i]),2)
	modelresult = MYFIT(wavelength, spectrum_scaled, spectrum_err_scaled, modelname, parinfo, compnum)
	IF (compnum EQ 1) THEN model = modelresult(0)+gauss1(wavelength, modelresult(1:3))
	IF (compnum EQ 2) THEN model = modelresult(0)+gauss1(wavelength, modelresult(1:3))+gauss1(wavelength, modelresult(4:6))
	IF (compnum EQ 3) THEN model = modelresult(0)+gauss1(wavelength, modelresult(1:3))+gauss1(wavelength, modelresult(4:6))+gauss1(wavelength, modelresult(7:9))
	IF (compnum EQ 4) THEN model = modelresult(0)+gauss1(wavelength, modelresult(1:3))+gauss1(wavelength, modelresult(4:6))+gauss1(wavelength, modelresult(7:9))+gauss1(wavelength, modelresult(10:12))

	model_err = MYMC(wavelength, model, spectrum_err, parinfo, compnum)

	fit_result[*,i] = modelresult
	fit_error[*,i] = model_err

	;Using sigma for velocity dispersion sigma = FWHM/2.355
	;Instrumentation Broadening from arc images (There were slight deviations between two gratings 3.74 vs 3.96, but close enough....)
	sig_inst = 1.64

	vel1 = ((((cd3_3 * modelresult[1]) + crval3) - o3_redshifted) / o3_redshifted) * 3.E5
	sigma1 = SQRT((modelresult[2]^2.) - (sig_inst^2.))
	IF (NOT (sigma1 GE 0.0)) THEN sigma1 = 0.0
	;;;;veldisp1 = ((sigma1 * cd3_3) / o3_redshifted) * 3.E5
	veldisp1 = ((sigma1 * cd3_3) / ((cd3_3 * modelresult[1]) + crval3)) * 3.E5

	vel_slice1[x[binidx]-1,y[binidx]-1] = vel1
	veldisp_slice1[x[binidx]-1,y[binidx]-1] = veldisp1	;Subtract instrumentation broadening in quadrature
	flux_slice1[x[binidx]-1,y[binidx]-1] = modelresult[3] / npix[i] / constscale

	IF (compnum EQ 1) THEN PRINTF, 1, xnode[i], ynode[i], fit_result[*,i], FORMAT='(F9.2,F9.2,F9.2,F9.2,F9.2,F11.1)' 

	IF (compnum EQ 2) THEN BEGIN

		vel2 = ((((cd3_3 * modelresult[4]) + crval3) - o3_redshifted) / o3_redshifted) * 3.E5
		sigma2 = SQRT((modelresult[5]^2.) - (sig_inst^2.))
		IF (NOT (sigma2 GE 0.0)) THEN sigma2 = 0.0
		;;;;veldisp2 = ((sigma2 * cd3_3) / o3_redshifted) * 3.E5
		veldisp2 = ((sigma2 * cd3_3) / ((cd3_3 * modelresult[4]) + crval3)) * 3.E5
		vel_slice2[x[binidx]-1,y[binidx]-1] = vel2
		veldisp_slice2[x[binidx]-1,y[binidx]-1] = veldisp2
		flux_slice2[x[binidx]-1,y[binidx]-1] = modelresult[6] / npix[i] / constscale

		PRINTF, 1, xnode[i], ynode[i], fit_result[*,i], FORMAT='(F9.2,F9.2,F9.2,F9.2,F9.2,F13.1,F9.2,F9.2,F13.1)' 

	ENDIF 

	IF (compnum EQ 3) THEN BEGIN

		vel2 = ((((cd3_3 * modelresult[4]) + crval3) - o3_redshifted) / o3_redshifted) * 3.E5
		sigma2 = SQRT((modelresult[5]^2.) - (sig_inst^2.))
		IF (NOT (sigma2 GE 0.0)) THEN sigma2 = 0.0
		;;;;veldisp2 = ((sigma2 * cd3_3) / o3_redshifted) * 3.E5
		veldisp2 = ((sigma2 * cd3_3) / ((cd3_3 * modelresult[4]) + crval3)) * 3.E5
		vel_slice2[x[binidx]-1,y[binidx]-1] = vel2
		veldisp_slice2[x[binidx]-1,y[binidx]-1] = veldisp2
		flux_slice2[x[binidx]-1,y[binidx]-1] = modelresult[6] / npix[i] / constscale

		vel3 = ((((cd3_3 * modelresult[7]) + crval3) - o3_redshifted) / o3_redshifted) * 3.E5
		sigma3 = SQRT((modelresult[8]^2.) - (sig_inst^2.))
		IF (NOT (sigma3 GE 0.0)) THEN sigma3 = 0.0
		;;;;veldisp3 = ((sigma3 * cd3_3) / o3_redshifted) * 3.E5
		veldisp3 = ((sigma3 * cd3_3) / ((cd3_3 * modelresult[7]) + crval3)) * 3.E5
		vel_slice3[x[binidx]-1,y[binidx]-1] = vel3
		veldisp_slice3[x[binidx]-1,y[binidx]-1] = veldisp3
		flux_slice3[x[binidx]-1,y[binidx]-1] = modelresult[9] / npix[i] / constscale


		PRINTF, 1, xnode[i], ynode[i], fit_result[*,i], FORMAT='(F7.1,F7.1,F7.1,F7.1,F7.1,F9.1,F7.1,F7.1,F9.1,F7.1,F7.1,F9.1)' 

	ENDIF 

	IF (compnum EQ 4) THEN BEGIN

		vel2 = ((((cd3_3 * modelresult[4]) + crval3) - o3_redshifted) / o3_redshifted) * 3.E5
		sigma2 = SQRT((modelresult[5]^2.) - (sig_inst^2.))
		IF (NOT (sigma2 GE 0.0)) THEN sigma2 = 0.0
		;;;;veldisp2 = ((sigma2 * cd3_3) / o3_redshifted) * 3.E5
		veldisp2 = ((sigma2 * cd3_3) / ((cd3_3 * modelresult[4]) + crval3)) * 3.E5
		vel_slice2[x[binidx]-1,y[binidx]-1] = vel2
		veldisp_slice2[x[binidx]-1,y[binidx]-1] = veldisp2
		flux_slice2[x[binidx]-1,y[binidx]-1] = modelresult[6] / npix[i] / constscale

		vel3 = ((((cd3_3 * modelresult[7]) + crval3) - o3_redshifted) / o3_redshifted) * 3.E5
		sigma3 = SQRT((modelresult[8]^2.) - (sig_inst^2.))
		IF (NOT (sigma3 GE 0.0)) THEN sigma3 = 0.0
		;;;;veldisp3 = ((sigma3 * cd3_3) / o3_redshifted) * 3.E5
		veldisp3 = ((sigma3 * cd3_3) / ((cd3_3 * modelresult[7]) + crval3)) * 3.E5
		vel_slice3[x[binidx]-1,y[binidx]-1] = vel3
		veldisp_slice3[x[binidx]-1,y[binidx]-1] = veldisp3
		flux_slice3[x[binidx]-1,y[binidx]-1] = modelresult[9] / npix[i] / constscale

		vel4 = ((((cd3_3 * modelresult[10]) + crval3) - o3_redshifted) / o3_redshifted) * 3.E5
		sigma4 = SQRT((modelresult[11]^2.) - (sig_inst^2.))
		IF (NOT (sigma4 GE 0.0)) THEN sigma4 = 0.0
		;;;;veldisp4 = ((sigma4 * cd3_3) / o3_redshifted) * 3.E5
		veldisp4 = ((sigma4 * cd3_3) / ((cd3_3 * modelresult[10]) + crval3)) * 3.E5
		vel_slice4[x[binidx]-1,y[binidx]-1] = vel4
		veldisp_slice4[x[binidx]-1,y[binidx]-1] = veldisp4
		flux_slice4[x[binidx]-1,y[binidx]-1] = modelresult[12] / npix[i] / constscale

;print, modelresult
		PRINTF, 1, xnode[i], ynode[i], fit_result[*,i], FORMAT='(F7.1,F7.1,F7.1,F7.1,F7.1,F9.1,F7.1,F7.1,F9.1,F7.1,F7.1,F9.1,F7.1,F7.1,F9.1)' 

	ENDIF 


	
;Set next startparm to this model
;startparm = modelresult


ENDFOR

CLOSE, 1


fitresultname = objname + 'fitresult.txt'
vel_slicename1 = objname + 'vel_slice1.fits'
veldisp_slicename1 = objname + 'veldisp_slice1.fits'
flux_slicename1 = objname + 'flux_slice1.fits'
vel_slicename2 = objname + 'vel_slice2.fits'
veldisp_slicename2 = objname + 'veldisp_slice2.fits'
flux_slicename2 = objname + 'flux_slice2.fits'
vel_slicename3 = objname + 'vel_slice3.fits'
veldisp_slicename3 = objname + 'veldisp_slice3.fits'
flux_slicename3 = objname + 'flux_slice3.fits'
vel_slicename4 = objname + 'vel_slice4.fits'
veldisp_slicename4 = objname + 'veldisp_slice4.fits'
flux_slicename4 = objname + 'flux_slice4.fits'





;WRITECOL, fitresultname, xnode, ynode, fit_result[0,*], fit_result[1,*], fit_result[2,*], fit_result[3,*]

WRITEFITS, vel_slicename1, vel_slice1
WRITEFITS, veldisp_slicename1, veldisp_slice1
WRITEFITS, flux_slicename1, flux_slice1

IF ((compnum EQ 2) OR (compnum EQ 3) OR (compnum EQ 4)) THEN BEGIN

WRITEFITS, vel_slicename2, vel_slice2
WRITEFITS, veldisp_slicename2, veldisp_slice2
WRITEFITS, flux_slicename2, flux_slice2



;Make flux weighted velocity map
;sumflux = flux_slice1 + flux_slice2
;ratio1 = flux_slice1 / sumflux
;ratio2 = flux_slice2 / sumflux
;nonzeroidx = WHERE(sumflux GT 0)
;avgvel = FLTARR(cubex,cubey)
;avgvel[nonzeroidx] = (vel_slice1[nonzeroidx] * ratio1[nonzeroidx]) + (vel_slice2[nonzeroidx] * ratio2[nonzeroidx])
;avgveldisp = FLTARR(cubex,cubey)
;avgveldisp[nonzeroidx] = (veldisp_slice1[nonzeroidx] * ratio1[nonzeroidx]) + (veldisp_slice2[nonzeroidx] * ratio2[nonzeroidx])

;avgvelname = objname + 'avgvelslice.fits'
;avgdispname = objname + 'avgdsisplice.fits'

;WRITEFITS, avgdispname, avgveldisp
;WRITEFITS, avgvelname, avgvel

ENDIF

IF ((compnum EQ 3) OR (compnum EQ 4)) THEN BEGIN

WRITEFITS, vel_slicename3, vel_slice3
WRITEFITS, veldisp_slicename3, veldisp_slice3
WRITEFITS, flux_slicename3, flux_slice3

ENDIF

IF (compnum EQ 4) THEN BEGIN

WRITEFITS, vel_slicename4, vel_slice4
WRITEFITS, veldisp_slicename4, veldisp_slice4
WRITEFITS, flux_slicename4, flux_slice4

ENDIF

END





PRO gmoscube_process, objnum

IF (objnum EQ 1) THEN BEGIN

; 3C303.1
; Had to turn off median clipping for this object since some lines are too narrow (low velocity dispersion)

objname = '3C303_1'
trimcubename = objname + 'trimcube.fits'
errcubename = objname + 'trimcube_err.fits'

b_begin = 4530
b_end = 4630
r_begin = 5250
r_end = 5350
;filename = '3C303_1nocal.fits'
filename = '3C303_1cal.fits'
o3_begin = 5070
o3_end = 5213
redshift = 0.2704

targetSN = 100.0
minsn = 30.0

parinfo_struct = {parinfostruct}
parinfo = REPLICATE(parinfo_struct, 7)

parinfo[*].value = [0.0, 70.06, 15.0, 2.E-16, 78.06, 5.0, 8.E-16]
parinfo[1].limited=[1,1]
parinfo[1].limits=[55.0,85.0]
parinfo[2].limited=[1,1]
parinfo[2].limits=[5.0,25.0]
parinfo[3].limited=[1,0]
parinfo[3].limits=[0.0,100000000.0]
parinfo[4].limited=[1,1]
parinfo[4].limits=[30.0,110.0]
parinfo[5].limited=[1,1]
parinfo[5].limits=[0.1,10.0]
parinfo[6].limited=[1,0]
parinfo[6].limits=[0.0,100000000.0]


;MKO3SUB, filename, objname, redshift, b_begin, b_end, r_begin, r_end, o3_begin, o3_end
;MAKECUBESLICE, trimcubename, objname, minsn
;DOVORONOI, objname, targetSN
CUBE_FITTING, trimcubename, errcubename, objname, parinfo, redshift, 2

ENDIF


IF (objnum EQ 2) THEN BEGIN

;PKS1306-09
;This object needs two components fit

objname = 'PKS1306_09'
trimcubename = objname + 'trimcube.fits'
errcubename = objname + 'trimcube_err.fits'

b_begin = 4800
b_end = 4900
r_begin = 5600
r_end = 5700
;filename = 'pks1306_09nocal.fits'
filename = 'pks1306_09cal.fits'
o3_begin = 5385
o3_end = 5540
redshift = 0.46685

targetSN = 100.0
minsn = 30.0

;startparm = [0.0, 58.25, 10.0, 6000.0, 93.67, 10.0, 1500.0]

parinfo_struct = {parinfostruct}
parinfo = REPLICATE(parinfo_struct, 7)

parinfo[*].value = [0.0, 58.25, 10.0, 5.0E-17, 93.67, 10.0, 4.2E-17]
parinfo[1].limited=[1,1]
parinfo[1].limits=[50.0,68.0]
parinfo[4].limited=[1,1]
parinfo[4].limits=[85.0,105.0]

;MKO3SUB, filename, objname, redshift, b_begin, b_end, r_begin, r_end, o3_begin, o3_end
;MAKECUBESLICE, trimcubename, objname, minsn
;DOVORONOI, objname, targetSN
CUBE_FITTING, trimcubename, errcubename, objname, parinfo, redshift, 2

ENDIF




IF (objnum EQ 3) THEN BEGIN

;3C 49
;four components....

objname = '3C49'
trimcubename = objname + 'trimcube.fits'
errcubename = objname + 'trimcube_err.fits'

b_begin = 3080
b_end = 3155
r_begin = 3500
r_end = 3580
;filename = '3C49nocalcomb.fits'
filename = '3C49calcomb_n.fits'
o3_begin = 3300
o3_end = 3450
redshift = 0.621

targetSN = 100.0
minsn = 50.0

parinfo_struct = {parinfostruct}
parinfo = REPLICATE(parinfo_struct, 13)
;parinfo = REPLICATE(parinfo_struct, 10)
;parinfo = REPLICATE(parinfo_struct, 7)
;parinfo = REPLICATE(parinfo_struct, 4)

;parinfo[*].value = [0.0, 65.2, 15.0, 1.06E-16, 90.2, 3.0, 1.06E-16]
;parinfo[*].value = [0.0, 70.2, 10.0, 1.06E-16, 70.2, 1.0, 5.06E-16, 88.2, 1.0, 5.06E-17]
parinfo[*].value = [0.0, 70.2, 10.0, 1.06E-16, 70.2, 1.0, 5.06E-16, 80.2, 1.0, 5.06E-17, 95.2, 1.0, 5.06E-17]

parinfo[0].limited=[1,0]
parinfo[0].limits=[0.0,100000.0]
parinfo[1].limited=[1,1]
parinfo[1].limits=[50.0,90.0]
parinfo[2].limited=[1,1]
parinfo[2].limits=[8.0,25.0]
parinfo[3].limited=[1,0]
parinfo[3].limits=[0.0,10000000.0]
parinfo[4].limited=[1,1]
parinfo[4].limits=[50.0,75.0]
parinfo[5].limited=[1,1]
parinfo[5].limits=[0.5,7.0]
parinfo[6].limited=[1,0]
parinfo[6].limits=[0.0,1.0]
parinfo[7].limited=[1,1]
parinfo[7].limits=[75.0,85.0]
parinfo[8].limited=[1,1]
parinfo[8].limits=[0.5,7.0]
parinfo[9].limited=[1,0]
parinfo[9].limits=[0.0,1.0]
parinfo[10].limited=[1,1]
parinfo[10].limits=[90.0,130.0]
parinfo[11].limited=[1,1]
parinfo[11].limits=[0.3,7.0]
parinfo[12].limited=[1,0]
parinfo[12].limits=[0.0,1.0]



;MKO3SUB, filename, objname, redshift, b_begin, b_end, r_begin, r_end, o3_begin, o3_end
;MAKECUBESLICE, trimcubename, objname, minsn
;DOVORONOI, objname, targetSN
CUBE_FITTING, trimcubename, errcubename, objname, parinfo, redshift, 4, fix_ctr_coord=[9,9]

ENDIF



IF (objnum EQ 4) THEN BEGIN

;PKS0023-26

objname = 'PKS0023'
trimcubename = objname + 'trimcube.fits'
errcubename = objname + 'trimcube_err.fits'

b_begin = 4900
b_end = 5030
r_begin = 5350
r_end = 5500
;filename = 'PKS0023nocalcomb.fits'
filename = 'PKS0023calcomb_n.fits'
o3_begin = 5190
o3_end = 5320
redshift = 0.32162

targetSN = 100.0
minsn = 65.0

parinfo_struct = {parinfostruct}
parinfo = REPLICATE(parinfo_struct, 7)

parinfo[*].value = [0.0, 62.7, 5.0, 5.7E-16, 69.0, 15, 2.E-15]
parinfo[0].limited=[1,0]
parinfo[0].limits=[0.0,100000.0]
parinfo[1].limited=[1,1]
parinfo[1].limits=[45.0,75.0]
parinfo[2].limited=[1,1]
parinfo[2].limits=[0.1,15.0]
parinfo[3].limited=[1,0]
parinfo[3].limits=[0.0,100000000.0]
parinfo[4].limited=[1,1]
parinfo[4].limits=[10.0,100.0]
parinfo[5].limited=[1,1]
parinfo[5].limits=[10.0,25.0]
parinfo[6].limited=[1,0]
parinfo[6].limits=[0.0,100000000.0]


;MKO3SUB, filename, objname, redshift, b_begin, b_end, r_begin, r_end, o3_begin, o3_end
;MAKECUBESLICE, trimcubename, objname, minsn
;DOVORONOI, objname, targetSN
CUBE_FITTING, trimcubename, errcubename, objname, parinfo, redshift, 2

ENDIF



IF (objnum EQ 5) THEN BEGIN

;PKS2135

objname = 'PKS2135'
trimcubename = objname + 'trimcube.fits'
errcubename = objname + 'trimcube_err.fits'

b_begin = 3230
b_end = 3300
r_begin = 3570
r_end = 3685
;filename = 'PKS2135nocalcomb.fits'
filename = 'PKS2135calcomb_n.fits'
o3_begin = 3435
o3_end = 3570
redshift = 0.63634	;Check redshift

targetSN = 100.0
minsn = 30.0

parinfo_struct = {parinfostruct}
parinfo = REPLICATE(parinfo_struct, 4)

parinfo[*].value = [0.0, 59.1, 17.0, 3.45E-15]
parinfo[3].limited=[1,0]
parinfo[3].limits=[0.0,100000000.0]

;MKO3SUB, filename, objname, redshift, b_begin, b_end, r_begin, r_end, o3_begin, o3_end
;MAKECUBESLICE, trimcubename, objname, minsn
;DOVORONOI, objname, targetSN
CUBE_FITTING, trimcubename, errcubename, objname, parinfo, redshift, 1

ENDIF




IF (objnum EQ 6) THEN BEGIN

;4C 39.56

objname = '4C3956'
trimcubename = objname + 'trimcube.fits'
errcubename = objname + 'trimcube_err.fits'

b_begin = 3890 
b_end = 3990  
r_begin = 4795
r_end = 4890 
;filename = '4C3956nocalcomb.fits'
filename = '4C3956calcomb_n.fits'
o3_begin = 4590
o3_end = 4720
redshift = 0.798	;Check redshift

targetSN = 100.0
minsn = 50.0


parinfo_struct = {parinfostruct}
parinfo = REPLICATE(parinfo_struct, 7)

parinfo[*].value = [0.0, 40.0, 3.0, 5E-16, 60.0, 20.0, 1E-15]
parinfo[0].limited=[1,0]
parinfo[0].limits=[0.0,100000.0]
parinfo[1].limited=[1,1]
parinfo[1].limits=[25.0,55.0]
parinfo[2].limited=[0,1]
parinfo[2].limits=[1.0,7.0]
parinfo[3].limited=[1,0]
parinfo[3].limits=[0.0,100000000.0]
parinfo[4].limited=[1,1]
parinfo[4].limits=[45.0,80.0]
parinfo[6].limited=[1,0]
parinfo[6].limits=[0.0,100000000.0]


;MKO3SUB, filename, objname, redshift, b_begin, b_end, r_begin, r_end, o3_begin, o3_end
;MAKECUBESLICE, trimcubename, objname, minsn
;DOVORONOI, objname, targetSN
CUBE_FITTING, trimcubename, errcubename, objname, parinfo, redshift, 2

ENDIF




IF (objnum EQ 7) THEN BEGIN

;3C 268.3

objname = '3C268_3'
;objname = 'cal3C268_3'
trimcubename = objname + 'trimcube.fits'
errcubename = objname + 'trimcube_err.fits'

b_begin = 4700
b_end = 4800
r_begin = 5400
r_end = 5500
;filename = '3C268_3nocal.fits'
filename = '3C268_3cal.fits'
o3_begin = 5260
o3_end = 5420
redshift = 0.37171

targetSN = 100.0
minsn = 30.0


parinfo_struct = {parinfostruct}
parinfo = REPLICATE(parinfo_struct, 7)

parinfo[*].value = [0.0, 72.0, 20.0, 4.15e-16, 73.02, 5.0, 8.8e-17]
parinfo[1].limited=[1,1]
parinfo[1].limits=[55.0,85.0]
parinfo[2].limited=[1,1]
parinfo[2].limits=[10.0,25.0]
parinfo[3].limited=[1,0]
parinfo[3].limits=[0.0,100000000.0]
parinfo[4].limited=[1,1]
parinfo[4].limits=[30.0,130.0]
parinfo[5].limited=[1,1]
parinfo[5].limits=[1.0,10.0]



;MKO3SUB, filename, objname, redshift, b_begin, b_end, r_begin, r_end, o3_begin, o3_end
;MAKECUBESLICE, trimcubename, objname, minsn
;DOVORONOI, objname, targetSN
CUBE_FITTING, trimcubename, errcubename, objname, parinfo, redshift, 2, fix_ctr_coord=[12,11]

ENDIF




IF (objnum EQ 8) THEN BEGIN

;3C 93.1

objname = '3C93_1'
;objname = 'cal3C93_1'
trimcubename = objname + 'trimcube.fits'
errcubename = objname + 'trimcube_err.fits'

b_begin = 4860
b_end = 4920
r_begin = 5225
r_end = 5285
;filename = '3C93_1nocalcomb.fits'
filename = '3C93_1calcomb.fits'
o3_begin = 5080
o3_end = 5210
redshift = 0.243

targetSN = 100.0
minsn = 50.0


parinfo_struct = {parinfostruct}
parinfo = REPLICATE(parinfo_struct, 7)

parinfo[*].value = [0.0, 60.0, 15.0, 1.e-17, 60.0, 5.0, 1.e-17]
parinfo[1].limited=[1,1]
parinfo[1].limits=[40.0,80.0]
parinfo[2].limited=[1,1]
parinfo[2].limits=[10.0,20.0]
parinfo[3].limited=[1,0]
parinfo[3].limits=[0.0,100000000.0]
parinfo[4].limited=[1,1]
parinfo[4].limits=[40.0,80.0]
parinfo[5].limited=[1,1]
parinfo[5].limits=[1.0,10.0]



;MKO3SUB, filename, objname, redshift, b_begin, b_end, r_begin, r_end, o3_begin, o3_end
;MAKECUBESLICE, trimcubename, objname, minsn
;DOVORONOI, objname, targetSN
CUBE_FITTING, trimcubename, errcubename, objname, parinfo, redshift, 2, fix_ctr_coord=[12,11]

ENDIF





IF (objnum EQ 9) THEN BEGIN

;3C 124
;z > 1, fitting O II 3727 line instead of O III 5007

objname = '3C124'
;objname = 'cal3C124'
trimcubename = objname + 'trimcube.fits'
errcubename = objname + 'trimcube_err.fits'

b_begin = 2600
b_end = 2770
r_begin = 2950
r_end = 3120
;filename = '3C124nocalcomb.fits'
filename = '3C124calcomb.fits'
o3_begin = 2770
o3_end = 2900
redshift = 1.083


targetSN = 100.0
minsn = 50.0


parinfo_struct = {parinfostruct}
parinfo = REPLICATE(parinfo_struct, 7)

parinfo[*].value = [0.0, 60.0, 10.0, 1.e-17, 70.0, 10.0, 1.e-17]
parinfo[1].limited=[1,1]
parinfo[1].limits=[40.0,65.0]
parinfo[2].limited=[1,1]
parinfo[2].limits=[5.0,22.0]
parinfo[3].limited=[1,0]
parinfo[3].limits=[0.0,100000000.0]
parinfo[4].limited=[1,1]
parinfo[4].limits=[65.0,90.0]
parinfo[5].limited=[1,1]
parinfo[5].limits=[5.0,22.0]
parinfo[6].limited=[1,0]
parinfo[6].limits=[0.0,100000000.0]



;MKO3SUB, filename, objname, redshift, b_begin, b_end, r_begin, r_end, o3_begin, o3_end
;MAKECUBESLICE, trimcubename, objname, minsn
;DOVORONOI, objname, targetSN
CUBE_FITTING, trimcubename, errcubename, objname, parinfo, redshift, 2, fix_ctr_coord=[12,11]

ENDIF




END
