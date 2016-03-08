pro map2hel, image, header, helarr, head2, delta, wdtx, wdty


;	INPUT
;		image		original image
;		header	    the corresponding header
;		ns, ew		heliographic coordinates of the centre of the area
;							to be extracted/mapped
;		wdt			width of the area (in degrees)


date=fxpar(header, 'DATE_OBS')
cd1=fxpar(header, 'CDELT1')
cd2=fxpar(header, 'CDELT2')
R=fxpar(header, 'Solar_R')
xc=fxpar(header, 'CENTER_X')
yc=fxpar(header, 'CENTER_Y')

nx=fxpar(header, 'NAXIS1')
ny=fxpar(header, 'NAXIS2')


;	Compute p, b0 and r to save time

;delta=.1

N1d=2d00*double(wdtx)/delta
N2d=2d00*double(wdty)/delta

N1=fix(N1d+.5)
N2=fix(N2d+.5)

helarr=intarr(N1, N2)
angle=pb0r(date, /earth)
rx=1.0/tan(angle(2)/60/!radeg)
clon=tim2carr(date)

 b0 = angle[1]/!radeg

ns=0.d00
ew=clon[0]


;	Take width equal to 30 degrees for time being

wdtx=double(wdtx)
wdty=double(wdty)

;	Image Dimensions
Sf=4.d00


delta1=2d00*double(wdtx)/N1 	; degrees per pixel
delta2=2d00*double(wdty)/N2 	; degrees per pixel

print, delta1, delta2

;	Prepare the header for writing the mapped image

head2=header

sxdelpar, head2, 'OBJECT_S'
sxdelpar, head2, 'OBS_MODE'
sxdelpar, head2, 'SOLAR_R'
sxdelpar, head2, 'CENTER_X'
sxdelpar, head2, 'CENTER_Y'
sxdelpar, head2, 'CRPIX1'
sxdelpar, head2, 'CRPIX2'
sxdelpar, head2, 'XSCALE'
sxdelpar, head2, 'YSCALE'

fxaddpar, head2, 'LATITUDE', double(ns), 'image cenre heliograph coordinates'
fxaddpar, head2, 'LONGITUD', double(ew), 'image cenre heliograph coordinates'
fxaddpar, head2, 'CTYPE1','DEGREES'
fxaddpar, head2, 'CDELT1', delta1, 'Scale along X_axis in degrees/pixel'
fxaddpar, head2, 'CDELT2', delta2, 'Scale along Y_axis in degrees/pixel'

fxaddpar, head2, 'CRVAL1', double(ew)-wdtx
fxaddpar, head2, 'CRVAL2', double(ns)-wdty

print, date, cd1, cd2, R, xc, yc

print, R*cd1, R*cd2

rad1=R*cd1
rad2=R*cd2



dx=0

; trying to compensate for b0 angle

b0offset=fix(angle[1]/delta+.5)

for i=0, N1-1 do $
	for j=0, N2-1 do begin

		lon=ew-wdtx-clon+i*delta1
		lat=ns-wdty+delta2*j +angle[1]

		arcsec=60*hel2arcmin_modified(lat, lon,  b0, angle[0], rx, angle[2])
	;	arcsec=60*hel2arcmin(lat, lon,  b0=angle[1], p=angle[0], r0= rx)


		if abs(arcsec[0]) gt rad1 or abs(arcsec[1]) ge rad2 then helarr[i, j]=0 $
			else begin
					xn=fix(arcsec[0]/cd1+xc+.5)
					yn=fix(arcsec[1]/cd2+yc+.5)
					helarr[i, j]=image[xn, yn]
				endelse
endfor


window, 21, xs=N1+2*dx, ys=N2+2*dx
tvscl, helarr

;writefits, 'Y:/meudon/processed/heliotest-cak1.fits', helarr, head2


end