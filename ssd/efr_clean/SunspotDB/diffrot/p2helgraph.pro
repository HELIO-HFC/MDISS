function p2helgraph, x, y, xc, yc, cdelt1, cdelt2, date


;
;		INPUTS
;			x, y 	pixel coordinates
;			xc, yc	disk center coordinates (pixel)
;			cdelti	scaling (arcsec per pixel)
;			date	date of observation

hg=fltarr(2)

clon=tim2carr(date)

pbr=pb0r(date, /earth)
print, pbr
hg=arcmin2hel((x-xc)*cdelt1/60., (y-yc)*cdelt2/60., date=date, $
				l0=clon, p=pbr[0], b0=pbr[1], r0=1.0/tan(pbr(2)/60/!radeg))

if hg[1] lt 0 then hg[1]=hg[1]+360.



return, hg
end