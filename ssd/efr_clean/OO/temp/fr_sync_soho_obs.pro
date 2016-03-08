function fr_sync_soho_obs, mg1, mg2

;	FUNCTION
;		synchronising object mg2
;			to the point of view of object mg1

;if not keyword_set(drot) then $
;	map2_1=drot_map(map2, ref_map=map1) $
;	  else map2_1=map2


m1=mg1->getstructure()
m2=mg2->getstructure()

clon1=tim2carr(m1.date)
clon2=tim2carr(m2.date)

nx=m1.nx
ny=m1.ny

;mgim=m1.data
;mgim0=map2_1.data


pix_x=indgen(1024)
pix_y=pix_x

; **** first define the pixels of interest, i.e.
;		pixel inside the solar disk
;		so we don't have to do calculations for
;		any of the off limb data

pCa=dblarr(2, m1.nx, m1.ny)
for i=0, ny-1 do pca(0, *, i)=pix_x
for i=0, nx-1 do pca(1, i, *)=pix_y

dist2=(pca(0, *, *)-m1.xc)*(pca(0, *, *)-m1.xc) + $
	 	(pca(1, *, *)-m1.yc)*(pca(1, *, *)-m1.yc)

dist2=reform(temporary(dist2), 1024, 1024)
locs=where(dist2 lt m1.R*m1.R)
npix=n_elements(locs)
; **** calculate arc & carrington coord
;		for all pixels
arcx = ((locs mod m1.nx) - m1.xc)*m1.cdx/60.d00
arcy = ((locs / m1.nx) - m1.yc)*m1.cdy/60.d00
hel1=arcmin2hel(arcx, arcy, date=m1.date, /soho)

arc2=60*hel2arcmin(hel1[0, *], hel1[1, *] + clon1[0] - clon2[0], $
					 date=m2.date, /soho)

x2=reform(arc2[0, *]/m2.cdx + m2.xc, npix)
y2=reform(arc2[1, *]/m2.cdy + m2.yc, npix)


image2=mg2->getimage()
image1=image2
image1(*)=-32768.0
image3=image1
locsx = long(y2+.5) * m2.nx + long(x2+.5)
image1(locs)=image2(locsx)

;for i=0l, npix-1 do begin
;	image3(locs(i)) = $
;		image2( fix(x2(i)+.5), fix(y2(i)+.5))
; endfor

;window, 1, xs=m2.nx, ys=m2.ny
;tvscl, bytscl(image1, max=2000, min=-2000)

return, image1

end