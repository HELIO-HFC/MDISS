pro fr_sync_soho_obs, mg1, mg2, drot=drot



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

; **** calculate arc & carrington coord
;		for all pixels
arcx = ((locs mod m1.nx) - m1.xc)*m1.cdx/60.d00
arcx = ((locs / m1.nx) - m1.yc)*m1.cdy/60.d00
hel1=arcmin2hel(arcx, arcy, date=m1.date, /soho)

stop
locs=where(m2.data gt -32000)


xp= locs mod nx
yp= locs / nx

arcx= (xp-xc)*dx/60.d00
arcy= (yp-yc)*dy/60.d00

hl=arcmin2hel(arcx, arcy, date=map2.time)
window, 0, xs=1024, ys=1024
tvscl, bytscl(mgim, max=1000, min=-1000)

window, 1, xs=1024, ys=1024
tvscl, bytscl(mgim, max=1000, min=-1000)

end