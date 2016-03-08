function extractMDI_heliographic, map, im, lon, lat, helio_size,$
	gauss=gauss, outmap=outmap, res=res


;	FUNCTION
;		to remap FD MDI observation into HelioGraphic plane

; *** size of the heilographic
;	  image n x n degrees
;		extra two degrees added in each direction
;		on each side

;helio_size=20.88
;helio_size=20.88

if not keyword_set(res) then res=.125

nx=fix(helio_size/res)
ny=fix(helio_size/res)

image=intarr(nx, ny)

;mreadfits, fn, index, im, header=hd

;im=map.data
time=map.date
xc=map.xc
yc=map.yc
cdx=map.cdx
cdy=map.cdy


lon0=tim2carr(time)
flare_arc=hel2arcmin(lat, lon-lon0,  $
			date=time, /soho)

;flare_x=xc+flare_arc[0]*60/cdx
;flare_y=yc+flare_arc[1]*60/cdy


; *** lower left corner coordinates
ll_lon=lon-lon0-helio_size/2.
ll_lat=lat-helio_size/2.

; *** upper right corner coordinates
ur_lon=lon-lon0+helio_size/2.
ur_lat=lat+helio_size/2.

;stop
lon_array=dblarr(nx, ny)
lat_array=dblarr(nx, ny)
for i=0, nx-1 do lon_array(i, *)=ll_lon+i*res
for i=0, ny-1 do lat_array(*, i)=ll_lat+i*res


arcarr=hel2arcmin(lat_array, lon_array, date=time, /soho)
xarr=xc+arcarr[0, *]*60/cdx
yarr=yc+arcarr[1, *]*60/cdy


in=size(im)
nx1=in(1) & ny1=in(2)
xarr=reform(fix(xarr+.5))
yarr=reform(fix(yarr+.5))


for i=0l,long(nx)*ny-1 do begin
	if xarr(i) ge 0 and xarr(i) lt nx1 and $
		yarr(i) ge 0 and yarr(i) lt ny1 then $
	image[i]=im(xarr(i), yarr(i)) else image[i]=0
endfor
;	if keyword_set(gauss) then image=gauss_smoothing(image, 2)
;	image=image(5:nx-6, 5:nx-6)
;	;window, 0, xs=160, ys=160
;	;tvscl, image
;	;print, flare_x, flare_y
;
;	dx=(nx-128)/2
;	dy=(ny-128)/2
;
;	if keyword_set(cut128) then image=image[dx:nx-dx-1, dy:ny-dy-1]

;outmap=map
;outmap.data=image
;stop
return, image

end

