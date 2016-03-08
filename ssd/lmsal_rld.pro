pro lmsal_rld, im, hd, output=output, qs_val=qs_val
;   PROCEDURE
;      to remove limb darkening from the image
;
;   INPUTS
;      im       image (white-light)
;      hd       header (TRACE)


cdx=fxpar(hd, 'CDELT1')      
cdy=fxpar(hd, 'CDELT2')   


xc=fxpar(hd, 'CRPIX1')
yc=fxpar(hd, 'CRPIX2')

in=size(im)
nx=in(1) & ny=in(2)

xi=findgen(nx)#(fltarr(ny)+1.)
yi=(fltarr(nx)+1.)#findgen(ny)

ra=sqrt((xi-xc)*(xi-xc)+(yi-yc)*(yi-yc))

phi=atan( yi-yc, xi-xc)

;phi=atan( xi-xc, yi-yc)
rs=dblarr(nx, ny)

ri=(min(ra)+findgen(nx)/(nx-1)*(max(ra)-min(ra)))#(fltarr(ny)+1)
pi=(fltarr(nx)+1.)#(min(phi)+findgen(ny)/(ny-1)*(max(phi)-min(phi)))

xx=ri*cos(pi)
yy=ri*sin(pi)

rs=interpolate(im, xx+xc, yy+yc, miss=-1.)

qsr=dblarr(nx, ny)
for i=0, nx-1 do begin
    line=reform(rs(i, *))
    locs=where(line ne -1)
    if locs(0) ne -1 then $
      qsr(i, *)=median(line(locs))
end

rr=(max(ra)-min(ra))/(nx-1)
pr=(max(phi)-min(phi))/(ny-1)

qs=interpolate(qsr, (ra-min(ra))/rr, (phi-min(phi))/pr, miss=0.)
;tvframe, qs, /bar, /asp


output=im-qs
output=output-min(output)

;tvframe, output, /bar, /asp

h=histogram(output, loc=hi, nbins=500)

locs=where(h eq max(h))
qs_val=hi(locs(0))
;stop
end


