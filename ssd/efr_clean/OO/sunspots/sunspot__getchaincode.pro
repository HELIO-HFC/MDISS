function sunspot:: getchaincode, code=code

; FUNCTION
;    writing penumbral chain code in arcsecs
;    for LMSAL catalogue


self->getbndrct, x0, y0, x1, y1

;print, x0, x1, y0, y1

dlt=5
nx0=x1-x0+1
ny0=y1-y0+1

nx1=nx0+2*dlt
ny1=ny0+2*dlt

imb=bytarr(self.nx, self.ny)
imb[*self.locs]=1
;if (*self.umbra)[0] ne -1 then imb[(*self.locs)[*self.umbra]]=2



;   xi=(findgen(2*self.nx)-2*self.xc)*self.cdx/2   
;   yi=(findgen(2*self.ny)-2*self.yc)*self.cdy/2
;    contour, congrid(imb, 2*self.nx, 2*self.ny) gt 0, xi, yi, path_xy=xy, /path_data_coord, level=[.99], /path_double
;stop
imc=imb[x0:x1, y0:y1]
imc1=imb[x0-dlt : x1+dlt, y0-dlt : y1+dlt]

in=size(imc1)
nx0=in(1) & ny0=in(2)

imc=imc1 & imc(*)=0
for i=1, nx0-2 do for j=1, ny0-2 do begin
   da=imc1(i-1:i+1, j) 
   da=da(*) 
   da=[da, imc1(i, j+1), imc1(i, j-1)]
   if (min(da) ne max(da)) and (imc1(i, j) eq 1) $
     then imc(i, j)=1
end 

locs=where(imc ne 0)
xp=x0-dlt+locs mod nx0
yp=y0-dlt+locs / nx0

nl=n_elements(locs)
xy1=dblarr(2, nl)

xy1(0, *)=(xp-self.xc)*self.cdx
xy1(1, *)=(yp-self.yc)*self.cdy
;stop
;		sc=5
;		window, 1, xs=sc*nx1, ys=sc*ny1, title='display'
;		tvscl, congrid(imc1, sc*nx1,sc*ny1)
;
;
;		window, 2, xs=sc*nx0, ys=sc*ny0, title='bounding'
;		tvscl, congrid(imc, sc*nx0,sc*ny0)
;


;stop
return, xy1
end


