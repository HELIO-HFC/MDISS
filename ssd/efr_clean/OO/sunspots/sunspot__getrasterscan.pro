function sunspot:: GetRasterScan, code=code

self->getbndrct, x0, y0, x1, y1

;print, x0, x1, y0, y1

dlt=5
nx0=x1-x0+1
ny0=y1-y0+1

nx1=nx0+2*dlt
ny1=ny0+2*dlt

imb=bytarr(self.nx, self.ny)
imb[*self.locs]=1
if (*self.umbra)[0] ne -1 then imb[(*self.locs)[*self.umbra]]=2

;stop
imc=imb[x0:x1, y0:y1]
imc1=imb[x0-dlt : x1+dlt, y0-dlt : y1+dlt]


;		sc=5
;		window, 1, xs=sc*nx1, ys=sc*ny1, title='display'
;		tvscl, congrid(imc1, sc*nx1,sc*ny1)
;
;
;		window, 2, xs=sc*nx0, ys=sc*ny0, title='bounding'
;		tvscl, congrid(imc, sc*nx0,sc*ny0)
;

;	end preliminaries
; *****************************

x=image2raster(imc, code=code)

;help, imc


;print, 'String Length', strlen(x)

;imx=raster2image(x, nx0)

;if max(imx-imc) ne 0 or  min(imx-imc) ne 0 then stop
;stop

return, x
end

