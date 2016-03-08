function FullDiskObs::rot2sync, refo

;	Function
;			rotate the whole observation
;			to the time & parameters of the referrence observation
;			(all observations assumed soho)



st=self->getstructure()
;help, st, /st

rfst=refo->getstructure()
image=*self.image
rfimage=refo->getimage()
rimage=intarr(rfst.nx, rfst.ny)

locs=where(rfimage ne 0)
rimage(LOCS)= max(rfimage)

window, 0, xs=rfst.nx, ys=rfst.ny
tvscl, *self.image

window, 1, xs=rfst.nx, ys=rfst.ny
tvscl, rfimage



;				n=n_elements(locs)
;				xind=fltarr(n)
;				yind=fltarr(n)
;
;				xind= ((locs mod rfst.nx)-rfst.xc)*rfst.cdx/(60)
;				yind= ((locs / rfst.nx)-rfst.yc)*rfst.cdy/(60)
;
;				rclon=tim2carr(rfst.date)
;				clon=tim2carr(self.date)
;
;				;print, rclon, clon
;				;print, pb0r(rfst.date, /soho)
;				Hlcoord=arcmin2hel(xind, yind, date=rfst.date, /soho)
;				print, min(hlcoord[1, *]), max(hlcoord[1, *])
;
;				hlcoord[1, *]=hlcoord[1, *] + rclon[0] - clon[0]
;;				print, min(hlcoord[1, *]), max(hlcoord[1, *])
;
;				arcoord=hel2arcmin(hlcoord[0, *], hlcoord[1, *], date=self.date, /soho)
;
;				;help, arcoord
;				sxind=(arcoord[0, *]*60)/self.cdx+self.xc
;				syind=(arcoord[1, *]*60)/self.cdy+self.yc
;kj
;				;stop
;				for i=0l, n-1 do begin
;;					 rimage[locs[i]]=image[fix(sxind[i]+.5), fix(syind[i]+.5)]
;				end
;

rimage=rotlocs(locs, *self.image, self.date, self.xc, self.yc, self.cdx,$
				self.cdy, rfst.date, rfst.nx, rfst.ny, rfst.xc, rfst.yc, $
					rfst.cdx, rfst.cdy )

window, 3, xs=st.nx, ys=st.ny
tvscl, rimage


return, rimage

end

;=======================================================


function rotlocs, locs, image, date, xc, yc, cdx, cdy, $
						date2, nx2, ny2, xc2, yc2, cdx2, cdy2

;	FUNCTION
;		2 rotate given pixel locations to
;		a set time and set of other parameters


clon=tim2carr(date)
clon2=tim2carr(date2)


n=n_elements(locs)
xind=fltarr(n)
yind=fltarr(n)

xind= ((locs mod nx2)-xc2)*cdx2/(60)
yind= ((locs / nx2)-yc2)*cdy2/(60)

Hlcoord=arcmin2hel(xind, yind, date=date2, /soho)
;print, min(hlcoord[1, *]), max(hlcoord[1, *])

hlcoord[1, *]=hlcoord[1, *] + clon2[0] - clon[0]
;print, min(hlcoord[1, *]), max(hlcoord[1, *])

arcoord=hel2arcmin(hlcoord[0, *], hlcoord[1, *], date=date, /soho)

;help, arcoord
sxind=(arcoord[0, *]*60)/cdx+xc
syind=(arcoord[1, *]*60)/cdy+yc

rimage=intarr(nx2, ny2)
rimage[locs]=image[fix(sxind+.5), fix(syind+.5)]


return, rimage
end