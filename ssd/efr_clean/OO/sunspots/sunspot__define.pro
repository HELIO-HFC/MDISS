FUNCTION Sunspot::Init, Fd, locs, umbra, n, mean0, res=res


struct_assign, Fd, self, /nozero


self.locs=ptr_new(locs)
self.n_umbras=n

;			self.nx=fd.nx
;			self.ny=fd.ny
;			self.cdx=fd.cdx
;			self.cdy=fd.cdy
;			self.xc=fd.xc
;			self.yc=fd.yc
;			self.R=fd.R
;			self.date=fd.date
;			self.origin=fd.origin

self.umbra=ptr_new(umbra)
self.meanInt=mean0
im=self->pix2hel(area=area, diameter=diam, res=res)

self->getbarycenter, fd, xm, ym

self.GCarcx=(xm-self.xc)*self.cdx
self.GCarcy=(ym-self.yc)*self.cdy

;print, self.GCarcx, self.GCarcy
self.helArea=area
self.helDiam=diam
	return, 1
END

;=======================================

PRO Sunspot::SetMagneticFlux, mgim

if not keyword_set(mgim) then mgim=mg->getimage()

locs=*self.locs
self.maxflux=max(mgim(locs))
self.minflux=min(mgim(locs))
self.totFlux=total(mgim(locs))
self.absTotFlux=total(abs(mgim(locs)))

umbra=*self.umbra
if umbra[0] ne -1 then begin
  ulocs=locs(umbra)
  self.umbflux=total(mgim(ulocs))
  self.absUmbFlux=total(abs(mgim(ulocs)))
  self.maxuflux=max(mgim(ulocs))
  self.minuflux=min(mgim(ulocs))
endif

END



;=======================================

PRO Sunspot::SetMgFlux, mg, absFlux=absFlux, mgim=mgim

if not keyword_set(mgim) then mgim=mg->getimage()


locs=*self.locs

plus=where(mgim[locs] gt 0)
minus=where(mgim[locs] lt 0)

if n_elements(plus) gt n_elements(minus) then $
		polarity=1 else polarity=-1
;help, plus
;help, minus

if not keyword_set(absFlux) then $
	self.totFlux= total(mgim[locs]) else $
	self.totFlux=total(abs(mgim[locs]))*total(mgim[locs])/abs(total(mgim[locs]))
self.maxFlux= max(abs(mgim[locs]))*self.totFlux/abs(self.totFlux)
;print, 'Total Flux:', total(mgim[locs])

self.umbFlux=0.d00
umbra=*self.umbra
if umbra[0] ne -1 then begin
	if not keyword_set(absFlux) then $
		self.umbFlux=total(mgim[locs[umbra]]) else $
		self.UmbFlux=total(abs(mgim[locs[umbra]]))*total(mgim[locs[umbra]])/abs(total(mgim[locs[umbra]]))
	self.maxuFlux= max(abs(mgim[locs[umbra]]))*self.umbFlux/abs(self.umbFlux)

endif
;print, self.umbFlux

END


;====================================================

PRO Sunspot::GetHelio

im=self->pix2hel(area=area, /display)
self.helArea=area

END



;-----------------------------------------------------------------------------

FUNCTION Sunspot:: GetUmbra


return, *self.umbra

end


;==============================================================================

PRO Sunspot:: GetData, helArea=helArea, pixSize=pixSize, totFlux=totFlux, helDiam=helDiam

helArea= self.helArea
totFlux= self.totFlux
helDiam=self.helDiam
pixSize=n_elements(*self.locs)

END
;==============================================================================
PRO Sunspot:: DisplayGC, fdo, scale=scale

if not keyword_set(scale) then scale=10
self->cropdisplay, fdo, scale=scale, x0=x0, y0=y0
	gx=fix(self.GCarcx/self.cdx+self.xc+.5)
	gy=fix(self.GCarcy/self.cdy+self.yc+.5)

print, gx, x0
print, gy, y0
	xyouts, (gx-x0)*scale, (gy-y0)*scale, '+', color=200, /device

END



;==============================================================================
PRO Sunspot::Cleanup

ptr_free, self.umbra
ptr_free, self.locs

END
;----------------------------------------------------------------------------
; Sunspot__DEFINE
;
; Purpose:
;  Defines the object structure for an orb object.
;
PRO Sunspot__define

    struct = { Sunspot, $
    			umbra:ptr_new(/allocate_heap), $
    			QuietSunInt:0.d00, $
    			meanInt:0.d00, $
    			n_umbras:0, $
    			totFlux:0.d00, $
    			absTotFlux:0.d00, $
    			umbFlux:0.d00, $
    			absUmbFlux:0.d00, $
    			maxFlux:0.d00,	$
    			minFlux:0.d00, $
    			maxuFlux:0.d00, $
    			minuFlux:0.d00, $
    			helDiam:0.d00, $
    			helArea:0.d00, $
    			GCarcx:0.d00, $
    			GCarcy:0.d00, $
    			inherits Feature $
             }
END