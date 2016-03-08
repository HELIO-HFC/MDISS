FUNCTION Feature::Init, Fd, locs


;if not keyword_set(locs) then return, 0

self.locs=ptr_new(locs)
;fd=FitsFDObs->getstructure()
self.nx=fd.nx
self.ny=fd.ny
self.cdx=fd.cdx
self.cdy=fd.cdy
self.xc=fd.xc
self.yc=fd.yc
self.R=fd.R
self.date=fd.date
self.origin=fd.origin
	return, 1
END



;=======================================================

FUNCTION Feature::PSize

return, n_elements(*self.locs)

END







;====================================================

PRO Feature::GetBaryCenter, fd, xmom, ymom, rev=rev


im=fd->getimage()
n=n_elements(*self.locs)
locs=*self.locs

if keyword_set(rev) then mx=max(im[locs])
;print, n
sum=0.D0
xmom=0.D0
ymom=0.D0
if not keyword_set(rev) then $
	for i=0l,n-1,1 do begin
		sum=sum+im[locs[i]]
		xmom=xmom+double(im[locs[i]])*double(locs[i] mod self.nx)
		ymom=ymom+double(im[locs[i]])*double(locs[i] / self.nx)
	end else $

	for i=0l,n-1,1 do begin
		sum=sum+mx-im[locs[i]]+1
		xmom=xmom+double(mx-im[locs[i]]+1)*double(locs[i] mod self.nx)
		ymom=ymom+double(mx-im[locs[i]]+1)*double(locs[i] / self.nx)
	end
xmom=xmom/sum
ymom=ymom/sum


END

;---------------------------------------------------------------------


FUNCTION Feature::GetDate

return, self.date

END

;----------------------------------------------------------------------------
; Feature__DEFINE
;
; Purpose:
;  Defines the object structure for an orb object.
;
PRO Feature__define

    struct = { Feature, $
    			locs:ptr_new(/allocate_heap), $
    			inherits FDObs $
             }
END