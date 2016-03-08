FUNCTION MgFitsFDObs::Init, location=location;, map=map


;
;
;	27-07-2004	added Solar_P0 angle compensation
;				removed map object references
;
if keyword_set(location) then $
	x=self->FitsFDObs::Init(location=location) $
	else x=self->FitsFDObs::Init(map=map)

;	correct SOLAR_P angle if necessary

	self->getheader, hd
	solar_p=fxpar(hd, 'SOLAR_P0')

	if (solar_p ne 0) then begin
		xc=fxpar(hd, 'X0')
		yc=fxpar(hd, 'Y0')
		im=*self.image
		im=rot(temporary(im), solar_p, 1, xc, yc)
		*self.image=im
	 endif
;stop
		disp=bytscl(*self.image, max=1000, min=-1000)
		self.display=ptr_new(disp)

	return, 1
END

;========================================================

PRO MgFitsFDObs::GetAR, Ararr, display=display, df=df, rect=rect

if not keyword_set(df) then df=4
if not keyword_set(rect) then $
	if not keyword_set(display) $
		then ARbin=MDI_AR(*self.image, *self.header, df=df) $
		else ARbin=MDI_AR(*self.image, *self.header, /display, df=df) $
	else ARbin=MDI_AR(*self.image, *self.header, df=df, /display, /rect)

LabelCountRegion, ARbin, n, ploc
print, 'Number of Active Regions:', n

Ararr=objarr(n)

for i=0, n-1 do begin
	Ararr[i]=obj_new('ActiveRegion', self, *ploc[i])
endfor
END

;=============================================================

PRO MgFitsFDObs::Cleanup

ptr_free, self.header
ptr_free, self.display
ptr_free, self.image

END

;----------------------------------------------------------------------------
; MgFitsFDObs__DEFINE
;
; Purpose:
;  Defines the object structure for an orb object.
;
PRO MgFitsFDObs__define

    struct = { MgFitsFDObs, $
    			display:ptr_new(/allocate_heap), $
    			inherits FitsFDObs $
             }
END
