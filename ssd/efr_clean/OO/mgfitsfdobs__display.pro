;----------------------------------------------------------------------------
; MgFitsFDObs::Display
;
; Purpose:
;  Retrieves the value of properties associated with the orb object.
;

		PRO MgFitsFDObs::Display, window=window
;
;		;	self->GetProperty, image=image
;
               XSIZE=800
               YSIZE=600
	if keyword_set(window) then window, 0, xs=XSIZE, ys=YSIZE ;xs=self.nx, ys=self.ny
			tvframe, *self.display, /asp
;		   tvscl, *self.display
;

;
		END
