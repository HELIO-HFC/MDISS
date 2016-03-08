;----------------------------------------------------------------------------
; FullDiskObs::Display
;
; Purpose:
;  Retrieves the value of properties associated with the orb object.
;

		PRO FullDiskObs::Display, window=window
;
;		;	self->GetProperty, image=image
                XSIZE=800
                YSIZE=600
		if keyword_set(window) then window, 0, xs=XSIZE, ys=YSIZE ;xs=self.nx, ys=self.ny
		;	tvscl,  *self.image
			tvframe, *self.image, /asp
;

;
		END
