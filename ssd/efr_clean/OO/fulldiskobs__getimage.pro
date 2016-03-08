;----------------------------------------------------------------------------
; FullDiskObs::Display
;
; Purpose:
;  Retrieves the value of properties associated with the orb object.
;

		FUNCTION FullDiskObs::GetImage
;
;		;	self->GetProperty, image=image


;
		return, *self.image
;
		END

