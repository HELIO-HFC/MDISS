;----------------------------------------------------------------------------
; FILE::GETSTRUCTURE
;
; Purpose:
;  Retrieves the value of properties associated with the orb object.
;

		FUNCTION FullDiskObs::GetStructure
;
;		;	self->GetProperty, image=image
;
			struct={FDObs, R:self.R, $
							xc:self.xc,$
							yc:self.yc, $
							cdx:self.cdx, $
							cdy:self.cdy, $
							nx:self.nx, $
							ny:self.ny,  $
							date:self.date, $



							origin:self.origin}
;
			return, struct
;
		END
