		FUNCTION FullDiskObs::Init, FitsFile

		if not arg_present(FitsFile) then FitsFile=obj_new('FtsFile')
		FitsFile->GetProperty, image=image, header=header
		st=FitsFile->GetObs()
		struct_assign, st, self, /nozero

		self.image=ptr_new(image)
		self.bitpix=fxpar(header, 'BITPIX')

		obj_destroy, fitsfile
		return, 1
;
		END



;---------------------------------------------------------------------

PRO FullDiskObs::CropDisplay, x0, y0,  x1, y1, scale=scale

im=*self.image

imC=im[x0:x1, y0:y1]
if not keyword_set(scale) then scale=1

window, 0, xs=(x1-x0)*scale, ys=(y1-y0)*scale
tvscl, congrid(imc, (x1-x0)*scale, (y1-y0)*scale)

END


;----------------------------------------------------------------------------
; FullDiskObs::Display
;
; Purpose:
;  Retrieves the value of properties associated with the orb object.
;

		FUNCTION FullDiskObs::GetDisplay
;
;		;	self->GetProperty, image=image

			 return, *self.image
;

;
		END




;----------------------------------------------------------------------------
; FullDiskObs::Display
;
; Purpose:
;  Retrieves the value of properties associated with the orb object.
;

		FUNCTION FullDiskObs::GetMgImage, max=max0, min=min0
;
;		;	self->GetProperty, image=image

		if not keyword_set(max0) then max0=1000
		if not keyword_set(min0) then min0=-1000
		 im=bytscl(*self.image, max=max0, min=min0)
;
		return, im
;
		END

;----------------------------------------------------------------------------
; FullDiskObs::ExtractHelio
;
; Purpose:
;
;

FUNCTION FullDiskObs::ExtractHelio, lon0, lon1, lat0, lat1, res, soho=soho

	nx=fix((lon1-lon0)/res+.5)
	ny=fix((lat1-lat0)/res+.5)
	imx=intarr(nx, ny)


	return, imx
;
END

;================================================================

PRO FullDiskObs::Cleanup


ptr_free, self.image
;stop
;ptr_free, self.locs

END


;----------------------------------------------------------------------------
; FullDiskObservation__DEFINE
;
; Purpose:
;  Defines the object structure for an orb object.
;
PRO FullDiskObs__define

    struct = { FullDiskObs,	$
    			inherits FDObs, $
    			image:ptr_new(/allocate_heap),	$
    			bitpix:0	$

             }
END