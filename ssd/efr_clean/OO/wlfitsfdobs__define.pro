function WLFitsFDObs::Init, inputfile=inputfile, location=location, map=map, filename=filename, $
                            unprocessed=unprocessed, header=header 

NULLVALUE=0
;if not keyword_set(File) then $
  if keyword_set(filename) then location=filename
	if keyword_set(location) then File=obj_new('ftsfile', location=location, map=map) $
				else File=obj_new('ftsfile', map=map)

		File -> GetProperty, header=header, image=image , location=location

    xxlocs=where(finite(image), comp=clocs)
    if clocs[0] ne -1 then image(clocs)=NULLVALUE

    	self.header=ptr_new(header)
    	self.image=ptr_new(image)
    	self.location=location

		st=File->GetObs()
		struct_assign, st, self, /nozero

		

	;	correct SOLAR_P angle if necessary

	;self->getheader, hd
	solar_p=fxpar(header, 'SOLAR_P0')

if keyword_set(unprocessed) then begin
   print, 'pre-processing.....'
   im=*self.image
   locs=where(im lt 0)
   if locs[0] ne -1 then im[locs]=0
   Qsim=qsmedian(im, self.r, self.xc, self.yc)
   result=getcontrastimage(im, qsim, self.r, self.xc, self.yc, /subtract)
 ;  stop
   *self.image=result
   print, 'done'
endif

self.bitpix=fxpar(header, 'BITPIX')
		self->setQuietSunInt

	if (solar_p ne 0) then begin

		im=*self.image
		im=rot(temporary(im), solar_p, 1, self.xc, self.yc)
		*self.image=im
	 endif

; *** cleaning bit added 2005-07-01
obj_destroy, file

	;	STOP
	return, 1
END




;-------------------------------------------------------------------------



PRO WLFitsFDObs::Cleanup

;stop
ptr_free, self.header
ptr_free, self.image

END

;----------------------------------------------------------------------------
; FitsFDObs__DEFINE
;
; Purpose:
;  Defines the object structure for an orb object.
;
PRO WLFitsFDObs__define

    struct = { WLFitsFDObs, $
    			QuietSunInt:0.d00, $
    			inherits FitsFDObs $
    		   }
END
