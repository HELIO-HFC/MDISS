FUNCTION FitsFDObs::Init, File=File, location=location, map=map

;if not keyword_set(File) then $
	if keyword_set(location) then File=obj_new('ftsfile', location=location, map=map) $
		else File=obj_new('ftsfile', map=map)

		File -> GetProperty, header=header, image=image , location=location


    	self.header=ptr_new(header)
    	self.image=ptr_new(image)
    	self.location=location

		st=File->GetObs()
		struct_assign, st, self, /nozero

		obj_destroy, File
		self.bitpix=fxpar(header, 'BITPIX')

	return, 1
END

;---------------------------------------------------------------------------

PRO FitsFDObs::GetHeader, head, filename=filename

head=*self.header
filename=self.location

END

;----------------------------------------------------------------------------
; FILE::GETFileName, Path, Extension
;
; Purpose:
;  Retrieves the value of properties associated with the orb object.
;
PRO FitsFDObs::GetFilename, FILENAME=filename, PATH=path, EXT=ext

	fn=self.location
    k0=5
	ext='.fits'
	i1=strpos(fn, ext, /REVERSE_SEARCH)
	if i1 eq -1 then begin
	ext='.fts'
	i1=strpos(fn, ext, /REVERSE_SEARCH)
	if i1 eq -1 then fname=-1 else k0=4

	endif


i2=strpos(fn, '\', /REVERSE_SEARCH)
if i2 eq -1 then i2=strpos(fn, '/', /REVERSE_SEARCH)
if i2 eq -1 then fname=-1

filename= strmid(fn, i2+1, i1-i2-1)
path= strmid(fn, 0, i2+1)

END

;-------------------------------------------------------------------------



; PRO FitsFDObs::Cleanup

;stop
; ptr_free, self.header
; ptr_free, self.image

; END


;----------------------------------------------------------------------------
; FitsFDObs__DEFINE
;
; Purpose:
;  Defines the object structure for an orb object.
;
PRO FitsFDObs__define

    struct = { FitsFDObs, $
    			location:'', 	$
    			header:ptr_new(/allocate_heap), 	$
    			inherits FullDiskObs $
             }
END
