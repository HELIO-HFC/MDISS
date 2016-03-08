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

