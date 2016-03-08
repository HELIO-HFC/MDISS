FUNCTION ftsfile::Init, LOCATION=location, IMAGE=image, HEADER=header, map=map

if not keyword_set(location) then $
	location=dialog_pickfile(path='/home/sz/IDLWorkspace/data/mdi/')
    self.location=location

if not keyword_set(image) or keyword_set(header) then $
			mreadfits, location, index, data, header=header
			index2map, index, data, map
	;	image=readfits(location, header)
    self.image=ptr_new(data)
    self.header=ptr_new(header)

	return, 1
END



;------------------------------------------------------------------------------
;----------------------------------------------------------------------------
; FILE::GETFileName, Path, Extension
;
; Purpose:
;  Retrieves the value of properties associated with the orb object.
;
PRO ftsfile::GetFilename, FILENAME=filename, PATH=path, EXT=ext

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


;================================================================================

FUNCTION ftsfile::GetObs


  ;  	self.header=ptr_new(header)
   ; 	self.image=ptr_new(image)
    ;	self.location=location
		GetHeaderInfo, *self.header, nx, ny, date, xc, yc, R, cdx, cdy


orkey= fxpar(*self.header, 'ORIGIN')
if orkey eq 'SOI Science Support Center' then origin='soho' else origin='ground'
RES={ FDObs, $
		R:R,		$
		xc:xc,			$
		yc:yc, 			$
		cdx:cdx,		$
		cdy:cdy,		$
		nx:nx,		$
		ny:ny,		$
		date:date, 		$
		origin:origin	$
		}


;struct_assign, res1, res
	return, RES
END

;================================================================

PRO ftsfile::Cleanup

ptr_free, self.header
ptr_free, self.image
;stop
;ptr_free, self.locs

END

;----------------------------------------------------------------------------
; File__DEFINE
;
; Purpose:
;  Defines the object structure for an orb object.
;
PRO ftsfile__define
    struct = { ftsfile, $
    			location:'', $
               	image:ptr_new(/allocate_heap), $
               	header:ptr_new(/allocate_heap) $
             }
END
