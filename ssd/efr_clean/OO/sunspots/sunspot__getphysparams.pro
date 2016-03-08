
;-------------------------------------------------------------

FUNCTION Sunspot::GetPhysParams, fd, min_value=min_value, $
					max_value=max_value, mdev=mdev, m2QS=m2QS, $
					print=print, direct_image=direct_image, umbra=umbra

; Verifying Sunspot Detection results
; against magnetogram data


if keyword_set(direct_image) then image=fd else $
		if not obj_isa(fd, 'fulldiskobs') then begin
				print, 'Incorrect object!'
				return, -1
			endif else	image=fd->getimage()
;stop

locs=*self.locs
if keyword_set(umbra) then begin
		umblocs=*self.umbra
		if umblocs[0] eq -1 then return, -1
		locs=locs[umblocs]
	endif


max_value=max(image[locs])
min_value=min(image[locs])
if n_elements(locs) ge 2 then $
	mom=moment(image[locs], mdev=mdev, /double) $
	else mom=image[locs]

m2QS=mom[0]/self.QuietSunInt
if keyword_set(print) then begin
	print, max_value, min_value $
	, mom[0], mdev, m2QS, self.QuietSunInt
endif
return, mom[0]

END