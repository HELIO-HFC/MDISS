function getquietsunint, image, nn=nn


MINVAL=100


if not keyword_set(nn) then nn=1000
hgrm=histogram(Image, nbin=nn, loca=locs)
 nh=n_elements(hgrm)

xx=where(locs gt MINVAL)
if xx(0) eq -1 then begin
	print, '************THIS SHOULDNT HAPPEN, PLEASE INVESTIGATE'
	return, -1
end
mx=max(hgrm(xx), ind)
qsint=locs(xx(ind(0)))

;if (nh eq 1) then return, -1
; mx=max(hgrm[1:nh-1])
; qsint0=where(hgrm[1:nh-1] eq mx)
; qsint=qsint0[0]+1
; print, 'Quiet Sun Intensity For The Image:', qsint

;stop
return, qsint
end