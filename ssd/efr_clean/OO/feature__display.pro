;========================================================

PRO Feature::Display, FullDiskObs, maxint=maxint, window=window
if arg_present(FullDiskObs) then begin

		image=FullDiskObs->getdisplay()
	info=size(image)
if info[1] eq self.nx and info[2] eq self.ny then begin
		bd=image
		if keyword_set(maxint) then	bd[*self.locs]=max(bd) $
			else bd[*self.locs]=0
	endif else begin
		bd=bytarr(self.nx, self.ny)
		bd[*self.locs]=1
	endelse
endif else begin
	bd=bytarr(self.nx, self.ny)
	bd[*self.locs]=1
     endelse
if keyword_set(window) then window, 0, xs=self.nx, ys=self.ny
tvscl, bd
;im_out=bd
;return, bd
END
