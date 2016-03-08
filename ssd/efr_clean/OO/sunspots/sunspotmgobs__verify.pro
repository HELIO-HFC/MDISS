
;-------------------------------------------------------------

FUNCTION SunspotMgObs::Verify

; Verifying Sunspot Detection results
; against magnetogram data

image=bytarr(self.nx, self.ny)

n=n_elements(*self.spotarr)
mgim=*self.mgdata
;print, n
for i=0, n-1 do begin

	locs=(*self.spotarr)[i]->feature::getlocs()
	ulocs=(*self.spotarr)[i]->getumbra()
	(*self.spotarr)[i]->getsdata, arcx, arcy
;	print, max(mgim[locs]), min(mgim[locs])
	mxflux=max(abs(mgim[locs]))
	if ( mxflux gt 130) $
;		or (mxflux gt 100 and abs(ArcX) gt 900) $
			then begin
		image[locs]=1
		if ulocs[0] ne -1 then image[locs[ulocs]]=2
	endif else begin

			continue
			window, 10, xs=self.nx, ys=self.ny
		    tvscl, *self.image

		    window, 11, xs=self.nx, ys=self.ny
		    tvscl, bytscl(*self.mgdata, max=1500, min=-1500)

			bim=bytarr(self.nx, self.ny)
			bim[locs]=1
			window, 12, xs=self.nx, ys=self.ny
				tvscl, bim

		;	stop
		endelse
endfor



return, image

END
