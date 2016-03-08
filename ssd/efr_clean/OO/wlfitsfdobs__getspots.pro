PRO WLFitsFDObs::GetSpots, Sparr, thr=thr, cent=cent, sbl=sbl, $
		nogs=nogs, no_show=no_show, one=one, res=res


;	KEYWORDS/PARAMETERS
;			res	-		resolution (degree per pixel) of
;						the heliographic image
;			no_show -	suppress all windows
;			nogs -		omit gaussian smoothing during detection
;
;			one	-		if set use the first detection method (check
;						(wl_morphological pro)

self->GetFileName, filename=fn, ext=ext, path=folder
self->setQuietSunInt
image=*self.image

;spots=WL_DetSpGS(image, self.nx, self.ny, self.xc, self.yc, self.R, self.QuietSunInt, threshold=10, cent=.5)

if not keyword_set(thr) then thr=10
if not keyword_set(cent) then cent=.1
if not keyword_set(res) then res=.1

;spots=WL_DetSpGS(image, self.nx, self.ny, self.xc, self.yc, self.R, self.QuietSunInt, threshold=thr, cent=cent,  /sbl)

;	spots=WL_DetSp(image, self.nx, self.ny, self.xc, self.yc, self.R, self.QuietSunInt)

spots=WL_DetSpGS3(image, self.nx, self.ny, self.xc, self.yc, self.R, self.QuietSunInt, $
				threshold=thr, no_show=no_show, $
				cent=cent,  /sbl, one=one)

;spots=readfits(folder+fn+'_detected'+ext)

im=self->getimage()
sp2=spots ge 1



		LabelCountRegion, sp2, n, ploc
;		print, 'Number of Sunspots:', n
;
		Sparr=objarr(n)

		for i=0, n-1 do begin
			umbra=where(spots(*ploc[i]) eq 2)
			umbraploc=ptr_new(umbra)
			if umbra[0] eq -1 then nu=0 else begin
				pl=*ploc[i]
				temp=bytarr(self.nx, self.ny)
				temp[pl[umbra]]=1
				LabelCountRegion, temp, nu, umbraploc
			;	print, nu
			endelse
			mean0=mean(im[*ploc[i]])
			Sparr[i]=obj_new('Sunspot', self, *ploc[i], umbra, nu, mean0, res=res)
		endfor
END
