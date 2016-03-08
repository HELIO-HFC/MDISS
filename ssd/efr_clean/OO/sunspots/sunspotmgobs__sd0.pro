;-------------------------------------------------------------

PRO SunspotMgObs::SD,  no_show=no_show, cent=cent, one=one,$
						res=res, thr=thr

image=*self.image
;stop

if not keyword_set(thr) then thr=10
if not keyword_set(cent) then cent=.1
if not keyword_set(res) then res=.1
spots=WL_DetSpGS3(image, self.nx, self.ny, self.xc, self.yc, self.R, self.QuietSunInt, $
				threshold=thr, no_show=no_show, $
				cent=cent,  /sbl, one=one)


;mgim=*self.mgdata

spot_image=bytarr(self.nx, self.ny)
mgim=gauss_smoothing(DOUBLE(*self.mgdata), 2)
sp2=spots ge 1



		LabelCountRegion, sp2, n, ploc
		count=0
print, n, '  candidates'
		Sparr=objarr(n)

		for i=0, n-1 do begin


			if abs(max(mgim[*ploc[i]])) > abs(min(mgim[*ploc[i]])) lt 100 then continue
			spot_image[*ploc[i]]=1
			umbra=where(spots(*ploc[i]) eq 2)
			umbraploc=ptr_new(umbra)
			if umbra[0] eq -1 then nu=0 else begin
				pl=*ploc[i]
				temp=bytarr(self.nx, self.ny)
				temp[pl[umbra]]=1
				spot_image[pl[umbra]]=2
				LabelCountRegion, temp, nu, umbraploc
			;	print, nu
			endelse
			mean0=mean(image[*ploc[i]])
			Sparr[count]=obj_new('Sunspot', self, *ploc[i], umbra, nu, mean0, res=res)
	;		print, count, max(mgim[*ploc[i]]), min(mgim[*ploc[i]])
			count=count+1
		endfor
if not keyword_set(no_show) then begin

	window, 24, xs=self.nx, ys=self.ny
	tvscl, spot_image
	print, count, '   spots'
end

self.n_spots=count
if (count ne 0) then $
	sp=sparr[0:count-1] $
  else sp=0

self.spotArr=ptr_new(sp)

;stop
END
