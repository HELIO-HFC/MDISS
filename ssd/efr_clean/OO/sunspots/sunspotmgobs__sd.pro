;-------------------------------------------------------------

PRO SunspotMgObs::SD,  no_show=no_show, cent=cent, one=one, display=display, $
						res=res, thr=thr, missvalscheck=missvalscheck, checkmg=checkmg


;
;	KEYWORDS
;			one			use "overdetect" option for WL detection
;						preferred option when using magnetogram data
;						to validate detection results

image=*self.image
;stop

if not keyword_set(thr) then thr=10
if not keyword_set(cent) then cent=.1
if not keyword_set(res) then res=.1
spots=WL_DetSpGS3(image, self.nx, self.ny, self.xc, self.yc, self.R, self.QuietSunInt, $
				threshold=thr, no_show=no_show, $
				cent=cent,  /sbl, one=one)


MAXIMUM_FEATURE_LENGTH_X=300
MAXIMUM_FEATURE_LENGTH_Y=300


;stop
spot_image=bytarr(self.nx, self.ny)
if keyword_set(missvalscheck) then mgim=*self.mgdata $
	else mgim=*self.mgdata;mgim=gauss_smoothing(DOUBLE(*self.mgdata), 5)
sp2=spots ge 1



		LabelCountRegion, sp2, n, ploc
		count=0
;print, n, '  candidates'
		if n ne 0 then	Sparr=objarr(n) else begin
			self.n_spots=n
			return
		 endelse

  for i=0, n-1 do begin

			maxfluxvalue=abs(max(mgim[*ploc[i]])) > abs(min(mgim[*ploc[i]]))

			locs=*ploc[i]
			irrad=total(image[locs])/(self.QuietSunInt*n_elements(locs))

                        minir=min(image[locs])




			if keyword_set(checkmg) then begin

				print, mgim[locs]
				im2=image
				im2[locs]=2*im2[locs]

				mg2=mgim
				mg2[locs]=1000

				print, irrad
				print, maxfluxvalue
				window, 1, xs=800, ys=600	;xs=1024, ys=1024
				tvframe, bytscl(mg2, max=1000, min=-1000), /asp
			;	wset, 2
				stop
      endif


                            xp=locs mod self.nx
                            yp=locs mod self.ny

                            if (max(xp) - min(xp)) ge MAXIMUM_FEATURE_LENGTH_X then continue
                            if (max(yp) - min(yp)) ge MAXIMUM_FEATURE_LENGTH_Y then continue

			if ((n_elements(locs) le 2) and irrad gt .98) then continue
			if (irrad lt .85)  then begin
							if (n_elements(locs) lt 10 and $
								 (maxfluxvalue lt 75)) then continue $
								else if (maxfluxvalue lt 40) then continue
					endif $
				else if (maxfluxvalue lt 100) then continue

                                 
			if keyword_set(missvalscheck) then $
				if (maxfluxvalue gt 32000) then continue

                        if keyword_set(missvalscheck) then $
				if (minir lt -30000) then continue



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
			count=count+1

                        ptr_free, umbraploc
  endfor
if not keyword_set(no_show) then begin
	window, 24, xs=800, ys=600 ;xs=self.nx, ys=self.ny
	tvframe, spot_image, /asp
;	tvscl, spot_image
	;print, count, '   spots'
end
;print, count, '   spots'
self.n_spots=count
if (count ne 0) then $
	sp=sparr[0:count-1] $
  else sp=0

self.spotArr=ptr_new(sp)

;stop
END
