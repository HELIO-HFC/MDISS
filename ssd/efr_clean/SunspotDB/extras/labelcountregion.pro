;===================================================

pro LabelCountRegion, binary, nreg, plocs, image, moms=moms0

LRIM=label_region(binary)

Lhgrm=histogram(temporary(LRIm), REVERSE_INDICES=r)
n=n_elements(temporary(lhgrm))

if n eq 1 then begin
	nreg=0
	return
endif
nreg=n-1
;print, 'Number of Regions:', n-1
plocs=ptrarr(n-1, /allocate_heap)
if arg_present(moms0) then begin
	moms=dblarr(n-1, 5)
	FOR i=1, n-1 DO BEGIN

	   *plocs[i-1] = r(r[i]:r[i+1]-1)
	   	 region=image[*plocs[i-1]]
	   	 if n_elements(*plocs[i-1]) lt 2 then begin
	   	 		pl=*plocs[i-1]
				moms[i-1, 4]=0
			   	 moms[i-1, 0]=image[pl[0]]
			   	 moms[i-1, 1]=0
			   	 moms[i-1, 2]=0
			   	 moms[i-1, 3]=0
			endif else begin
			   	 mom=moment(region, mdev=md )
			   	 moms[i-1, 4]=md
			   	 moms[i-1, 0]=mom[0]
			   	 moms[i-1, 1]=sqrt(mom[1])
			   	 moms[i-1, 2]=mom[2]
			   	 moms[i-1, 3]=mom[3]
			endelse
			moms0=moms
	endfor
endif else FOR i=1, n-1 DO *plocs[i-1] = r(r[i]:r[i+1]-1)

end


;===================================================

;			pro LabelCountRegion, map, binary, nreg, plocs, moms=moms0
;
;			LRIM=label_region(binary)
;
;			Lhgrm=histogram(LRIm, REVERSE_INDICES=r)
;			n=n_elements(lhgrm)
;
;			nreg=n-1
;			;print, 'Number of Regions:', n-1
;			plocs=ptrarr(n-1, /allocate_heap)
;			if arg_present(moms0) then begin
;				moms=dblarr(n-1, 5)
;				FOR i=1, n-1 DO BEGIN
;
;				   *plocs[i-1] = r(r[i]:r[i+1]-1)
;				   	 region=map.image[*plocs[i-1]]
;				   	 if n_elements(*plocs[i-1]) lt 2 then begin
;				   	 		pl=*plocs[i-1]
;							moms[i-1, 4]=0
;						   	 moms[i-1, 0]=map.image[pl[0]]
;						   	 moms[i-1, 1]=0
;						   	 moms[i-1, 2]=0
;						   	 moms[i-1, 3]=0
;						endif else begin
;						   	 mom=moment(region, mdev=md )
;						   	 moms[i-1, 4]=md
;						   	 moms[i-1, 0]=mom[0]
;						   	 moms[i-1, 1]=sqrt(mom[1])
;						   	 moms[i-1, 2]=mom[2]
;						   	 moms[i-1, 3]=mom[3]
;						endelse
;						moms0=moms
;				endfor
;			endif else FOR i=1, n-1 DO *plocs[i-1] = r(r[i]:r[i+1]-1)
;
;			end