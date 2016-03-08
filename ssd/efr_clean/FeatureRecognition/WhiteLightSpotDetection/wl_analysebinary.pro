pro WL_AnalyseBinary, image, binary, nreg, plocs, moms

;	PROCEDURE
;		calculates the locations of the regions from binary and
;		deduces regions' statistical properties wrt image
;
;		if region consists of a single pixel, the mean is set equal
;		to the intensity of the pixel

locs=where(binary ne 0)


LRIM=label_region(binary)

Lhgrm=histogram(LRIm, REVERSE_INDICES=r)
n=n_elements(lhgrm)

nreg=n-1
;print, 'Number of Regions:', n-1
moms=dblarr(n-1, 5)
plocs=ptrarr(n-1, /allocate_heap)

FOR i=1, n-1 DO BEGIN

   *plocs[i-1] = r(r[i]:r[i+1]-1)
   	 region=image[*plocs[i-1]]

	if n_elements(*plocs[i-1]) gt 2 then begin
	   	 mom=moment(region, mdev=md )
	   	 moms[i-1, 4]=md
	   	 moms[i-1, 0]=mom[0]
	   	 moms[i-1, 1]=sqrt(mom[1])
	   	 moms[i-1, 2]=mom[2]
	   	 moms[i-1, 3]=mom[3]
	 end else begin
	 ;	help, *plocs[i-1]
	  if n_elements(*plocs[i-1]) eq 1 then $
	 	moms[i-1, 0]=image[*plocs[i-1]] $
	 		else moms[i-1, 0]=total(image[*plocs[i-1]])/2.
	 endelse

endfor

end