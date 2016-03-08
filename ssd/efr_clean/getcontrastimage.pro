;===================================================================================================

function getcontrastimage, image, QSImage, R, xc, yc, Inorm=Inorm, subtract=subtract


;	FUNCTION
;			Producing cleaned image by dividing/subtracting original by quiet sun
;
;		INPUT
;				image 			input image
;				QSImage			Quiet Sun Image
;				R				Disk Radiuse
;				xc, yc			disk centre coordinate
;	 (optional)	Inorm			normalizing parameter for division method
;
;		KEYWORD
;				subtract		cleaned image is taken as difference between the original
;									and quiet sun
;
;		RESULT
;				cleaned image


info=size(image)
nx=info[1]
ny=info[2]
CImage=QSImage
if not keyword_set(Inorm) then Inorm=max(QSImage)

if not keyword_set(subtract) then begin
		for i=0, nx-1 do $
			for j=0, ny-1 do begin
				IF (CImage(i,j)  GT 0) THEN CImage(i,j) $
				 = FIX(DOUBLE(image(i,j)) / DOUBLE(CImage(i,j))   * Inorm +0.5)
			endfor
			endif else begin
				locs=where(QSImage eq 0)
				CImage=image-QSImage
				mci=min(CImage)
				xarr=replicate(mci, nx, ny)

				BLAS_AXPY, CImage, -1, xarr
				CImage[locs]=0
			endelse

return, CImage
end


