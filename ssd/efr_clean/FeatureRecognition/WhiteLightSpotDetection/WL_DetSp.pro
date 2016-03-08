function WL_DetSp, image, nx, ny, xc, yc, Radius, qsint, no_show=no_show

;	PROCEDURE
;		detecting and analyzing/classifying sunspots for MDI
;		white light dopplergrums combined with magnetograms

;if not keyword_set(fn) then WL_ReadInData, image, nx, ny, xc, yc, Radius, qsint, imdi, bimdi, wls, mdate $
; else WL_ReadInData, image, nx, ny, xc, yc, Radius, qsint, imdi, bimdi, wls, mdate, fn=fn
binary=WL_Morphological(image, xc, yc, Radius)

WL_AnalyseBinary, image, binary, n, plocs, moms

imageResult=intarr(nx, ny)
imageResult2=intarr(nx, ny)

for i=0, n-1 do begin
	WL_RegProcess, image, *plocs[i], moms[i, *], qsint, imageResult, imageResult2
endfor


lres=where(imageResult ne 0)
lres2=where(imageResult2 ne 0)

binary=bytarr(nx, ny)
binary(lres)=1

imageumbra=WL_DetectUmbra (image, binary, imdi, bimdi, nx, ny, qsint)

image_copy=image
image_copy[lres]=0;max(image)
;image_copy[lres]=image[lres]*2

image_copy2=image
image_copy2[lres2]=image[lres2]*2



;print, 'Sunspot Area', n_elements(lres), '  pixels'
;print, 'Sunspot Area to Full Disk Ratio:', double(n_elements(lres))/(2*!pi*Radius*Radius)*100, ' percent'


;bimdiCopy=bimdi
;bimdiCopy[lres]=255b

if not keyword_set(no_show) then begin
window, 13, xs=nx, ys=ny, title='Detected Sunspots'
tvscl, image_copy

;window, 21, xs=nx, ys=ny, title='Magnetogram with Sunspots superimposed'
;tvscl, bimdiCopy

window, 23, xs=nx, ys=ny, title='Umbras'
tvscl, imageumbra+binary


;stop
endif
detected=imageumbra+binary


return, detected
end

