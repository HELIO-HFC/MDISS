function WL_DetSpGS2, image, nx, ny, xc, yc, Radius, qsint, $
				no_show=no_show, threshold=threshold, cent=cent, $
				nogs=nogs, sbl=sbl, one=one

;	PROCEDURE
;		detecting and analyzing/classifying sunspots for MDI
;		white light dopplergrums combined with magnetograms
;
;
;	INPUT
;		image		a preprocessed (flat) 'white light' image
;					(SOHO/MDI continuum IntensityGram)
;		nx, ny		image size parameters
;		xc, yc		Solar disk centre pixel coordinates
;		Radius		Solar disk radius
;
;	KEYWORDS
;
;		/no_show	do not display results
;		/nogs		do not use gaussian smoothing prior to edge
;					detection
;		/sbl		use Sobel edge detection
;					(as opposed to morphological gradient)
;
;	OPTIONAL PARAMETERS
;
;		threshold
;		cent


;if not keyword_set(fn) then WL_ReadInData, image, nx, ny, xc, yc, Radius, qsint, imdi, bimdi, wls, mdate $
; else WL_ReadInData, image, nx, ny, xc, yc, Radius, qsint, imdi, bimdi, wls, mdate, fn=fn
image1=image

if not keyword_set(nogs) then image=gauss_smoothing(image, 2)
;STOP

if not keyword_set(sbl) then $
	binary=WL_Morphological(image, xc, yc, Radius, threshold=threshold, cent=cent) $
		else binary=WL_Morphological2(image, xc, yc, Radius, $
			threshold=threshold, cent=cent,  /sbl, edges=edges)




WL_AnalyseBinary, image, binary, n, plocs, moms

binx=binary

imageResult=intarr(nx, ny)
imageResult2=intarr(nx, ny)

for i=0, n-1 do begin
	WL_RegProcess3, image, *plocs[i], moms[i, *], qsint, imageResult, imageResult2
endfor


lres=where(imageResult ne 0)
lres2=where(imageResult2 ne 0)

binary=bytarr(nx, ny)
if keyword_set(one) then binary(lres)=1 else binary(lres2)=1



imageumbra=WL_DetectUmbra (image, binary, imdi, bimdi, nx, ny, qsint)

lumb=where(imageumbra ne 0)
image_copy=image
;image_copy[lres]=0;max(image)
image_copy[lres]=image[lres]*2

image_copy2=image
image_copy2[lres2]=image[lres2]*2
image_copy2[lres2]=max(image)
if lumb[0] ne -1 then image_copy2[lumb]=min(image)



if not keyword_set(no_show) then begin

		window, 13, xs=nx, ys=ny;, title='Detected Sunspots'
		tvscl, image_copy2

		window, 21, xs=nx, ys=ny;, title='Magnetogram with Sunspots superimposed'
		tvscl, image_copy

		xxx=imageumbra+binary

		N=100
			for i=0, N-1 do begin
				xi=xc+Radius*cos((2*!pi/N)*i)
				yi=yc+Radius*sin((2*!pi/N)*i)

				xxx[fix(xi+.5), fix(yi+.5)]=max(xxx)
			end
		window, 23, xs=nx, ys=ny;, title='Umbras'
		tvscl, xxx


;stop
endif
detected=imageumbra+binary

image=image1
return, detected
end

