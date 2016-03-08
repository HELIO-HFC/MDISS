pro limb_efit, image, xc1, yc1, R1, R2, theta, stdev, stdevGeo, algerr, $
						e_specific=e_specific, verbose=verbose, binary, binary2

;		PROCEDURE
;					Fitting an Ellipse to the Full Disk Solar Image
;		INPUT:
;			image		the full disk solar image (assumed)
;		OUTPUT:
;			xc1, yc1 	ellipse centre pixel coordinates
;			R1, R2  	ellipse axis
;			theta		angle
;			stdev		standard deviation
;			stdevGeo	standard deviation geometric
;			algerr		algebraic error
;			binary		cleaned binary image for the last fit (to be removed)
;			binary2		initial binary image (to be removed?)
;
;
;		KEYWORDS
;			e_specific 	if set non_svd (ellipse specific) algorithm is used for fitting
;		(useful when testing):

;
;		KNOWN ISSUES:
;			not reliable when working with 8-bit data due to current implementation of Canny
;				( to be resolved soon)
;
;	 	LimbEllipseFit, image [,/e_specific [, xc1, xc2, theta, stdev, stdevGeo, algerr, binary, binary2]
;
;
;		CALLLS:
;			frclean_clim, frclean_gcentre, frclean_clbin, frclean_checkgap, frclean_dhistclean,
;			gauss_smoothing, canny, ellipsefit, ellipsefitsvd
;



info=size(image)
nx=info[1]
ny=info[2]
type=info[3]


; checking input image data type (see KNOWN ISSUES)
; if not byte or integer, then scaling to Byte
if type ne 1 and type ne 2 then image=bytscl(image)
image_copy=image


cenx=nx/2.
ceny=ny/2.
box_size=(nx < ny) /16

HistImage=histogram(image)
n=n_elements(HistImage)
AveHistImage=fltarr(n)

;	Generate a Smoothed Histogram

Ni=n/128

if Ni lt 4 then Ni=4
for i=0, Ni-1, 1 do AveHistImage[i]=HistImage[i]

for i=Ni, n-Ni-1, 1 do begin
	s=float(HistImage[i])
	for j=1, Ni, 1 do s=s+float(HistImage[i+j])+float(HistImage[i-j])
	AveHistImage[i]=s/(2.*Ni+1)
end

sum=0ul
count=0l

; Get average intensity value around the centre of the image

for i=fix(cenx+.5)-box_size, fix(cenx+.5)+box_size do $
	for j=fix(ceny+.5)-box_size, fix(ceny+.5)+box_size do begin
		sum=sum+image[i, j]
		count=count+1
	endfor

average=sum/count

ave_value=long(average+.5)

;	Analyse Histogram to get the cutoff intensity around limb

ind=fix(average+.5)
thr_value=0

value=AveHistImage[ind]

while ind gt 2 and (HistImage[ind-1] lt HistImage[ind] or HistImage[ind] gt 0.05*value) do ind=ind-1
thr_value=ind

if keyword_set(verbose) then begin
	print, 'Threshold Value:', thr_value
	print, 'Average Value:', ave_value
endif

; If threshold value found is equal to zero then stop
if thr_value lt 5 then begin
	RLT=dialog_message('Failed. Could not find threshold cutoff value', /error)
	retall
end

; Clean Image by LabelRegion
frclean_clim, image, thr_value, imageB, count

; (OR) Clean Image by rastering
;			CleanImage, image, thr_value, cenx, ceny, nx, ny, imageB, count



R=sqrt(double(count)/!pi)
; calculate the gravity centre of the cleaned image
grv_cen= frclean_gcentre (imageB)
if keyword_set(vervose) then begin
	print, 'Count:', count
	print, 'Radius:', R
	print, 'Gravity Centre:', grv_cen
endif
xc=grv_cen[0]
yc=grv_cen[1]

; produce the binary image suitable for the ellipse fit

if type eq 2 then begin
			CannyUpperThreshold=400
			CannyDelta=40
			CannyLowerThreshold=80
endif else begin
				CannyUpperThreshold=50
				CannyDelta=5
				CannyLowerThreshold=5
endelse

times=0
im=gaussiansmoothing(image_copy, 2l)

CannyUpperThreshold=.95
CannyLowerThreshold=.9
;stop
; Ellipse Fittiing Loop:
jump0:	binary=canny(im, high=CannyUpperThreshold, $
	low=CannyLowerThreshold)


		frclean_dhistclean, binary, xc, yc, R, R1, R2, flag
		times=times+1
		if flag eq 1 and times lt 10 then begin
			if keyword_set(vervose) then $
				print, 'Not Enough Points, Changing Threshold'
			CannyUpperThreshold=CannyUpperThreshold-40l
		;	print, times, canny_threshold
			goto, jump0
		endif
		if flag eq 1 then begin
			RLT=dialog_message('Failed. Could not determine Canny Threshold', /error)
			retall
		end
		binary=frclean_clbin(binary, xc, yc, R1, R2 )
		gap=frclean_checkgap(binary, xc, yc)
		if keyword_set(vervose) then $
			print, 'gap=', gap
		if gap gt 25 and times lt 10 then begin
			if keyword_set(vervose) then $
				print, 'Big Gap, Changing Threshold'
			CannyUpperThreshold=CannyUpperThreshold-40l
		;	print, times, canny_threshold
			goto, jump0
		endif

; fit the ellipse to the binary image using SVD fit procedure
binary2=binary


if keyword_set(e_specific) then $
	frclean_ellipsefit, binary2, xc, yc,  xc1, yc1, R1, R2, theta, stdev, stdevGeo, algerr $
	else frclean_ellipsefitsvd, binary2, xc, yc,  xc1, yc1, R1, R2, theta, stdev, stdevGeo, algerr
end










