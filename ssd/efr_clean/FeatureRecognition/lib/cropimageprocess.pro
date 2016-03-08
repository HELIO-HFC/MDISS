pro CropImageProcess, Cimage, image, xmin, ymin, qsint, mom, imageResult, imageResult2

;	PROCEDURE
;			processing cropped sunspot regions
;
;	INPUT
;		Cimage		cropped image
;		image		original image
;		xmin, ymin	cropped image lower left corner coord in original pixels
;		qsint		quiet sun intensity
;		mom			region's statistical moment array

info=size(image)
nx=info[1]
ny=info[2]
info=size(Cimage)


nxc=info[1]
nyc=info[2]

imCrop=Cimage
if nxc le 2 or nyc le 2 then begin
								Region=where(imcrop ge 0 and imcrop lt 0.91*qsint)
								xloc=Region mod nxc
								yloc=Region /nxc
								imageResult[xloc+xmin, yloc+ymin]=1
								return
							endif


Sc=10

;window, 1, xs=Sc*nxc, ys=Sc*nyc, xp=0, yp=0
;tvscl, congrid(imCrop, Sc*nxc, Sc*nyc)




;examine the cropped image

;print, 'Quiet Sun Intensity:', qsint
;print, 'Mean Intensity:',mom[0]
;print, 'St. Deviation:', mom[1]
;print, 'Mean Abs Deviation:', mom[4]

; threshold the image at mean intensity and label the regions
locs=where(imCrop ge 0 and imCrop lt mom[0]-mom[4]/4.)
imc1=bytarr(nxc, nyc)
imc1[locs]=1


;window, 2, xs=Sc*nxc, ys=Sc*nyc, xp=Sc*nxc, yp=Sc*nyc
;tvscl, congrid(imc1, Sc*nxc, Sc*nyc)

LRimc1=label_region(imc1)
hlrimc1=histogram(LRimc1, reverse_indices=r)
nh=n_elements(hlrimc1)

imc2=intarr(nxc, nyc)
imc3=intarr(nxc, nyc)
;print, 'Number of Candidate Sunspots in the Region:', nh-1
for i=1, nh-1 do begin
	plcs=r(r[i]:r[i+1]-1)

	;calculate seed location
	mnlc=where(imCrop[plcs] eq min(imCrop[plcs]))

	xp=plcs[mnlc[0]] mod nxc
	yp=plcs[mnlc[0]] / nxc

	;calculate Upper Threshold

	UpperThreshold=.79*Qsint > mom[0]
	UpperThreshold2=.93*Qsint> mom[0]
	Region=Search2d(imCrop, xp, yp, 0, UpperThreshold)
	Region1=Search2d(imCrop, xp, yp, 0, UpperThreshold2)

	imc2[Region]=imCrop[Region]
	imc3[Region1]=imCrop[Region1]

	xloc=Region mod nxc
	yloc=Region /nxc
	imageResult[xloc+xmin, yloc+ymin]=1

	xloc1=Region1 mod nxc
	yloc1=Region1 /nxc
	imageResult2[xloc1+xmin, yloc1+ymin]=1

endfor

;window, 3, xs=Sc*nxc, ys=Sc*nyc, xp=2*Sc*nxc, yp=2*Sc*nyc
;tvscl, congrid(imc2, Sc*nxc, Sc*nyc)
;stop
end