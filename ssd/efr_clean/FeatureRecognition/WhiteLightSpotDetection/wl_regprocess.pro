;========================================================

pro WL_RegProcess, image, ploc, mom, qsint, imageResult, imageResult2, display=display

;	PROCEDURE to analyze/find sunspots locally
;
;	used by detspot7morph.pro and wldetspot.pro (?)

;	INPUT
;		image		input image
;		ploc		(array of (?)) region indices
;		mom			region statistical data
;		qsint		quiet sun intensity
;	OUTPUT
;		imageResult

;discard small regions (less then 5 pixels)
;if n_elements(ploc) lt 5 then return

;get image size
info=size(image)
nx=info[1]
ny=info[2]

;get number of regions and deduce pixel coordinates
N=n_elements(ploc)
xpos=ploc mod nx
ypos=ploc / nx

xmin=min(xpos)
xmax=max(xpos)
ymin=min(ypos)
ymax=max(ypos)

	;create buffer images and import the region
imtmp=make_array(nx, ny, /int, value=-1)
bintmp=bytarr(nx, ny)
imtmp[ploc]=image[ploc]
bintmp[ploc]=1

	;crop the region and display cropped and buffer images
imCrop=imtmp[xmin:xmax, ymin:ymax]


nxc=xmax-xmin+1
nyc=ymax-ymin+1

if nxc le 2 or nyc le 2 then begin
								Region=where(imcrop ge 0 and imcrop lt 0.91*qsint)
								xloc=Region mod nxc
								yloc=Region /nxc
								imageResult[xloc+xmin, yloc+ymin]=1

			endif else begin


				Sc=10
			if keyword_set(display) then begin
				window, 1, xs=Sc*nxc, ys=Sc*nyc, xp=0, yp=0
				tvscl, congrid(imCrop, Sc*nxc, Sc*nyc)

					;examine the cropped image

				print, 'Quiet Sun Intensity:', qsint
				print, 'Mean Intensity:',mom[0]
				print, 'St. Deviation:', mom[1]
				print, 'Mean Abs Deviation:', mom[4]
			end
					; threshold the image at mean intensity and label the regions
				locs=where(imCrop ge 0 and imCrop lt mom[0]-mom[4]/4.)
				imc1=bytarr(nxc, nyc)
				imc1[locs]=1

			if keyword_set(display) then begin
				window, 2, xs=Sc*nxc, ys=Sc*nyc, xp=Sc*nxc, yp=Sc*nyc
				tvscl, congrid(imc1, Sc*nxc, Sc*nyc)
			end
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

					UpperThreshold=.8*Qsint > mom[0]
					UpperThreshold2=.93*Qsint> mom[0]
					Region=Search2d(imCrop, xp, yp, 0, UpperThreshold)
					Region1=Search2d(imCrop, xp, yp, 0, UpperThreshold2)

		;			Region= where(imcrop le UpperThreshold)


					imc2[Region]=imCrop[Region]
					imc3[Region1]=imCrop[Region1]

					xloc=Region mod nxc
					yloc=Region /nxc
					imageResult[xloc+xmin, yloc+ymin]=1

					xloc1=Region1 mod nxc
					yloc1=Region1 /nxc
					imageResult2[xloc1+xmin, yloc1+ymin]=1

				endfor
		endelse

end