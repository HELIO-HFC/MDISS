;========================================================

pro RegProcess, image, ploc, mom, qsint, imageResult

;	PROCEDURE to analyze/find sunspots locally
;
;	used by detspot7morph.pro and wldetspot.pro (?)

;	INPUT
;		image		input image
;		ploc		(array of (?)) region locations
;		mom			region statistical data
;		qsint		quiet sun intensity
;	OUTPUT
;		imageResult

;get image size
info=size(image)
nx=info[1]
ny=info[2]

;get number of regions and deduce pixel coordinates
N=n_elements(ploc)
xpos=ploc mod nx
ypos=ploc / nx

xmin=min(xpos)-5
xmax=max(xpos)+5
ymin=min(ypos)-5
ymax=max(ypos)+5

;create buffer iimages and import the region
imtmp=intarr(nx, ny)
bintmp=bytarr(nx, ny)
imtmp[ploc]=image[ploc]
bintmp[ploc]=1

;crop the region and display cropped and buffer images
imCrop=imtmp[xmin:xmax, ymin:ymax]

Sc=5
window, 0, xs=1024, ys=1024
tvscl, hist_equal(imtmp, percent=.2)

window, 1, xs=Sc*(xmax-xmin+1), ys=Sc*(ymax-ymin+1), xp=0, yp=0
tvscl, congrid(hist_equal(imCrop, percent=.5), Sc*(xmax-xmin+1), Sc*(ymax-ymin+1))


LRim=label_region(median(imtmp lt 0.975*mom[0], 2))
hLRim=histogram(LRim, REVERSE_INDICES=r)


nr=n_elements(hLRim)


;result0=dialog_message([' Mean:'+string(mom[0])+' Variance:'+$
;							string(mom[1]), 'Skewness:'+string(mom[2])+', Kurtosis:'+$
;							string(mom[3])+', Mean Abs Dev:'+string(mom[4]), $
;							'Is there a Spot?', 'Number of Seeds/Spots in the Region:'+string(nr-2)], /question)

for i=2, nr-1 do begin
	locs=r(r[i]:r[i+1]-1)
	mnl=min(image[locs])
	mnlpos=where(image[locs] eq mnl)
	seed=locs[mnlpos[0]]
	xp=seed mod nx
	yp=seed / nx
	if 0.92*qsint gt image[xp, yp] then upperThr=0.92*qsint $
		else begin
				cf=(fix(image[xp, yp]/mom[0]*100)+1)/100.
				upperThr=cf*mom[0]
		endelse
 	region=search2d(image, xp, yp, 0, upperThr)
;	region=search2d(image, xp, yp, 0, 0.995*mom[0])

	if where(region eq 0) eq -1 and mom[2] lt 0 and $
				mom[4] gt 10 then begin



; ***** now check that the regions should not overgrow the neighbourhood
		btmp=bytarr(nx, ny)
		btmp[Region]=1
		tmp=bintmp+btmp
		Region=where(tmp eq 2)

		imageResult[Region]=1
		imtmp2=image
		imtmp2[region]=image[Region]/2
;		window, 6, xs=1024, ys=1024, title='Region Grown'
;		tvscl, imtmp2
;		wait, 1
	end

endfor


end