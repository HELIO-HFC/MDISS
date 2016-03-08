function MDI_AR, image, header, display=display, DF=DF, rect=rect


if not arg_present(image) then begin
	fn=dialog_pickfile(path='Y:\mdi_data\mg_april_02')
	image=readfits(fn, header)
endif

GetHeaderInfo, header, nx, ny, date, xc, yc, R, cdx, cdy

;info=size(image)
;nx=info[1]
;ny=info[2]

;im=image
image=double(image)
im=gauss_smoothing(image, 3)

bim=bytscl(im, max=1500, min=-1500)



impos=im gt 70
imneg=bytarr(nx, ny)

locs=where(im gt -30000 and im lt -70)
imneg[locs]=1

imneg=RemLimbBinary (imneg, xc, yc, R, nx)

;window, 1, xs=1024, ys=1024, title='Positive Thresh'
;tvscl, impos
;window, 2, xs=1024, ys=1024, title='Negative Thresh'
;tvscl, imneg

if keyword_set(display) then begin
	window, 0, xs=1024, ys=1024, title='Input Smoothed'
	tvscl, bim
	window, 3, xs=nx, ys=ny
	tvscl, 2*impos+imneg
endif
;window, 4, xs=nx, ys=ny
;tvscl, impos+imneg

if not keyword_set(DF) then DF=5

if not keyword_set(rect) then $
	circ=shift(dist(2*df+1), df, df) le df $
	else circ=replicate(1, df, df)

impos1=dilate(impos, circ)
imneg1=dilate(imneg, circ)

ComAdd=impos1+imneg1



locs=where(ComAdd gt 1)

help, locs

ComAddx=intarr(nx, ny)
ComAddx[locs]=1

ARpr=intarr(nx, ny)

LR=label_region(ComAddx)

hLR=histogram(LR, reverse_indices=r)

n=n_elements(hLR)

print, 'There are possibly ', n-1, ' Active Regions in this picture'

FOR i=1, n-1 DO BEGIN

   seedloc = r(r[i]:r[i+1]-1)
	xp=seedloc[0] mod nx
	yp=seedloc[0] / nx

	lcs=where(ARpr[seedloc] eq 1)
	if lcs[0] eq -1 then begin
		reg=search2d(ComAdd, xp, yp, 1, 2)
		ARpr[reg]=1
	endif
endfor
if keyword_set(display) then begin
window, 5, xs=nx, ys=ny
tvscl, ComAdd
window, 6, xs=nx, ys=ny
tvscl, ARpr
endif
return, ARpr
end