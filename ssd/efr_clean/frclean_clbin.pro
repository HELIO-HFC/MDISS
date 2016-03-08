;===========================================================================

function frclean_clbin, binary, xc, yc, r1, r2

;	FUNCTION
;		cleaning the input binary image by setting all pixels
;		outside r1 and r2 to 0
;	INPUT
;		image
;		r1, r2	cleaning annuli radiuses
;		xc, yc	centre (must be float)

xc=float(xc)
yc=float(yc)

rmax=float(r1>r2)
rmin=float(r1<r2)

info=size(binary)
nx=info[1]
ny=info[2]
im=binary
locs=where(binary)

n=n_elements(locs)
if n lt 2 then return, binary
count=0
for k=0l, n-1, 1 do begin
	j=fix(locs[k]/nx)
	i=locs[k] mod nx
	dist=sqrt((xc-i)*(xc-i)+(yc-j)*(yc-j))
	if dist gt rmax or dist lt rmin then begin
		 im[locs[k]]=0
		 count=count+1
	endif
endfor

;print, 'Number of Points Removed:', count
return, im


end


