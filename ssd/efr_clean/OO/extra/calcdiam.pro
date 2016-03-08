;================================================

function CalcDiam, imc=imc, imb=imb0

; FUNCTION/PROCEDURE to CalCulate diameter of the feature




imb=byte(imb0)

in=size(imb)
nx=in[1]
ny=in[2]

; kernel for 8-neighbourhood
		;kernel=replicate(1, 3, 3)
		;inval=9
; kernel for 4-neighbourhood
kernel=[[0, 1, 0], [1, 1, 1], [0, 1, 0]]
inval=5

;help, kernel
;print, max(imb)

if not arg_present(imc) then begin
	;			chaincode, imb, cn, fp, imc=imc
				;im0=bytarr(nx+2, ny+2)
				im0= replicate(1, nx+2, ny+2)
				im0[1:nx, 1:ny]=imb
				im0=convol(im0, kernel)
				imc0=im0[1:nx, 1:ny]
				locs=where(imc0 eq inval)
				if locs[0] ne -1 then begin
						imc1=imb
						imc1[locs]=0

		;				if max(imc-imc1) ne 0 then stop
						imc=imc1
					endif else imc=imb

	endif
;imb=contourfill(imc)

in=size(imc)
nx=in[1]
ny=in[2]
clocs=where(imc ne 0)
n=n_elements(clocs)
xp= clocs mod nx
yp= clocs / nx


;				D0=dblarr(n, n)
;
;				for i=0, n-1 do $
;					for j=0, n-1 do $
;
;					D0[i, j]=sqrt(double((xp[i]-xp[j])*(xp[i]-xp[j]))+ $
;									double((yp[i]-yp[j])*(yp[i]-yp[j])))
;
;				diameter= max(D0)
;				locs=where(D0 eq diameter)
;				;help, locs
;				i0= locs mod n
;				j0= locs / n
;
				;print, 'First Pixel:', xp[i0[0]], yp[i0[0]]
;				;print, 'Second Pixel:', xp[j0[0]], yp[j0[0]]
;				print, 'Diameter:', diameter

xmax=max(xp)
xmin=min(xp)

ymax=max(yp)
ymin=min(yp)

Dmax=double(xmax-xmin)^2+double(ymax-ymin)^2
if xmax-xmin gt ymax-ymin then begin

	ylen=double(ymax-ymin)^2
endif else begin
	xp0=xp
	xp=yp
	yp=xp0

	ylen=double(ymax-ymin)^2
endelse

	xi=sort(xp)

	Dest= double(xp[xi[0]]- xp[xi[n-1]])^2+double(yp[xi[0]]- yp[xi[n-1]])^2


;	print, 'Dest 1:', sqrt(Dest)
	s0=0
	s1=n-1

	count=0
	while (xp[xi[s1]]-xp[xi[s0]])^2 gt Dest - ylen do begin


		while (xp[xi[s1]]-xp[xi[s0]])^2 gt Dest - ylen do begin

				D0=double(xp[xi[s0]]- xp[xi[s1]])^2+double(yp[xi[s0]]- yp[xi[s1]])^2
				if D0 gt Dest then Dest=D0
				s1=s1-1
				count=count+1
				if count gt n*n then stop
		endwhile
		s0=s0+1
		s1=n-1
	endwhile
;	print, 'Dest 2:', sqrt(Dest)


;	if abs(sqrt(Dest)-Diameter) gt 1.d-06 then stop

;print, n*n
;print, count

;return, diameter

return, sqrt(Dest)

end