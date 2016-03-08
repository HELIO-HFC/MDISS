pro ChainCode, imb, chain, firstp, n8=n8, imc=imc

; fill the gaps in the regions (experimental(
	imz= ((watershed(imb) gt 1)+imb) ge 1
	imb0=imb
	imb=imz
;
in=size(imb)
nx=in[1]
ny=in[2]
locs=where(imb ne 0)
n=n_elements(locs)
xp= locs mod nx
yp= locs / nx

imc=bytarr(nx, ny)

if keyword_set(n8) then $
	for i=0l, n-1 do begin
		value=imb[xp[i]-1, yp[i]-1]*imb[xp[i], yp[i]-1] $
				* imb[xp[i]+1, yp[i]-1] * imb[xp[i]-1, yp[i]] $
				* imb[xp[i]+1, yp[i]] * imb[xp[i]-1, yp[i]+1] $
				* imb[xp[i], yp[i]+1] * imb[xp[i]+1, yp[i]+1]

		if value eq 0 then imc[xp[i], yp[i]]=1
	endfor else $
	for i=0l, n-1 do begin
		value=imb[xp[i]-1, yp[i]]*imb[xp[i]+1, yp[i]] $
			 * imb[xp[i], yp[i]-1]*imb[xp[i], yp[i]+1]
		if value eq 0 then imc[xp[i], yp[i]]=1
	endfor

;window, 1, xs=nx, ys=ny
;tvscl, imc

clocs=where(imc ne 0)
n=n_elements(clocs)
xp= clocs mod nx
yp= clocs / nx

firstp=intarr(2)
Chain=intarr(10*n)


i0=where(xp eq min(xp))
firstp[0]=xp[i0[0]]
firstp[1]=yp[i0[0]]

cp=intarr(2)
cp=firstp
source=2
count=0

;
;		1	2	3
;		0	X	4
;		7	6	5

movearray=intarr(8, 2)
movearray[0, 0]=-1
movearray[0, 1]=0

movearray[1, 0]=-1
movearray[1, 1]=1

movearray[2, 0]=0
movearray[2, 1]=1

movearray[3, 0]=1
movearray[3, 1]=1

movearray[4, 0]=1
movearray[4, 1]=0

movearray[5, 0]=1
movearray[5, 1]=-1

movearray[6, 0]=0
movearray[6, 1]=-1

movearray[7, 0]=-1
movearray[7, 1]=-1
while (cp[0] ne firstp[0] or cp[1] ne firstp[1] and count lt 10*n) or count eq 0 do begin

		direction=(source+6) mod 8
		tmp=cp+movearray[direction, *]
		count0=0
		while imc[tmp[0], tmp[1]] eq 0 do begin
			direction=(direction+1) mod 8
			tmp=cp+movearray[direction, *]
			count0=count0+1
			if count0 gt 8 then if n ne 1 then begin
				stop
				tvim, im

					endif else begin
						chain[count]=-1
						count=1
						goto, L0
				endelse
		end
		cp=tmp
		chain[count]=direction
		source=direction
		count=count+1

endwhile

L0:	chain0=chain[0:count-1]
	chain=chain0
	imb=imb0

;	tvim, chain2p(chain, firstp, nx, ny)-imc
;stop
end

;===================================================

function chain2p, ch, fp, nx, ny

im=bytarr(nx, ny)
movearray=intarr(8, 2)
movearray[0, 0]=-1
movearray[0, 1]=0

movearray[1, 0]=-1
movearray[1, 1]=1

movearray[2, 0]=0
movearray[2, 1]=1

movearray[3, 0]=1
movearray[3, 1]=1

movearray[4, 0]=1
movearray[4, 1]=0

movearray[5, 0]=1
movearray[5, 1]=-1

movearray[6, 0]=0
movearray[6, 1]=-1

movearray[7, 0]=-1
movearray[7, 1]=-1

n=n_elements(ch)
im[fp[0], fp[1]]=1
cp=fp

print, n
if n gt 1 then $
	for i=0, n-1 do begin

		cp=cp+movearray[ch[i], *]
	;	if im[cp[0], cp[1]] ne 0 then stop
		im[cp[0], cp[1]]=1
	endfor

return, im
end

;========================================================

function ContourFill, imc

in=size(imc)
nx=in[1]
ny=in[2]
im=bytarr(nx, ny)
im[*]=1
res=search2d(imc, 0, 0, 0, 0)
im[res]=0

;window, 0, xs=nx, ys=ny
;tvscl, im

return, im
end
