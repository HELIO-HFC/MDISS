function GetSpotList, wlfn=wlfn, count=n0

fn=findfile('Y:\mdi_data\mdi_whitel_april_02\*', count=k)

;print, k

list=strarr(k*4)
wlfn=strarr(k*4)
n0=0

for i0=2, k-1 do begin

	fn2=findfile(fn[i0]+'Processed*.dat', count=l)
;	help, fn2
	list[n0:n0+l-1]=fn2
;	WLfn=strarr(l)
for i1=0, l-1 do begin
	ln= strlen(fn2[i1])
	WLfn[n0+i1]= strmid(fn2[i1], 0, ln-9)+'.fits'
end

	n0=n0+l
endfor

list=temporary(list[0:n0-1])
wlfn=temporary(wlfn[0:n0-1])
return, list
end
