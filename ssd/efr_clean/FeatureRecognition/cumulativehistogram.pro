function CumulativeHistogram, image ;, nc

h=histogram(image)
n=n_elements(h)
hC=h
hCdif=h
hCdif2=h
hCdif3=h
sum=0

;	 build inverted cumulative histogram

for i=0l, n-1, 1 do begin
	sum=sum+h[n-1-i]
	hC[n-1-i]=sum
	endfor

return, hc
end

;=======================================================

function getpixelcoord, locs, nx

n=n_elements(locs)
IndSpot=intarr(2, n)
for i=0l, n-1, 1l do begin
	IndSpot[0, i] = locs[i] mod nx
	IndSpot[1, i] = locs[i] / nx
endfor

return, IndSpot
end

;=======================================================

function CumHist, hist

sum=0
n=n_elements(hist)
hc=intarr(n)
for i=0, n-1, 1 do begin
	sum=sum+hist[n-1-i]
	hC[n-1-i]=sum
	endfor

return, hc
end