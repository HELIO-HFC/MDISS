pro GaussianSmoothing, image, n

;Generate Kernel

if n gt 7 then n=7
if n lt 2 then n=2
w=2*n+1
a=5.0/w*n
smal=exp(-a*a/2.0)

mask2=intarr(n, n)
mask=intarr(15)
mask[0]=1
sum=1


for i = 1, w-1 do begin
		a = 5.0/w*(n-i) ;
		mask[i] = exp(-a*a/2.0)/smal + 0.5 ;
		sum = sum + mask[i] ;
	endfor

print, mask


for i=0, n-1 do $
	for j=0, n-1 do begin

endfor
end