function GaussianSmoothing, image, n

;Generate Kernel

if n gt 7 then n=7
if n lt 2 then n=2
w=2*n+1
a=5.0/w*n
smal=exp(-a*a/2.0)

;mask2=intarr(n, n)
mask=dblarr(15)
mask[0]=1
sum=1

image2=double(image)
for i = 1, w-1 do begin
		a = 5.0/w*(n-i) ;
		mask[i] = exp(-a*a/2.0)/smal ;+ 0.5 ;
		sum = sum + mask[i] ;
	endfor

in=size(image)
nx=long(in(1)) & ny=in(2)
image2(*)=0
for j=n, nx-1-n do $
   for i=n, ny-1-n do begin
	value=0
	ijn=i*nx+j-n
	for k=0, w-1 do value=value+mask[k]*image[ijn+k]
	image2[j, i]=value/sum	
endfor
image3=image2
image3(*)=0

for j=n, nx-1-n do $
   for i=n, ny-1-n do begin
	value=0
	ijn=(i-n)*nx+j
	for k=0, w-1 do value=value+mask[k]*image2[ijn+k*nx]
	image3[j, i]=value/sum	
endfor
return, image3
end
