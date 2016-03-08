pro CorrelateSpotData, data


data1=ChooseArcx(data, 500)
data2=ChooseNumbras(data, 1, 2)
data3=ChooseUmbraSize(data, 1, 5)

pixSizeArray=data2[*, 5]
helSizeArray=data2[*, 7]

umbPixSizeArray=data2[*, 6]
helPixSizeArray=double(data2[*, 7])*double(data2[*, 6])/data2[*, 5]

totFluxArray=abs(data2[*, 8])/(4*!pi)
maxFluxArray=abs(data2[*, 10])/(4*!pi)

irradArray=(data[*, 12]-data[*, 13])/data[*, 13]

print, correlate(data2[*, 7], abs(data2[*, 8])/(4*!pi))

end

;======================================================================

function ChooseArcx, data, aSecLimX

; FUNCTION:
;
;		returns the array of sunspots with absolute value of
;		barycenter x-coordinate less than aSecLimX
info=size(data)
n1=info[1]
nz=info[2]
n2=0
data1=dblarr(n1, nz)

for i=0, n1-1 do $
	if abs(data[i, 0]) lt aSecLimX and abs(data[i, 0]) gt .1d-5 $
		and abs(data[i, 1]) lt 600 then begin
		data1[n2, *]=data[i, *]
		n2=n2+1
	endif

data2=data1[0:n2-1, *]

return, data2
end

;======================================================================

function ChooseNumbras, data, Nmin, Nmax

; FUNCTION:
;
;		returns the array of sunspots with the number of umbras
;		greater than N

if Nmax lt Nmin then begin
	print, 'N max should be greater or equal to N min'
	return, data
end
info=size(data)
n1=info[1]
nz=info[2]
n2=0
data1=dblarr(n1, nz)

for i=0, n1-1 do $
	if abs(data[i, 4]) le Nmax $
		and abs(data[i, 4]) ge Nmin $
		and abs(data[i, 0]) gt .1d-5 then begin
		data1[n2, *]=data[i, *]
		n2=n2+1
	endif

data2=data1[0:n2-1, *]

return, data2
end

;======================================================================

function ChooseUmbraSize, data, Nmin, Nmax

; FUNCTION:
;
;		returns the array of sunspots with umbras of size
;		greater or equal to Nmin and less or equal to Nmax

info=size(data)
n1=info[1]
nz=info[2]
n2=0
data1=dblarr(n1, nz)

for i=0, n1-1 do $
	if data[i, 4] gt 0 $
		and abs(data[i, 6]) le Nmax $
		and abs(data[i, 6]) ge Nmin $
		and abs(data[i, 0]) gt .1d-5 then begin
		data1[n2, *]=data[i, *]
		n2=n2+1
	endif

data2=data1[0:n2-1, *]

return, data2
end

;======================================================================

function ChooseHemisphere, data, north=north, south=south

; FUNCTION:
;
;		returns the array of sunspots with positive or negative
;		latitude

info=size(data)
n1=info[1]
n2=0
data1=dblarr(n1, 14)

if keyword_set(north) then $
	for i=0, n1-1 do $
			if data[i, 3] ge 0 $
			and abs(data[i, 0]) gt .1d-5 then begin
			data1[n2, *]=data[i, *]
			n2=n2+1
		endif

if keyword_set(south) then $
	for i=0, n1-1 do $
		if	data[i, 3] lt 0 $
			and abs(data[i, 0]) gt .1d-5 then begin
			data1[n2, *]=data[i, *]
			n2=n2+1
		endif

data2=data1[0:n2-1, *]

return, data2
end

;======================================================================

pro DataPlot, x, y, nofit=nofit

window, 0, xs=1000, ys=700

; Plot the axis
plot, x, y, title='Full Disk: area versus absolute value of magnetic flux', $
		ytitle='Flux, Gauss', xtitle='Size, Degrees', $
		/nodata

xyouts, x, y, 'x', alignment=.5

print, correlate(x, y)
if not keyword_set(nofit) then begin
		lf=linfit(x, y)
		print, lf

		ind=indgen(n_elements(x))
		oplot, ind, lf[0]+lf[1]*ind, thick=2
	endif
end


pro fitting, x, y

degree=5
cf= poly_fit(x, y, degree)

z=poly(x, cf)
ind=sort(x)
window, 0, xs=1024, ys=700
plot, x[ind], z[ind]

xyouts, x, y, 'x', alignment=.5
print, cf
end