;============================================================================

pro frclean_efiterr, binary, xc, yc, vector, binary2, stdev, stdevGeo, abs_err, geoErr, display=display

;FUNCTION
;	calculating the absolute error of the fit
;
;	INPUT
;		binary		binary image with the data
;		xc, yc		coordinate centre
;		vector		ellipse parameters
;
;	OUTPUT
;		binary2		the cleaned binary with the outliers removed
;		stdev		fits standard deviation (algebraic)
;		stdevGeo	fits standard deviation (geometric)
;		algerr		fits algebraic error




locs=where(binary)
N=n_elements(locs)
d=dblarr(6,N)
err=dblarr(N)
geo_err=dblarr(N)
info=size(binary)
nx=info[1]

abs_err=0
geoErr=0
for i=0l,N-1,1 do begin

; calculate alebraic error

	yi1=fix(locs[i]/nx)
	xi1=locs[i]-long(yi1)*nx
	yi=double(yi1)-yc
	xi=double(xi1)-xc

	d[0,i]=xi*xi
	d[1,i]=xi*yi
	d[2,i]=yi*yi
	d[3,i]=xi
	d[4,i]=yi
	d[5,i]=1.d00

	err[i]=0
	for k=0, 5, 1 do err[i]=err[i]+vector[k]*d[k, i]
	abs_err=abs_err+abs(err[i])
	;	geometric error

	if xi ne 0 then begin
		cof=yi/xi
		A1=vector[0]+vector[1]*cof+vector[2]*cof*cof
		B1=vector[3]+vector[4]*cof
		C1=vector[5]
		x1=(-B1+sqrt(B1*B1-4*C1*A1))/(2*A1)

		if x1*xi gt 0 then x=x1 else x=(-B1-sqrt(B1*B1-4*C1*A1))/(2*A1)
		y=cof*x
		geo_err[i]=sqrt((xi-x)*(xi-x)+(yi-y)*(yi-y))

	endif else begin
		A1=vector[2]
		B1=vector[4]
		C1=vector[5]
		y1=(-B1+sqrt(B1*B1-4*C1*A1))/(2*A1)
		if y1*yi ge 0 then y=y1 else y=(-B1-sqrt(B1*B1-4*C1*A1))/(2*A1)
		geo_err[i]=abs(yi-y)
	endelse
geoErr=geoErr+geo_err[i]*geo_err[i]

if err[i] lt 0 then geo_err[i]=-geo_err[i]

endfor

geoErr=sqrt(geoErr/N)
stdevGeo=stddev(geo_err, /double)
resultGeo=moment(geo_err, /double, mdev=dev)


result=moment(err, /double, mdev=dev)
stdev=stddev(err, /double)

;res=dialog_message('Standard Deviation (%):'+string(100.*stdev)+'   Mean: ' + string( result[0] ) + '  Variance: ' + string( result[1] ) + $
;   '  Skewness: ' + string ( result[2] ) + '  Kurtosis: ' + string (result[3]), /information)


if keyword_set(display) then begin
	print, 'Geometric Error:', geoErr, '  pixels'
	print, 'Standard Deviation (pixels):', stdevGeo & $
	PRINT, 'Mean: ', resultGeo[0] & PRINT, 'Variance: ', resultGeo[1] ;& $
;	  PRINT, 'Skewness: ', resultGeo[2] & PRINT, 'Kurtosis: ', resultGeo[3]
	print

	print, 'Algebraic Absolute Error:',abs_err
	print, 'Standard Deviation (%):', 100*stdev & $
	PRINT, 'Mean: ', result[0] & PRINT, 'Variance: ', result[1] ;& $
;	  PRINT, 'Skewness: ', result[2] & PRINT, 'Kurtosis: ', result[3]
	print
	print
endif


;cleaning the input binary image by removing pixels with bigger errors

thr_locs=where(abs(geo_err) ge max(geo_err)/1.1)
;thr_locs=where(abs(err) ge max(err)/1.05)
binary2=binary
binary2[locs[thr_locs]]=0

end

