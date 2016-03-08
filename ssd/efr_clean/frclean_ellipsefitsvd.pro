;=============================================================================



pro frclean_ellipsefitsvd,  binary, xc, yc, xc1, yc1, R1, R2, theta, stdev, stdevGeo, algerr

;	PROCEDURE
;			fitting the ellipse to the binary image
;
;	INPUT
;		binary		cleaned binary image
;		xc, yc		gravity centre coordinates in pixel coordinates
;
;	OUTPUT
;		xc1, yc1	ellipse centre coordiinates
;		R1, R2		ellipse axis
;		theta		ellipse angle
;		stdev		fits standard deviation (algebraic)
;		stdevGeo	fits standard deviation (geometric)
;		algerr		fits algebraic error


;fit ellipse to the original image
binary_copy=binary
vector=FitEllipseSVD(binary,  xc, yc)

;and obtain ellipse parameters
frclean_ealgdata, xc, yc, vector, xc1, yc1, R1, R2, theta, /noprint

; change coordinate centre to the centre of the first fit
xc=xc1
yc=yc1
count=0
stdev=1.

; while stdev is large remove the furthest oullier and fit the ellipse
while stdev ge .3e-02 and count lt 200 do begin
		count=count+1

		vector=FitEllipseSVD(binary,  xc, yc)

		frclean_efiterr, binary, xc, yc, vector, binary2, stdev

		binary=binary2
endwhile

; perform the final fit and obtain ellipse parameters
		vector=FitEllipseSVD(binary,  xc, yc)
		frclean_ealgdata, xc, yc, vector, xc1, yc1, R1, R2, theta


; calculate the final fit errors and print it onscreen
frclean_efiterr, binary, xc, yc, vector, binary2, stdev, stdevGeo, algerr, geoErr ;, /display

end
