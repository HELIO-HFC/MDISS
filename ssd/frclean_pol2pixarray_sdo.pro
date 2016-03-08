
function frclean_pol2pixarray_sdo, image, pixel_array, R1, xc, yc


;	FUNCTION
;			transferring polar coordinate image to rectangular
;
;		INPUT
;				image			?????
;				pixel_array		image
;				R1				disk Radius
;				xc, yc			disk centre coordinates
;
;		RESULT
;				rectangular array corresponding to the initial disk of radius R
;				centred at xc, yc



info=size(image)

image1x = info[1]
image1y = info[2]


imbg = MAKE_ARRAY(image1x,image1y, /float, VALUE = 0)
Rmax=fix(R1)
FOR i=0, image1x-1 DO BEGIN
	FOR j=0, image1y-1 DO BEGIN

		ixc = (i - float(xc))
		jyc = (j - float(yc))
		rr = jyc^2 + ixc^2
		r  = long(SQRT(rr)+0.5)

		IF r LT Rmax THEN BEGIN

			theta = ATAN(jyc,ixc)

			IF (theta LT 0)THEN theta = ((2*!DPI) + theta)
			;theta1 = LONG((theta *1023)/(2*!DPI) + 0.5)
			;	theta1 = LONG((theta *image1x)/(2*!DPI) + 0.5)
			theta1 = fix((theta *(image1x-1))/(2*!DPI) + 0.5)
			imbg(i,j) = pixel_array(theta1,r)



       ENDIF ;r LE Rmax-1 (i.e. inside the disk)
	ENDFOR
ENDFOR

return, imbg
end

