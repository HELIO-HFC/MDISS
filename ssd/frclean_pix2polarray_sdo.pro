
function frclean_pix2polarray_sdo, image, R, xc, yc

;	FUNCTION
;			transferring rectangular image to polar coordinate
;
;		INPUT
;				image			image
;				R				disk Radius
;				xc, yc			disk centre coordinates
;
;		RESULT
;				rectangular array corresponding to the initial disk of radius R
;				centred at xc, yc

info=size(image)

image1x = info[1]
image1y = info[2]

Rmax=fix(R)
Rmax2 = Rmax*Rmax
R2PI = fix (2 * !DPI * Rmax +0.5)

pixel_array = MAKE_ARRAY(image1x, Rmax, /float, VALUE = 0)

FOR i=0, Rmax-1 DO BEGIN
	FOR j=0, image1x-1 DO BEGIN
		; i  r
		; j is the widith of Array
		;	theta = j * ((2*!DPI)/image1x)
		theta = j * ((2*!DPI)/(image1x-1))
		x = i * COS(theta)
		y = i * SIN(theta)
		X= fix(x+xc+0.5)>(0l)<(image1x-1l)
		Y= fix(y+yc+0.5)>(0l)<(image1y-1l)

		pixel_array(j,i) = image(X,Y)
	ENDFOR
ENDFOR

return, pixel_array

end

;==========================================================================