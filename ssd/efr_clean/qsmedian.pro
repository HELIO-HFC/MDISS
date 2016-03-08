;===================================================================================================

function qsmedian, image, R, xc, yc

;	PROCEDURE
;			Constructing Quiet Sun from the full disk image using median
;
;		INPUT
;				image 			input image
;				R				Disk Radius
;				xc, yc			disk centre coordinate
;
;
;		RETURN
;				MedianImage		image cleaned via Median
;



p_arr_median2=frclean_pix2polarray_sdo(image, R, xc,yc)

info=size(image)
nx=info[1]
ny=info[2]
imin=size(p_arr_median2)


;recontructed_original=recontruct_from_pixel( image, p_arr_median2, R, xc,yc)

minx=imin[1]
miny=imin[2]


for i=0, miny-1 do begin
	value=median(p_arr_median2[*, i], /even)
	p_arr_median2[*,i]=value
endfor

reconstructed_median2=frclean_pol2pixarray_sdo(image, p_arr_median2, R, xc,yc)
;stop
return, reconstructed_median2

end

