;----------------------------------------------------------------------------
; FILE::GETFileName, Path, Extension
;
; Purpose:
;  Creates a map object by rereading fitsfile at location
;
FUNCTION FitsFDObs::GetMap, test=test, soho=soho

	fn=self.location
if not keyword_set(test) then begin
	mreadfits, fn, index, data
	index2map, index, data, map
 endif else begin
	xc=((self.nx-1)/2. - self.xc) * self.cdx
	yc=((self.ny-1)/2. - self.yc) * self.cdy
	map=make_map(*self.image, soho=soho, xc=xc, yc=yc, $
					dx=self.cdx, dy=self.cdy, time=self.date)
 endelse
;stop
return, map

END