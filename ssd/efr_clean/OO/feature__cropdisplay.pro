;=========================================================

PRO Feature::CropDisplay, FD, scale=scale, exact=exact, x0=x0, y0=y0

if not keyword_set(scale) then scale=1
self->GetBndRct, x0, y0, x1, y1
if keyword_set(exact) then $1
	FD->CropDisplay, x0, y0,  x1, y1, scale=scale $
	else FD->CropDisplay, x0-10, y0-10,  x1+10, y1+10, scale=scale

x0=x0-10
y0=y0-10
END
