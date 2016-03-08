pro GetHeaderInfo, header, nx, ny, date, xc, yc, R, cdx, cdy, map


nx=fxpar(header, 'NAXIS1')
ny=fxpar(header, 'NAXIS2')
t_start=fxpar(header, 'T_START')
t_stop=fxpar(header, 'T_STOP')
date=fxpar(header, 'DATE_OBS')

xc=fxpar(header, 'CENTER_X')
yc=fxpar(header, 'CENTER_Y')

xc=fxpar(header, 'X0')
yc=fxpar(header, 'Y0')


R=fxpar(header, 'R_SUN')
if R eq 0 then R=fxpar(header, 'SOLAR_R')
cdx=fxpar(header, 'CDELT1')
cdy=fxpar(header, 'CDELT2')

;if arg_present(map) then map={image:intarr(nx, ny), detected:intarr(nx, ny), date:date, xc:xc, yc:yc, cdx:cdx, cdy:cdy, R:R, nx:nx, ny:ny}
end