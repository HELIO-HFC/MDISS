pro PrMeudSyn, ew, ns, header0, header, image0, image


ny=fxpar(header0, 'NAXIS2')
date=fxpar(header0, 'DATE_OBS')
dateN=fxpar(header, 'DATE_OBS')

Rold=fxpar(header0, 'SOLAR_R')
Rnew=fxpar(header, 'SOLAR_R')
Xold=fxpar(header0, 'CENTER_X')
Yold=fxpar(header0, 'CENTER_Y')

Xnew=fxpar(header, 'CENTER_X')
Ynew=fxpar(header, 'CENTER_Y')

Cd1o=fxpar(header0, 'CDELT1')
Cd2o=fxpar(header0, 'CDELT2')

Cd1n=fxpar(header, 'CDELT1')
Cd2n=fxpar(header, 'CDELT2')

if date ne dateN then return
print, 'OLD    Radius:', Rold, '    Center:', Xold, Yold
print, 'NEW    Radius:', Rnew, '    Center:', Xnew, Ynew

crot=tim2clon(date)
clon=tim2carr(date)
print, clon, crot, date


arcsec=60*hel2arcmin(ns, ew-clon, date=date)
print, arcsec
print, Cd1o, cd2o

xpo=fix(arcsec[0]/cd1o+xold+.5)
ypo=fix(arcsec[1]/cd2o+ny-yold+.5)
print, xpo, ypo

xpn=fix(arcsec[0]/cd1n+xnew+.5)
ypn=fix(arcsec[1]/cd2n+ynew+.5)

imx=image0[xpo-50:xpo+50, ypo-50:ypo+50]
imxn=image[xpn-50:xpn+50, ypn-50:ypn+50]
window, 10, xs=100, ys=100, title='Original'
tvscl, imx

window, 12, xs=100, ys=100, title='Processed'
tvscl, imxn

window, 11, xs=1024, ys=1024
tvscl, image0
xyouts, xpo, ypo, 'X', alignment=.5,/device, color=0

window, 13, xs=1024, ys=1024
tvscl, image
xyouts, xpn, ypn, 'X', alignment=.5, /device, color=0

end