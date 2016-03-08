pro MDIsync, mdi, halpha, Result, dir=dir0, soho=soho, earth=earth

	; Read In MDI image

fn=dialog_pickfile(path='Y:\mdi_data\mdi_whitel_april_02')
;fn=dialog_pickfile(path='Y:\mdi_data\mg_april_02')
;fn=dialog_pickfile(path='Y:\meudon\processed\')
image=readfits(fn, header)
bim=bytscl(image, max=1000, min=-1000)
GetHeaderInfo, header, nx, ny, date, xc, yc, R, cdx, cdy

	; Create mdi structure
mdi={image:image, display:bim, R:0.d00, nx:nx, ny:ny, $
		xc:xc, yc:yc, header:header, date:date, $
		cdx:cdx, cdy:cdy}



if R eq 0 then R=420.d00
mdi.R=R

;mdi.display=mdi.image

	; Read in H-alpha image
fn2=dialog_pickfile(path='Y:\mdi_data\mdi_whitel_april_02')
;fn2=dialog_pickfile(path='Y:\mdi_data\mg_april_02')
;fn2=dialog_pickfile(path='Y:\meudon\processed\');Halpha\20021114\subtract')
image=readfits(fn2, header)
GetHeaderInfo, header, nx, ny, date, xc, yc, R, cdx, cdy


halpha={image:image, R:0.d00, nx:nx, ny:ny, $
		xc:xc, yc:yc, header:header, date:date, $
		cdx:cdx, cdy:cdy}
if R eq 0 then R=420.d00
halpha.R=R

;-------------------------------------------------

	;	Create filename for saving the processed image
len=strlen(mdi.date)
nl=strpos(mdi.date, 'T')
timex=strmid(mdi.date, nl, 9)
print, timex
pos= strsplit(timex, ':', length=len)
k=n_elements(pos)
time=strmid(timex, pos[0], pos[1]-pos[0]-1)

for i=1, k-1 do begin

		time=time+strmid(timex, pos[i], len[i])

endfor
print, time


path1='Y:\mdi_data\mg_april_02\fd_M_96m_'
len=strlen(fn)
nl=strpos(fn, '000')


fname=strmid(fn, nl, len-nl)
foldername=strmid(fn, nl-9, 9)

print, fname


len=strlen(fn2)
nl=strpos(fn2, '.fits')

fname2=strmid(fn2, nl-34, 15)
print, fname2

print, path1+foldername+fname2+'.'+time+'.'+fname
writepath=path1+foldername+fname2+'.'+time+'.'+fname



;-----------------------------------------------------


Result=intarr(mdi.nx, mdi.ny)

print, 'Magnetogram:', mdi.date
print, 'Halpha:', halpha.date

dfr=abs(anytim(mdi.date)-anytim(halpha.date))
print, 'The Difference is ', abs(anytim(mdi.date)-anytim(halpha.date)), '  seconds'

;if dfr gt 4000 then stop
if not keyword_set(earth) then mangle=pb0r(mdi.date, /soho) else mangle=pb0r(mdi.date, /earth)
mrx=1.0/tan(mangle(2)/60/!radeg)
mclon=tim2carr(mdi.date)
mb0 = mangle[1]/!radeg
;if keyword_set(earth) then begin
		mdi.cdx=mangle[2]*60/mdi.R
		mdi.cdy=mangle[2]*60/mdi.R
;endif
;use_earth_view
if keyword_set(soho) then hangle=pb0r(halpha.date, /soho) else hangle=pb0r(halpha.date, /earth)
hrx=1.0/tan(hangle(2)/60/!radeg)
hclon=tim2carr(halpha.date)
hb0 = hangle[1]/!radeg

;if not keyword_set(soho) then begin
		halpha.cdx=hangle[2]*60/halpha.R
		halpha.cdy=hangle[2]*60/halpha.R
;endif
imRoi=bytarr(mdi.nx, mdi.ny)
for i=0, mdi.nx-1 do $
	for j=0, mdi.ny-1 do $
		if (mdi.xc-i)*(mdi.xc-i)+(mdi.yc-j)*(mdi.yc-j) lt mdi.R*mdi.R then $
			imRoi[i, j]=1
locs=where(imRoi ne 0)
;locs=where(mdi.image gt -10000)
n=n_elements(locs)
xp=locs mod mdi.nx
yp=locs / mdi.nx
use_earth_view

arcx=(double(xp-mdi.xc)*mdi.cdx)/60.d00
arcy=(double(yp-mdi.yc)*mdi.cdy)/60.d00
if keyword_set(earth) then helio=arcmin2hel(arcx, arcy, date=mdi.date) else helio=arcmin2hel(arcx, arcy, date=mdi.date, /soho)
helio[1, *]=helio[1, *]+mclon[0]
helio[0, *]=helio[0, *]

if keyword_set(soho) then $
	harc=hel2arcmin(helio[0,*], helio[1,*]-hclon[0], date=halpha.date, /soho) $
		else harc=hel2arcmin(helio[0,*], helio[1,*]-hclon[0], date=halpha.date)
hxp=harc[0, *]*60.d00/halpha.cdx+halpha.xc
hyp=harc[1, *]*60.d00/halpha.cdy+halpha.yc

for i=0l, n-1 do begin
	Result[locs[i]]=halpha.image[fix(hxp[i]+.5), fix(hyp[i]+.5)]
endfor


window, 0, xs=mdi.nx, ys=mdi.ny
tvscl, result

window, 1, xs=mdi.nx, ys=mdi.ny
tvscl, halpha.image

window, 2, xs=mdi.nx, ys=mdi.ny
tvscl, mdi.image
CheckResize, mdi.image, result

;		im=result
;		im2=mdi.display
;		locs=where(mdi.image lt 5000)
;		;locs=where(abs(mdi.image) gt 700 and abs(mdi.image) lt 3000)
;		im[locs]=0
;		im2[locs]=0
;		window, 3, xs=mdi.nx, ys=mdi.ny
;		tvscl,im
;
;		window, 4, xs=mdi.nx, ys=mdi.ny
;		tvscl,im2

;writefits, path1+foldername+fname2+'.'+time+'.'+fname, Result, halpha.header
writefits, path1+'.'+time+'.'+fname, Result, halpha.header
;window, 5, xs=mdi.nx, ys=mdi.ny
;tvscl,mdi.image-result
end


pro CheckResize, mdim,  result


window, 5, xs=1024, ys=1024
locs=where(mdim lt 6000. and mdim gt 0)
;locs2=where(mdim lt 4000. and mdim gt 0)
im=result
im[locs]=max(result)
tvscl, im

window, 6, xs=1024, ys=1024
tvscl, mdim-result
end