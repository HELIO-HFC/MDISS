pro MDI_spotsmagnetic

fn=dialog_pickfile(path='Y:\mdi_data\mdi_whitel_april_02')
image=readfits(fn, header)

getFileName,fn, filename, extension
foldername=getFolderName(fn)
GetHeaderInfo, header, nx, ny, date, xc, yc, R, cdx, cdy

im_or={image_str,  image:ptr_new(/allocate_heap), header:ptr_new(/allocate_heap), location:'' }
IniImage_str, im_or, image, header, fn

temp={fits_str, nx:0, ny:0, date:'', xc:0.d00, yc:0.d00, R:0.d00, cdx:0.d00, cdy:0.d00, foldername:'', filename:'', extension:''};, inherits image_str };}
wl_or={file_str, image: ptr_new(/allocate_heap), fitsinfo: ptr_new(/allocate_heap)}
;wl_or={detected_str, image:image_str, fits: fits_str, detected:ptr_new(/allocate_heap)}

wl_or.image=ptr_new(im_or)

im2=obj_new('image_str')
im2.image=ptr_new(image)
im2.location=fn
stop
IniFile_str, wl_or, image, header, nx, ny, date, xc, yc, R, cdx, cdy, fn, filename, extension, foldername

window, 0, xs=wl_or.nx, ys=wl_or.ny
tvscl, *wl_or.image

print, wl_or.filename
print, wl_or.foldername
print, wl_or.extension

fn2=wl_or.foldername+wl_or.filename+'_detected'+wl_or.extension

image=readfits(fn2, header)
GetHeaderInfo, header, nx, ny, date, xc, yc, R, cdx, cdy


window, 1, xs=wl_or.nx, ys=wl_or.ny
tvscl, image
;stop
end

;================================================

pro IniImage_str, im, image, header, fn

im.image=ptr_new(image)
im.header=ptr_new(header)
im.location=fn

end

;=======================================================

pro IniFile_str, struct, image, header, nx, ny, date, xc, yc, R, cdx, cdy, fn, filename,extension,  foldername

struct.image=ptr_new(image)
struct.header=ptr_new(header)
struct.location=fn
struct.nx=nx
struct.ny=ny
struct.date=date
struct.xc=xc
struct.yc=yc
struct.R=R
struct.cdx=cdx
struct.cdy=cdy
struct.filename=filename
struct.foldername=foldername
struct.extension=extension
end
;=====================================================

function getFolderName, fn

k0=5
i1=strpos(fn, '.fits', /REVERSE_SEARCH)
if i1 eq -1 then begin
	i1=strpos(fn, '.fts', /REVERSE_SEARCH)
	if i1 eq -1 then return, i1 else k0=4
endif


i2=strpos(fn, '\', /REVERSE_SEARCH)
if i2 eq -1 then i2=strpos(fn, '/', /REVERSE_SEARCH)
if i2 eq -1 then return, i2

fname= strmid(fn, 0, i2+1);, i1-i2+k0-1)
return, fname
end

;================================================

pro getFileName, fn, fname, ext

k0=5
ext='.fits'
i1=strpos(fn, ext, /REVERSE_SEARCH)
if i1 eq -1 then begin
	ext='.fts'
	i1=strpos(fn, ext, /REVERSE_SEARCH)
	if i1 eq -1 then fname=-1 else k0=4
							;endelse
endif


i2=strpos(fn, '\', /REVERSE_SEARCH)
if i2 eq -1 then i2=strpos(fn, '/', /REVERSE_SEARCH)
if i2 eq -1 then fname=-1

fname= strmid(fn, i2+1, i1-i2-1)

end