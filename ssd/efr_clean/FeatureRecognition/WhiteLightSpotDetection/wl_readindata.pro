;===================================================================

pro WL_ReadInData, image, nx, ny, xc, yc, Radius, qsint, imdi, bimdi, wls, mdate, fn=fn

;	PROCEDURE
;		reading in MDI white light image (chosen manually) and corresponding MDI
;		magnetogram (automatically) as well as determining the quiet sun intensity level
;		based on white light image histogram
;
;
;	OUTPUT
;		image			white light image
;		nx, ny			image size (as in the header)
;		xc, yc, Radius	disk center coordinates and radius
;		qsint			quiet sun intensity
;		imdi			magnetogram array
;		bimdi			magnetogram array rescaled for viewing
;
;	NOTES
;			magnetogram parameters (nx, ny, xc, yc, radius) are assumed the same as white light
if not keyword_set(fn) then fn=dialog_pickfile(path='Y:\mdi_data\mdi_whitel_april_02\')

image=readfits(fn, header)


nx=fxpar(header, 'NAXIS1')
ny=fxpar(header, 'NAXIS2')
t_start=fxpar(header, 'T_START')
t_stop=fxpar(header, 'T_STOP')
date=fxpar(header, 'DATE_OBS')
xc=fxpar(header, 'CENTER_X')
yc=fxpar(header, 'CENTER_Y')
Radius=fxpar(header, 'R_SUN')
cdx=fxpar(header, 'CDELT1')
cdy=fxpar(header, 'CDELT2')


wls={image:image, nx:nx, ny:ny, date:date, cdx:cdx, cdy:cdy, R:Radius, xc:xc, yc:yc, detected:image}

 hgrm=histogram(Image)
 nh=n_elements(hgrm)
 mx=max(hgrm[1:nh-1])
 qsint0=where(hgrm[1:nh-1] eq mx)
 qsint=qsint0[0]+1
 print, 'Quiet Sun Intensity For The Image:', qsint

return
;	***************************************************************
;		MDI LINK UP
;*******************************************************

; open MDI image

print, fn

;path1='Y:\mdi_data\mg_april_02\fd_M_96m_'

path1='Y:\mdi_data\MDI-magnetograms\mg2002-07\lev1.8\fd_M_96m_'
len=strlen(fn)
;nl=strpos(fn, 'Processed')
nl=strpos(fn, 'PR')


fname=strmid(fn, nl, len-nl)
foldername=strmid(fn, nl-9, 9)


case fname of
	'PR0000.fits': name='0000.fits'
	'PR0001.fits': name='0004.fits'
	'PR0002.fits': name='0007.fits'
	'PR0003.fits': name='0011.fits'
end

;case fname of
;	'Processed0000.fits': name='0000.fits'
;	'Processed0001.fits': name='0004.fits'
;	'Processed0002.fits': name='0007.fits'
;	'Processed0003.fits': name='0011.fits'
;end
fnmdi=path1+foldername+name
print, fnmdi

msgres='Yes'

;	fnmdi=dialog_pickfile(path='Y:\mdi_data\mg_april_02')
	imdi=readfits(fnmdi, mheader)
	GetHeaderInfo, mheader, mnx, mny, mdate, mxc, myc, mR


if anytim(mdate)-anytim(date) ge 60 then  begin
back:					msgres=dialog_message(['Magnetogram Date:'+mdate, $
						'White Light Image Date:'+date, $
						'Bad Date, Choose Another One?'], /question)
					if msgres eq 'Yes' then begin
						fnmdi=dialog_pickfile(path='Y:\mdi_data\mg_april_02')
						imdi=readfits(fnmdi, mheader)
						GetHeaderInfo, mheader, mnx, mny, mdate, mxc, myc, mR
						if anytim(mdate)-anytim(date) gt 60 then goto, back
						endif
					endif



print, 'Magnetogram Info:', mdate

bimdi=bytscl(imdi, max=1000, min=-1000)

end