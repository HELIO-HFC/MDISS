pro efr_CleanImage, im, result,  delta=delta, subtract=subtract, noplot=noplot

;  what is delta??  -- seems like a fudge factor
if not keyword_set(delta) then delta=0.0 else delta=fix(delta)
;print, delta

;  Normalising factor for division in getcontrastimage
norm_fact = 100		; if not specified, set to QuietSun maximum intensity 
 
;  Parameters to standardize the image to 1024 by 1024 pixels, 
;  making Sun's disk circular of radius 420 pixels, centred at coordinates 511.5, 511.5 
nx = 1024 & ny = 1024
rsun = 420.0 & xcen = 511.5 & ycen = 511.5

;
;  Added required library routines using linkimage
;------------------------------------------------------------------------------------------

efr_cleaning_startup

;
;    Correct any distortions in the shape of the solar disk
;------------------------------------------------------------------------------------------

;  Fit ellipse to the input image (im), results of the fits are xc, yc, R1, R2, theta 
;						(others are optional, see the code)

limb_efit, im, xc, yc, R1, R2, theta

;  Standardize the image using the fit parameters and standardization parameters as input

imst = frclean_standardize(im, xc, yc, R1, R2, theta, nx, ny, rsun+delta, xcen, ycen)

;
;    Produce a cleaned image by using the quiet sun image to "correct" the original image
;------------------------------------------------------------------------------------------

;  Generate the quiet sun image based on standardized image

Qsim = qsmedian(imst, rsun, xcen, ycen)

;  Then clean the image using the quiet sun image, either by division or subtraction

if keyword_set(subtract) then result = getcontrastimage(imst, qsim, rsun, xcen, ycen, /subtract) $
    else result = getcontrastimage(imst, qsim, rsun, xcen, ycen, Inorm=norm_fact)  
;OR result=getcontrastimage(imst, qsim, 420., 512, 512.)


;  plot the resulting image...
if not keyword_set(noplot) then begin
  window, 0, xs=1024, ys=1024, title='Cleaned Image'
  tvscl, result
endif

end
