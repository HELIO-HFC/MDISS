pro efr_cleaning_startup, try_os=try_os

;	PROCEDURE
;
;	loads gauss_smoothing and canny function into IDL
;
;	NB: currently only tested for windows, MacOS (darwin)

os_arch = !version.OS + '_' + !version.ARCH
if keyword_set(try_os) then os_arch = try_os
message,'Loading EFR library shared objects for OS/Arch:  '+os_arch,/cont

linkdir = concat_dir('$EGSO_FRLIB',os_arch)
files = file_list(linkdir,'*')

routine = 'gscanny.dll'
linkname = concat_dir(linkdir,routine)

;  this probably only true for Windows
if file_exist(linkname) then begin
  print,linkname
  linkimage, 'gauss_smoothing', linkname, 1, 'GaussSmoothing'
  linkimage, 'canny', linkname, 1, 'CannyHT'

;  this may be the more normal case
endif else if n_elements(files) ge 2 then begin
  f_gauss = concat_dir(linkdir,'gauss.so') & print,f_gauss
  if file_exist(f_gauss) then linkimage, 'gauss_smoothing', f_gauss, 1, 'GaussSmoothing'  $
    else message,'Problem adding: '+f_gauss
  f_canny = concat_dir(linkdir,'canny.so') & print,f_canny
  if file_exist(f_canny) then linkimage, 'canny', f_canny, 1, 'CannyHT'  $
    else message,'Problem adding: '+f_gauss

endif else begin
  message, 'Problem adding EGSO Feature Recognition library routine', /cont
  message, 'Shared Object NOT implemented for OS/ARCH:  '+os_arch

endelse

end
