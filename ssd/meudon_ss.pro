pro meudon_ss, fnc, fnm, fn, display=display

;         PROCEDURE 
;               to run sunspot detection on SOHO/MDI intensity and magnetogram data
;
;         INPUTS
;               fnc        locations of MDI level 2 intensity images (limb darkening removed)
;               fnm        locations of corresponding MDI magnetograms
;
;
;         OUTPUT
;               fn          list of written xml files corresponding to each individual detected sunspot



XSIZE=800
YSIZE=700
;get_mdi_files, t0, t1, files_cont=fnc, files_fd=fnm
fn=''
n=n_elements(fnc)

;stop
if keyword_set(display) then window, 1, xs=1000, ys=800, xp=XSIZE+100, yp=YSIZE
for i=0l, n-1 do begin

 ; **** read in white light and magnetogram data
   wl=obj_new('wlfitsfdobs', loc=fnc(i))
   mg=obj_new('mgfitsfdobs', loc=fnm(i))

  if keyword_set(display) then begin
      window, 14, xs=XSIZE, ys=YSIZE, tit='Continuum'+wl->getdate()
      wl->display
      window, 15, xs=XSIZE, ys=YSIZE, tit='Magnetogram', xp=XSIZE, yp=0
      mg->display
      wset, 1
  endif

 ; **** run the detection
  if keyword_set(display) then  sp=obj_new('sunspotmgobs', wl, mg) $
  else  sp=obj_new('sunspotmgobs', wl, mg, /no_show)


 ; **** write detection result to xml file (LMSAL version)
ss_xml, sp, fn0
 k=sp->n_spots()
 
 ; **** clean objects and pointers
obj_destroy, sp
obj_destroy, wl
obj_destroy, mg
heap_gc;, /verb


if k ne 0 then fn=[fn, fn0]
if keyword_set(display) and (k ne 0) then begin
    if i mod 8 eq 0 then stop else wait, 3.
endif
endfor
if n_elements(fn) gt 1 then fn=fn(1:*)
end
