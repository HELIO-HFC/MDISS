FUNCTION SunspotMgObs::Init, wl, mg, sp, no_show=no_show, one=one, $
			mgdat=mgdat, res=res, fn=fn, local=local, drot=drot, $
			show_dates=show_dates,  missvalscheck=missvalscheck, time_dif=time_dif, $
			checkmg=checkmg, year=year, obsmg=obsmg, getobs=getobs, setmgvalues=setmgvalues

;	initiate sunspot obserevation object
;		containing original wl observation,
;		array of sunspots detected,
;		and synchronised magnetogram data
;
;	KEYWORDS
;
;		drot			use differential rotation when synchronising magnetogram
;
;		one				use "overdetect" option for WL detection
;						preferred option when using magnetogram data
;						to validate detection results

;stop
if keyword_set(local) then path='C:\SERHII\MDI\continuum' $
	else path='Y:\ScratchSpace\DATA\MDI\continuum'
if keyword_set(local) then mgpath='C:\SERHII\MDI\magnetograms' $
	else mgpath='Y:\ScratchSpace\DATA\MDI\magnetograms'


if n_params() eq 0 $
	and not keyword_set(fn) then wl=obj_new('wlfitsfdobs'); , map=wlmap)

if keyword_set(fn) then wl=obj_new('wlfitsfdobs', location=fn)

if n_params() le 1 or keyword_set(fn) then begin
	if not obj_isa(wl, 'wlfitsfdobs') then return, -1


date=wl->getdate()

;stop
if keyword_set(show_dates) then $
	print, 'WhiteLightDate:', date

	if keyword_set(mgdat) then begin
		self.mgdata=ptr_new(mgdat)
		n=n_elements(sp)
		self.spotArr=ptr_new(sp)
		struct_assign, wl, self, /nozero
		return, 1
	endif

;fn2=MgFind( date, obsmg=obsmg, lac=lac, local=local, year=year, getobs=getobs)

; if not keyword_set(lac) then $
; 	mg=obj_new('mgfitsfdobs', location=mgpath+fn2[0]) $
;	else mg=obj_new('mgfitsfdobs', location=fn2[0])

mgdate=mg->getdate()
date=wl->getdate()

if keyword_set(show_dates) then $
	print, 'MagnetogramDate:', mgdate, '.       Location:', fn2
endif

mgdate=mg->getdate()
date=wl->getdate()
jd1=anytim2jd(date)
jd2=anytim2jd(mgdate)
time_dif=abs(double(jd1.int)+jd1.frac - $
				double(jd2.int)-jd2.frac)

if (time_dif ge .5) then begin
	if not keyword_set(no_show) then $
		print, 'No magnetogram found within one day of observation'
		return, -1
	endif
;spst=sp[0]->getstructure()
wlst=wl->getstructure()
mgst=mg->getstructure()


;stop
;wl->getspots, sp, no_show=no_show, one=one, res=res

;if  comparefdstr(spst, wlst) ne 2 then return, 0

;wlmap=wl -> getmap()
;mgmap=mg -> getmap()

;stop
if  comparefdstr(mgst, wlst) ne 2 then $
	if not keyword_set(drot) then rotmg= fr_sync_soho_obs(wl, mg) $
		else begin
			wlmap=wl -> getmap(/test, /soho)
			mgmap=mg -> getmap(/test, /soho)
	;		stop
			dmgm=drot_map(mgmap,ref_map=wlmap)
			rotmg=dmgm.data
		endelse $
	else rotmg=mg->getimage()

;stop

self.mgdata=ptr_new(rotmg)
;n=n_elements(sp)
;self.spotArr=ptr_new(sp)
struct_assign, wl, self, /nozero

self.image=ptr_new(wl->getimage())

self->sd, no_show=no_show,  missvalscheck=missvalscheck, $
		checkmg=checkmg, one=one

if keyword_set(setmgvalues)  then begin
    mgim=*self.mgdata 
  for i=0, self.n_spots-1 do begin
      (*self.spotarr)(i)->setmagneticflux, mgim
     ; stop
   endfor
endif


; *** cleaning bits - added 2005-07-01
;stop
if n_params() eq 0 $
	or keyword_set(fn) then begin
	obj_destroy, wl
	obj_destroy, mg
endif

return, 1

END



;-------------------------------------------------------------

PRO SunspotMgObs::SpotDisplay, n, image=image

sp=*self.spotArr

print, 'aaa'
stop
if (sp[0] ne -1) and (n lt n_elements(sp)) then begin
	sp[n]->cropdisplay, self, scale=10
	if keyword_set(image) then begin
			locs=sp[n]->getlocs()
			imcopy=image
			imcopy[locs]=2*imcopy[locs]
			window, 1, xs=800, ys=600; xs=self.nx, ys=self.ny
			tvframe, imcopy, /asp ;tvscl, imcopy
		endif
 endif
END


;-------------------------------------------------------------

FUNCTION SunspotMgObs::date



return, self.date

END

;-------------------------------------------------------------

FUNCTION SunspotMgObs::GetImage


return, *self.image

END


;--------------------------------------------------------------

FUNCTION SunspotMgObs::DataArray


sp=*self.spotArr
n=n_elements(sp)

Data=dblarr(n, 14)
	for i=0, n-1 do begin
		sp[i]->SpotInfo, Bcx, Bcy, Bclon, Bclat, $
									n_umb, pixSize, upixSize, helArea, $
									totFlux, umbFlux, meanInt, QuietSunInt, $
									 maxFlux, maxuFlux,  /noprint
			data[i, 0]=Bcx
			data[i, 1]=Bcy
			data[i, 2]=Bclon
			data[i, 3]=Bclat
			data[i, 4]=n_umb
			data[i, 5]=pixSize
			data[i, 6]=upixSize
			data[i, 7]=helArea
			data[i, 8]=totFlux
			data[i, 9]=umbFlux
			data[i, 10]=maxFlux
			data[i, 11]=maxuFlux
			data[i, 12]=meanInt
			data[i, 13]=QuietSunInt
			;data[i, 14]=QuietSunInt
	endfor
return, data

END


;================================================================

PRO SunspotMgObs::DisplayALL, show_on=show_on, umbra=umbra, bin=bin, over=over

n=n_elements(*self.spotarr)

image=intarr(self.nx, self.ny)

;print, n
for i=0, n-1 do begin

	locs=(*self.spotarr)[i]->getlocs()
	ulocs=(*self.spotarr)[i]->getumbra()
	image[locs]=1
	if ulocs[0] ne -1 then image[locs[ulocs]]=2

endfor
bin=image
if not arg_present(show_on) then tvframe, image, /asp $ ;tvscl, image

	else	begin

			locs=where(image eq 1)
			ulocs=where(image eq 2)
			im=show_on->getdisplay()
			im[locs]=min(im)
			im[ulocs]=max(im)
			tvframe,  im, /asp ;tvscl, im


endelse
;ptr_free, self.umbra
;ptr_free, self.locs

END

;================================================================

PRO SunspotMgObs::Cleanup

;stop
if self.n_spots ne 0 then begin
;	n=n_elements(*self.spotarr)
	n=self->n_spots()
	for i=0, n-1 do if obj_valid((*self.spotarr)(i)) then obj_destroy, (*self.spotarr)(i)
 endif else begin
 	ptr_free, self.spotarr
;	 stop
 endelse
ptr_free, self.spotArr
ptr_free, self.mgData
ptr_free, self.header
ptr_free, self.image
;stop
;ptr_free, self.locs

END
;----------------------------------------------------------------------------
; SunspotMgObs__DEFINE
;
; Purpose:
;  Defines the object structure for an orb object.
;
PRO SunspotMgObs__define

    struct = { SunspotMgObs, $
    			 n_spots:0, $
    			 spotArr: ptr_new(0,/allocate_heap), $
    			 mgdata: ptr_new(/allocate_heap), $
    			 inherits wlfitsfdobs $
             }
END
