FUNCTION SunspotObs::Init, sp, wl

spst=sp[0]->getstructure()
wlst=wl->getstructure()
if  comparefdstr(spst, wlst) ne 2 then return, 0

n=n_elements(sp)
self.spotArr=ptr_new(sp)
struct_assign, wl, self, /nozero

return, 1
END

;-------------------------------------------------------------

PRO SunspotObs::SpotDisplay, n

sp=*self.spotArr

sp[n]->cropdisplay, self, scale=10

END

;-------------------------------------------------------------

FUNCTION SunspotObs::N_spots

sp=*self.spotArr

return, n_elements(sp)

END

;-------------------------------------------------------------

FUNCTION SunspotObs::date



return, self.date

END
;-------------------------------------------------------------

FUNCTION SunspotObs::GetSpot, n

sp=*self.spotArr

return, sp[n]

END

;-------------------------------------------------------------

FUNCTION SunspotObs::GetImage


return, *self.image

END


;--------------------------------------------------------------

FUNCTION SunspotObs::DataArray


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
	endfor
return, data

END


;================================================================

PRO SunspotObs::DisplayALL, show_on=show_on, umbra=umbra

n=n_elements(*self.spotarr)

image=intarr(self.nx, self.ny)
;print, n
for i=0, n-1 do begin

	locs=(*self.spotarr)[i]->getlocs()
	ulocs=(*self.spotarr)[i]->getumbra()
	image[locs]=1
	if ulocs[0] ne -1 then image[locs[ulocs]]=2

endfor

if not arg_present(show_on) then tvscl, image $

	else	begin

			locs=where(image eq 1)
			ulocs=where(image eq 2)
			im=show_on->getdisplay()
			im[locs]=max(im)
			im[ulocs]=min(im)
			tvscl, im


endelse
;ptr_free, self.umbra
;ptr_free, self.locs

END

;================================================================

PRO SunspotObs::Cleanup

;ptr_free, self.umbra
;ptr_free, self.locs

END
;----------------------------------------------------------------------------
; SunspotObs__DEFINE
;
; Purpose:
;  Defines the object structure for an orb object.
;
PRO SunspotObs__define

    struct = { SunspotObs, $
    			 spotArr: ptr_new(/allocate_heap), $
    			 inherits FullDiskObs $
             }
END