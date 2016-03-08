;-----------------------------------------------------------------

PRO Sunspot::SpotInfo, Bcx, Bcy, Bclon, Bclat, $
						n_umb, pixSize, upixSize, helArea, $
						totFlux, umbFlux, meanInt, QuietSunInt, maxFlux, maxuFlux,$
						 noprint=noprint, lun=lun, diam=diam

clon=tim2carr(self.date)

case self.origin of

'soho':		GChel=arcmin2hel(self.GCarcx/60, self.GCarcy/60, date=self.date, /soho)
'ground':	GChel=arcmin2hel(self.GCarcx/60, self.GCarcy/60, date=self.date)

endcase


diam=self.heldiam
umb=*self.umbra
if umb[0] eq -1 then UmbSize=0 else UmbSize=n_elements(umb)

Bcx= self.GCarcx
Bcy= self.GCarcy
Bclon = GChel[1]+clon[0]
Bclat= GChel[0]

pixSize=n_elements(*self.locs)
upixSize=UmbSize
n_umb= self.n_umbras
helArea= self.helArea
totFlux= self.totFlux
umbFlux= self.umbFlux
meanInt=self.meanInt
QuietSunInt= self.QuietSunInt
maxFlux=self.maxFlux
maxuFlux=self.maxuFlux


if not keyword_set(noprint) then begin
	;	STOP
		if not keyword_set(lun) then begin
			print, 'Barrycenter ArcCoodinates:', self.GCarcx, self.GCarcy
			print, 'Barrycenter HelioCoodinates:', GChel[1]+clon[0], $
								 						GChel[0]
			print, 'Original Pixel Size', n_elements(*self.locs)
			print, 'Umbra Pixel Size:', UmbSize
			print, 'Number Of Umbras:', self.n_umbras
			print, 'Area:', self.helArea, '    square degrees'
			print, 'Diameter:', self.helDiam, '     degrees'

			print, 'Total Flux:', self.totFlux
			print, 'Umbral Flux:', self.umbFlux
			print, 'Mean Intensity', self.meanInt
			print, 'QuietSunIntensity', self.QuietSunInt
			print, 'Max Flux:', self.maxFlux
			print, 'Max Umbral Flux:', self.maxuFlux
		endif else begin
			printf, lun, 'Barrycenter ArcCoodinates:', self.GCarcx, self.GCarcy
			printf, lun, 'Barrycenter HelioCoodinates:', GChel[1]+clon[0], $
								 						GChel[0]
			printf, lun, 'Original Pixel Size', n_elements(*self.locs)
			printf, lun, 'Umbra Pixel Size:', UmbSize
			printf, lun, 'Number Of Umbras:', self.n_umbras
			printf, lun, 'Area:', self.helArea, '    square degrees'

			printf, lun, 'Total Flux:', self.totFlux
			printf, lun, 'Umbral Flux:', self.umbFlux
			printf, lun, 'Mean Intensity', self.meanInt
			printf, lun, 'QuietSunIntensity', self.QuietSunInt
			printf, lun, 'Max Flux:', self.maxFlux
			printf, lun, 'Max Umbral Flux:', self.maxuFlux
		 	endelse
	endif
END