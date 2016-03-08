;-------------------------------------------------------------

FUNCTION SunspotMgObs::N_spots

;sp=*self.spotArr

;if (sp[0] ne -1) then rs=n_elements(sp) else rs=0
return, self.n_spots

END
