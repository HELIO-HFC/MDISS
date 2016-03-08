;-----------------------------------------------------------------------------

FUNCTION Sunspot:: GetSCounts, n_umb=n_umb, ucount=ucount, area=area

; Returns SampleCount i.e. # of pixels
;	also umbral pixel size and #umbras

if (*self.umbra)[0] ne -1 then $
		ucount=n_elements(*self.umbra) $
		else ucount=0
n_umb=self.n_umbras
area=self.helArea

return, n_elements(*self.locs)

end