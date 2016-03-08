;==============================================================================
PRO Sunspot::GetGC, arcx, arcy, BClon, BClat, diameter=diameter

arcx=self.GCarcx
arcy=self.GCarcy

;print, self.date
clon=tim2carr(self.date)

case self.origin of

'soho':		GChel=arcmin2hel(self.GCarcx/60, self.GCarcy/60, date=self.date, /soho)
'ground':	GChel=arcmin2hel(self.GCarcx/60, self.GCarcy/60, date=self.date)

endcase

Bclon = GChel[1]+clon[0]
Bclat= GChel[0]
diameter=self.helDiam

END
