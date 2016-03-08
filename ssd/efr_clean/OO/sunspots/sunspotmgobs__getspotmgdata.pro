;-------------------------------------------------------------

FUNCTION SunspotMgObs::GetSpotMgData, n

;	FUNCTION
;		returns the array of magnetic data
;		extracted from the n-th sunspot
;
;	data[0]		Absolute Total Flux
;	data[1]		Spot Max Flux
;	data[2]		Spot Min Flux
;	data[3]		Spot Total Flux
;	data[4]		Absolute Umbral Flux
;	data[5]		Umbral Max Flux
;	data[6]		Umbral Min Flux
;	data[7]		Umbral flux
;
if n gt n_elements(sp) then return, -1

data=dblarr(8)

sp=(*self.spotarr)[n]



return, *self.image

END
