pro getMgDBData, obs

n=n_elements(obs.filename)

path='Y:\mdi_data\MDI-whitelight'

count=1
fname=path+'/SpotsMgJuly.txt'
get_lun, lun
openw, lun, fname
for i=0, n-1 do begin

	 sp=ReadSunspotObs(path+obs.spotlocspath[i])
	 Nsp= n_elements(sp)
	 print, Nsp
;	 wl=obj_new('wlfitsfdobs', location=path+obs.filename[i])
	for j=0, Nsp-1 do begin

		sp[j]->SpotInfo, Bcx, Bcy, Bclon, Bclat, $
						n_umb, pixSize, upixSize, helArea, $
						totFlux, umbFlux, meanInt, QuietSunInt, maxFlux, maxuFlux,$
						 /noprint

		printf, lun, FORMAT= $
		'(I8, ", ", I5, ", ", I5, 2( ",", A65 ), 3(",", I8), 12(",", D15.5))', $;, ",", I10, ",", D0.10,  3(",",I5 ), 5(",", F9.4), 2(",", I6))', $
			 count, i, j, obs.filename[i], obs.spotlocspath[i], pixSize, n_umb, upixsize,$
			 bclat, bclon, bcx, bcy, helarea,  $
			 meanInt/QuietSunInt, totFlux, umbFlux, maxFlux, maxuFlux
		count=count+1

	;	sp[j]->cropdisplay, wl
	;	stop
	endfor
endfor

close, lun
free_lun, lun

end