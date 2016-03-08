pro GetTrackData, sobsarr, spnx, n, offset=offset

flux=dblarr(n)
helSize=dblarr(n)
uhelSize=dblarr(n)
maxflux=dblarr(n)
irradiance=dblarr(n)

timex=dblarr(n, 7)
time=dblarr(n)
time2=dblarr(n)
for i=0, n-1 do begin

		dt=sobsArr[i+offset]->date()
		s=anytim2jd(dt)
		time2[i]=s.int+s.frac
		timex[i, *]=anytim(dt, /ex)
		time[i]=anytim(sobsArr[i+offset]->date())
		sp=sobsArr[i+offset]->getspot(spnx[i])

		sp->SpotInfo, Bcx, Bcy, Bclon, Bclat, $
									n_umb, pixSize, upixSize, helArea, $
									totFlux, umbFlux, meanInt, QuietSunInt, $
									 maxFlux0, maxuFlux,  /noprint

		if Bcx gt 700 then break
		maxflux[i]=maxFlux0
		flux[i]=totFlux
		helSize[i]=helArea
		uhelSize[i]=double(uPixSize)*helArea/pixSize
		irradiance[i]=double(meanInt)/QuietSunInt

;		data=sobsArr[i]->dataarray()
;		maxflux[i]=data[spnx[i], 10]
;		flux[i]=data[spnx[i], 8]
;		helSize[i]=data[spnx[i], 7]
;		uhelSize[i]=double(data[spnx[i], 6])*data[spnx[i], 7]/data[spnx[i], 5]
;		irradiance[i]=double(data[spnx[i], 12])/data[spnx[i], 13]
end

;tm=julday(timex[0, 5], timex[0, 4], $
;							timex[0, 6], timex[0, 0], $
;							timex[0, 1], timex[0, 2])
print, timex[0, *]

print,i
dummy2 = LABEL_DATE(DATE_FORMAT=['%M %D, %H %A'] , offset=0)

mytime = TIMEGEN(START=julday(timex[0, 5], timex[0, 4], $
							timex[0, 6], timex[0, 0], $
							timex[0, 1], timex[0, 2]), $
			FINAL=julday	(timex[i-1, 5], timex[i-1, 4], $
							timex[i-1, 6], timex[i-1, 0], $
							timex[i-1, 1], timex[i-1, 2]));, $
;res=label_date(DATE_FORMAT=time2[*])
;print, time2
window, 10, xs=1024, ys=700
plot, time2[0:i-1], abs(flux[0:i-1]), yticks=1, xstyle=1, XTICKFORMAT = 'LABEL_DATE'
oplot, time2[0:i-1], abs(maxflux[0:i-1])*100 , linestyle=1
oplot, time2[0:i-1], helSize[0:i-1]*10000, linestyle=2
oplot, time2[0:i-1], uhelSize[0:i-1]*100000l, linestyle=3
;oplot, time, irradiance[0:i-1]*250000l, linestyle=4


;window, 11, xs=1024, ys=700
;plot, time, abs(flux[0:i-1]), xtitle='Time', yticks=1, xstyle=1;, XTICKFORMAT = 'LABEL_DATE'
;oplot, time, abs(maxflux[0:i-1])*100 , linestyle=1
;oplot, time, helSize[0:i-1]*10000, linestyle=2
;oplot, time, uhelSize[0:i-1]*100000l, linestyle=3

;stop
end



