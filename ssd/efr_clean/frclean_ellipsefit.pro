;================================================================

pro frclean_ellipsefit, binary, xc, yc, xc1, yc1, R1, R2, theta, stdev, stdevGeo, algerr

binary_copy=binary
vector=fitellipselsq(binary,  xc, yc)
print, 'first fit'

frclean_ealgdata, xc, yc, vector, xc1, yc1, R1, R2, theta

count=0
stdev=1.
while stdev ge .3e-02 and count lt 200 do begin
		count=count+1

		vector=fitellipselsq(binary,  xc, yc)
		frclean_efiterr, binary, xc, yc, vector, binary2, stdev
		binary=binary2
endwhile

		vector=fitellipselsq(binary,  xc, yc)
		frclean_ealgdata, xc, yc, vector, xc1, yc1, R1, R2, theta
		color=max(binary)/2


frclean_efiterr, binary, xc, yc, vector, binary2, stdev, stdevGeo, algerr, geoErr

;stdev=stdev
;stdevGeo=stdevGeo
;algerr=algerr


end
;=============================================================================
