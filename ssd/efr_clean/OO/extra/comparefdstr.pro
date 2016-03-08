function CompareFDStr, st01, st02

result=2

if (st01.nx ne st02.ny) $
	or (st01.ny ne st02.ny) $
	or abs(st01.xc-st02.xc) gt .1 $
	or abs(st01.xc-st02.xc) gt .1 $
	or abs(st01.R-st02.R) gt .1   $
	or abs(st01.cdx-st02.cdx) gt .1 $
	or abs(st01.cdy-st02.cdy) gt .1 then result=0 $
else $
	if abs(anytim(st01.date)-anytim(st02.date)) gt 3*600 then result=1

return, result

end