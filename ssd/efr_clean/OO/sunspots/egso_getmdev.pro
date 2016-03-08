;===================================================================

function egso_getmdev, x

; function
;	 to return abs mean deviation
;			(deals with one element arrays
;

if n_elements(x) le 1 then md=0 else m=moment(x, mdev=md)

return, md

end
