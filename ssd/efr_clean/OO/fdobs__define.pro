;----------------------------------------------------------------------------
; FDObs__DEFINE
;
; Purpose:
;  Defines the object structure for an orb object.
;
PRO FDObs__define

    struct = { FDObs, $
    			R:0.d00,	$
    			xc:0.d00,	$
    			yc:0.d00,	$
    			cdx:0.d00,	$
    			cdy:0.d00,	$
    			nx:0,		$
    			ny:0,		$
    			date:'',	$
    			origin:''	$
             }
END