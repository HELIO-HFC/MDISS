;+
; NAME:
;		printl
; 
; PURPOSE:
;		Print a message on the terminal screen without skipping line. 
; 
; INPUTS:
;		text = Scalar string wich contains the message to print.
;
; KEYWORD PARAMETERS:
;		None.
;
; OUTPUTS:
;		None.
;
; HISTORY:
;		Written by X.Bonnin, 10-MAR-2006
;		06-MAY-2010, X.Bonnin:	Renamed to printl.
;
;-



PRO printl, text

if (n_params() eq 1) then begin
	CASE !version.os_family OF
    	'Windows': cr = string("15b)+string("12b)
    	'MacOS' : cr = string("15b)
    	'unix' : cr = string("15b)
	ENDCASE
	
	txt = string(text[0])
	form="($,a,a)" 
	print,form = form,txt,cr
endif else print,'Call is: printl,text'

END