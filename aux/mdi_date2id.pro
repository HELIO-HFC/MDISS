FUNCTION mdi_date2id,date

;+
; NAME:
;		mdi_date2id
;
; PURPOSE:
;   	Convert SOHO/MDI date of observation into its corresponding file index.
;
; CATEGORY:
;		I/O
;
; GROUP:
;		None.
;
; CALLING SEQUENCE:
;		IDL>num = mdi_date2id(date)
;
; INPUTS:
;       date - SOHO/MDI date of observation (scalar of string type).
;              Date format is "YYYY-MM-DD".
;
; OPTIONAL INPUTS:
;       None.
;
; KEYWORD PARAMETERS:
;      	None.
;
; OUTPUTS:
;		index - SOHO/MDI file index				
;
; OPTIONAL OUTPUTS:
;		None.		
;		
; COMMON BLOCKS:
;		None.		
;
; SIDE EFFECTS:
;		None.
;
; RESTRICTIONS/COMMENTS:
;		None.
;
; CALL:
;       None.
;
; EXAMPLE:
;		print,mdi_date2id('1997-01-01')
;               1461
;		
; MODIFICATION HISTORY:
;		Written by:		X.Bonnin (LESIA).
;
;-


if (n_params() lt 1) then begin
    message,/INFO,'Call is:'
    print,'date = mdi_date2id,date'
endif

index0 = 1234
date0 = '1996-05-19'
jd0 = julday(05,19,1996,12,00,00) 

jd = fix(strsplit(date[0],'-',/EXTRACT))
jd = julday(jd[1],jd[2],jd[0],12,0,0)

index = long((jd - jd0) + index0)

return,index
END