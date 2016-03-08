FUNCTION mdi_id2date,index

;+
; NAME:
;		mdi_id2date
;
; PURPOSE:
;   	Convert a SOHO/MDI file index number into its corresponding date of observation.
;
; CATEGORY:
;		I/O
;
; GROUP:
;		None.
;
; CALLING SEQUENCE:
;		IDL>date = mdi_id2date(num)
;
; INPUTS:
;       index - SOHO/MDI files index
;
; OPTIONAL INPUTS:
;       None.
;
; KEYWORD PARAMETERS:
;      	None.
;
; OUTPUTS:
;		date - Corresponding SOHO/MDI date of observation (date format is 'YYYY-MM-DD').				
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
;		print,mdi_id2date(1461)
;               '1997-01-01'
;		
; MODIFICATION HISTORY:
;		Written by:		X.Bonnin (LESIA).
;
;-



if (n_params() lt 1) then begin
    message,/INFO,'Call is:'
    print,'date = mdi_id2date,num'
endif

index0 = 1234
date0 = '1996-05-19'
jd0 = julday(05,19,1996,12,00,00) 

jd = jd0 + (index[0] - index0)
caldat,jd,mm,dd,yyyy,hh,nn,ss
date = string(yyyy,format='(i4.4)')+'-'+string(mm,format='(i2.2)')+'-'+string(dd,format='(i2.2)')

return,date
END