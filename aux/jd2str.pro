;+
; NAME:
;       JD2STR
;
; PURPOSE:
;       Compute string format time from Julian day.
;
; CATEGORY:
;		Time
;
; CALLING SEQUENCE:
;       str = jd2str(jd)
;
; INPUTS:
;		jd - Julian day.
;
; OPTIONAL INPUTS:
;		format - Scalar of integer type specifying the format of the input time:
;					0 - str = 'YYYY-MM-DD HH:NN:SS' (default)
;					1 - str = 'YYYY-MM-DDTHH:NN:SS'
;				    2 - str = 'YYYY-MM-DDTHH:NN:SS.SSS'
;					3 - str = 'YYYY-MM-DD HH:MM:SS.SSS'
;
; KEYWORD PARAMETERS:
;		None.
;
; OUTPUTS:
;       str - Time in the string format ['YYYY-MM-DD HH:NN:SS']  
;
; OPTIONAL OUTPUTS:
;		None.
;
; COMMON BLOCKS:
;		None.
;
; RESTRICTIONS/COMMENTS:
;		None.
;
; MODIFICATION HISTORY:
;		Written by X.Bonnin,	21-JUL-2010.       
;
;-
 
FUNCTION jd2str, jd, format=format

if (n_params() lt 1) then begin
	message,/INFO,'Call is:'
	print,'Result = jd2str(jd,format=format)'
	return,''
endif

if (~keyword_set(format)) then form=0 else form = fix(format[0])
 
n = n_elements(jd)  
str = strarr(n)
 
for i=0l,n-1l do begin 
	caldat,jd[i],mm,dd,yy,hh,nn,ss 
 
 	sec = string(ss,format='(i2.2)')
 	if (form eq 0) or (form eq 2) then sep = 'T' else sep = ' '
 	if (form eq 2) or (form eq 3) then sec = sec + '.' + string((ss-long(ss))*1.e3,format='(i3.3)')
 
	str[i] = string(yy,format='(i4.4)') + '-' + string(mm,format='(i2.2)') + '-' + string(dd,format='(i2.2)') $
		  + sep +string(hh,format='(i2.2)') + ':' + string(nn,format='(i2.2)') + ':' + sec
endfor
if (n eq 1) then str = str[0] 
 
return,str 
END