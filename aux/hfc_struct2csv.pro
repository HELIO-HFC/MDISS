FUNCTION hfc_struct2csv,struct, $
		     		    header=header,$
                        separator=separator, $
					    error=error

;+
; NAME:
;		hfc_struct2csv
;
; PURPOSE:
; 		Convert the input structure into
;		a scalar or array of string type,
;		where values are separated using 
;		semi colon (;).
;
; CATEGORY:
;		I/O
;
; GROUP:
;		HFC
;
; CALLING SEQUENCE:
;		IDL>string_out = hfc_struct2csv(struct) 
;
; INPUTS:
;		struct - structure to convert.		
;	
; OPTIONAL INPUTS:
;		separator - Specify a separator character between fields.
;                   Default is ";".
;
; KEYWORD PARAMETERS:
;		None.
;
; OUTPUTS:
;		string_out - Scalar or array of string type containing
;					 the structure values separated by semi-colon
;					 characters.
;		
; OPTIONAL OUTPUTS:
;		header - Scalar of string type containing 
;				 the name of the input structure tags.
;				 (As for string_out, names are separated by semi-colons.)  
;		error  - Equal to 1 if an error occurs, 0 else.
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
;		None.
;
; EXAMPLE:
;		None.		
;
; MODIFICATION HISTORY:
;		Written by X.Bonnin, 28-OCT-2011.
;									
;-


error = 1
if (n_params() lt 1) then begin
	message,/INFO,'Call is:'
	print,'string_out = hfc_struct2csv(struct,header=header,separator=separator,error=error)'
	return,'' 
endif

if not (keyword_set(separator)) then sep = ';' else sep = strtrim(separator[0],2)

header = strjoin(strupcase(tag_names(struct)),sep)

nrow = n_elements(struct)
ntags = n_tags(struct)
tags = tag_names(struct)

string_out = strarr(nrow)
for i=0l,nrow-1l do begin
	string_out_i = ''
	for j=0,ntags-1 do begin
        flag = execute('struct_i = strtrim(struct(i).'+tags(j)+',2)')
        if (struct_i eq "") then struct_i = "NULL"
		string_out_i = string_out_i + struct_i + sep
    endfor
	string_out[i] = strmid(string_out_i,0,strlen(string_out_i)-1)
endfor

error = 0
return,string_out
END