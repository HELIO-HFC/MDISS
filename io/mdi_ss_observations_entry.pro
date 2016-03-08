;+
; NAME:
;		mdi_ss_observations_entry
;
; PURPOSE:
;   	Returns observations structure entry for observations ASCII file
;
; CATEGORY:
;		I/O
;
; GROUP:
;		MDISS
;
; CALLING SEQUENCE:
;		IDL>struct_out = mdi_ss_observations_entry(header, struct_in)
;
; INPUTS:
;       header     - Input SOHO/MDI fits file header (structure type).
;       struct_in  - Structure containing meta-data for the observations output file.
;       
;
; OPTIONAL INPUTS:
;       indice - Indice to define the meta-data origin:
;					indice=1 --> MDI Ic header (default)
;					indice=2 --> MDI M  header
;	    
; KEYWORD PARAMETERS:
;       None.
;
; OUTPUTS:
;		struct_out - Return the struct_in structure updated with fields found in the header.				
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
;		anytim2jd
;		tim2carr
;       tag_exist
;
; EXAMPLE:
;		None.
;		
; MODIFICATION HISTORY:
;		Written by:		S.Zharkov (MSSL).
;
;       20-NOV-2011, X.Bonnin:  Modify inputs/outputs arguments.
;-

FUNCTION mdi_ss_observations_entry, header, struct_in, indice

if (n_params() lt 2) then begin
    message,/INFO,'Call is:'
    print,'struct_out = mdi_ss_observations_entry(header, struct_in[, indice])'
    return,struct_in    
endif

struct_out = struct_in
if (not keyword_set(indice)) then indice = 1

if (indice eq 2) then begin
	;M
	if (tag_exist(header,'T_OBS')) then date_obs = header.t_obs $
	else message,'No tag T_OBS found in the fits file header!'
endif else begin
	;Ic
	if (tag_exist(header,'DATE_OBS')) then date_obs = header.date_obs $
	else message,'No tag DATE_OBS found in the fits file header!'
endelse

struct_out.date_obs = date_obs
struct_out.date_end = date_obs
jd=anytim2jd(date_obs)
CRnum=fix(tim2carr(date_obs, /dc))
struct_out.jdint = jd.int
struct_out.jdfrac = jd.frac
struct_out.c_rotation = CRnum

if (tag_exist(header,'NAXIS1')) then struct_out.naxis1 = header.naxis1 
if (tag_exist(header,'NAXIS2')) then struct_out.naxis2 = header.naxis2
if (tag_exist(header,'CDELT1')) then struct_out.cdelt1 = header.cdelt1     
if (tag_exist(header,'CDELT2')) then struct_out.cdelt2 = header.cdelt2 
if (tag_exist(header,'X0')) then struct_out.center_x = header.x0
if (tag_exist(header,'Y0')) then struct_out.center_y = header.y0
if (tag_exist(header,'R_SUN')) then struct_out.r_sun = header.r_sun
if (tag_exist(header,'BSCALE')) then struct_out.bscale = header.bscale
if (tag_exist(header,'BZERO')) then struct_out.bzero = header.bzero
if (tag_exist(header,'COMMENT')) then begin
    comment = strjoin(header.comment)
    struct_out.comment = strjoin(strsplit(comment,';',/EXTRACT),',')
endif
if (tag_exist(header,'QUALITY')) then struct_out.quality = header.quality
if (tag_exist(header,'BITPIX')) then struct_out.bitpix = header.bitpix

return,struct_out
end
