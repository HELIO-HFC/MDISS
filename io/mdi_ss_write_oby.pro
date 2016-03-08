;+
; NAME:
;		mdi_ss_write_oby
;
; PURPOSE:
; 		Write Observatory information in a csv format file.
;
; CATEGORY:
;		I/O
;
; GROUP:
;		SSMDI
;
; CALLING SEQUENCE:
;		mdi_ss_write_oby,fnoby
;
; INPUTS:
;		fnoby - Scalar of string type containing the full path for the output ascii file.	
;	
; OPTIONAL INPUTS:
;		run_date - Scalar of string type containing the current date and time in ISO 8601 format.
;
; KEYWORD PARAMETERS:
;		/SILENT - Quiet mode.	
;
; OUTPUTS:
;		None.	
;
; OPTIONAL OUTPUTS:
;		error - Equal to 1 if an error occurs, 0 otherwise.
;
; COMMON BLOCKS:		
;		None.
;	
; SIDE EFFECTS:
;		None.
;
; RESTRICTIONS/COMMENTS:
;		None
;		
; CALL:
;		anytim
;
; EXAMPLE:
;		None.		
;
; MODIFICATION HISTORY:
;		Written by:		X.Bonnin, 18-OCT-2011.
;
;-


PRO mdi_ss_write_oby,fnoby, $
					 run_date=run_date,$
					 error=error, $
					 SILENT=SILENT

error = 1
if (n_params() lt 1) then begin
	message,/CONT,'Call is:'
	print,'mdi_ss_write_oby,fnoby,error,/SILENT'
	return
endif	
SILENT = keyword_set(SILENT)
if (~keyword_set(run_date)) then run_date = anytim(!stime, /ccsds)

;Observations information
obs_info = {Observatory:'SOHO',Instrument:'MDI',Telescope:['continuum','magnetogram'],$
			Units:['Counts','Gauss'],ChannelID:'SOHO/MDI',MeanWavel:'676.78',$
			WavelUnit:'nm',wavename:'Ni I',spectral_name:['visible','line_of_sight magnetic field'],$
			obs_type:['Remote-sensing'],comment:['Ic_6h','M_96m']}


; *** Write observatory output file for the current session ****
if (~SILENT) then print,'Writing observatory csv file...'

;Write Observatory file
openw,lun,fnoby,/GET_LUN

;Fields: 
printf,lun,'ID_OBSERVATORY;OBSERVAT;INSTRUMENT;TELESCOPE;UNITS;WAVEMIN;WAVEMAX;WAVENAME;WAVEUNIT;SPECTRAL_NAME;OBS_TYPE;COMMENT;RUN_DATE'
;Continuum values
printf,lun,strjoin(['1',obs_info.Observatory,obs_info.Instrument,obs_info.Telescope[0],obs_info.Units[0],$
				   obs_info.MeanWavel,obs_info.MeanWavel,obs_info.wavename,obs_info.WavelUnit,obs_info.spectral_name[0],$
				   obs_info.obs_type,obs_info.comment[0],run_date],';')
;Magnetogram values
printf,lun,strjoin(['2',obs_info.Observatory,obs_info.Instrument,obs_info.Telescope[1],obs_info.Units[1],$
				   obs_info.MeanWavel,obs_info.MeanWavel,obs_info.wavename,obs_info.WavelUnit,obs_info.spectral_name[1],$
				   obs_info.obs_type,obs_info.comment[1],run_date],';')

close,lun
free_lun,lun

if (~SILENT) then print,fnoby+' has been saved.'
if (~SILENT) then print,'Writing observatory csv file OK'

error = 0
END
