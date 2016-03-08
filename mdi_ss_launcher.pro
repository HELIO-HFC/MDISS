PRO mdi_ss_launcher,starttime=starttime,endtime=endtime,$
				    sample=sample, $
				    data_dir=data_dir,$
					fnc=fnc, fnm=fnm, $
				    output_dir=output_dir,$
					write_png=write_png, $ 
				    WRITE_CSV=WRITE_CSV,CLEAN_DATA=CLEAN_DATA,$
				    DOWNLOAD_DATA=DOWNLOAD_DATA,$
				    SILENT=SILENT,UNPROCESSED=UNPROCESSED,$
                    NO_SHOW=NO_SHOW,HELP=HELP

;+
; NAME:
;		mdi_ss_launcher
;
; PURPOSE:
;		Launch the SunSpot recogniton code for SOHO/MDI data 
;		(developped by Serguei Zharkov in the frame of the 
;		EGSO project) in a given time range.
;		Data can be loaded from local disk, or 
;		downloaded from MEDOC server.
;
; CATEGORY:
;		Image processing
;
; GROUP:
;		MDISS
;
; CALLING SEQUENCE:
;		IDL> mdiss_ss_launcher,starttime=starttime,endtime=endtime
;
; INPUTS:
;		starttime - Scalar of string type containing the starting date 
;				    of SOHO/MDI data to process.
;                   Date format is "YYYY-MM-DD".
;		endtime	  - Scalar of string type containing the ending date 
;				    of SOHO/MDI data to process.
;                   Date format is "YYYY-MM-DD".
;
; OPTIONAL INPUTS:
;		sample		- Define the cadence of SOHO/MDI images 
;					  (by default highest cadence of Ic images).		
;		data_dir 	- The full path name of the directory where the
;					  SOHO/MDI fits files are saved. 
;				      Default is current one.
;		output_dir	- The full path name of directory where the output files produced by SDOSS code will be saved.
;					  Default is current one.
;		fnc			- Vector of string type providing the list of 
;					  input SOHO/MDI Ic fits file(s) to process. (Must be fd_Ic_6h_01d fits file(s)).
;		fnm			- Vector of string type providing the list of 
;					  input SOHO/MDI M fits file(s) to process. (Must be fd_M_96m_01d fits file(s)).
;		write_png   - Write output png image files:
;						write_png = 0 --> no output png file. (Default.)
;						write_png = 1 --> write only pre-processed continuum images (without detection results).
;						write_png = 2 --> Like 1, but also write the corresponding magnetograms.
;						write_png = 3 --> write pre-processed continuum images with detection results.
;						write_png = 4 --> Like 1, but also write the corresponding magnetograms with detection results.
;
; KEYWORD PARAMETERS: 
;       /UNPROCESSED    - Run pre-processing on input data.
;		/CLEAN_DATA	    - Remove SOHO+MDI data file from local disk after processing.
;						  (USE WITH CAUTION!)
;		/WRITE_CSV	    - Write detection results into csv format files 
;		/DOWNLOAD_DATA  - If set, download data from MEDOC server if no file
;						  is found in the local disk.
;						  Downloaded files will be saved in data_dir directory.
;		/SILENT		    - Quiet mode.
;       /NO_SHOW        - No display observations and results.
;       /HELP           - Display help.
;		
; OUTPUTS:
;		None.				
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
;		- The Solar SoftWare (SSW) must be loaded.
;		- MDISS auxiliary IDL routines must be compiled.
;		- An internet access is required to download SOHO/MDI data from MEDOC.
;
;
; CALL:
;		mdi_ss_finding_observations
;		mdi_date2id
;       mdi_wget
;		jd2str
;		tim2jd
;       anytim
;
; EXAMPLE:
;		None.
;		
; MODIFICATION HISTORY:
;		Written by:		X.Bonnin, 10-JAN-2012. 
;
;-

;[1]:Initializing 
;[1]:============
launch_time = systime(/sec)
run_date = !stime
if (keyword_set(HELP)) then begin
	message,/INFO,'Call is:'
	print,'IDL>mdi_ss_launcher,starttime=starttime,endtime=endtime, $'
	print,'                    sample=sample,$'
	print,'                    data_dir=data_dir,$'
	print,'                    fnc=fnc, fnm=fnm, $'
	print,'	                   output_dir=output_dir,$'
	print,'                    write_png=write_png, $'
    print,'                    /UNPROCESSED, /NO_SHOW, $'
    print,'                    /WRITE_CSV,/CLEAN_DATA,$'
	print,'                    /DOWNLOAD_DATA,/SILENT,/HELP'
	return
endif

WRITE_CSV = keyword_set(WRITE_CSV)
CLEAN = keyword_set(CLEAN_DATA)
SILENT = keyword_set(SILENT)
DOWNLOAD = keyword_set(DOWNLOAD_DATA)
UNPROCESSED = keyword_set(UNPROCESSED)
NO_SHOW = keyword_set(NO_SHOW)

iflag = keyword_set(fnc) 
mflag = keyword_set(fnm)

cd,current=current_dir
if (~keyword_set(data_dir)) then dat_dir = current_dir else dat_dir = strtrim(data_dir[0],2)
if (~keyword_set(output_dir)) then begin
	out_dir = 'products'
	if (~file_test(out_dir,/DIR)) then spawn,'mkdir '+out_dir
endif else out_dir = strtrim(output_dir[0],2)

if not (iflag) or not (mflag) then begin
    if (not keyword_set(endtime)) then $
	   endtime = (strsplit(anytim(!stime, /ccsds),'T.',/EXTRACT))[0]
        jend = anytim2jd(endtime) & jend = jend.int + jend.frac
        nend = mdi_date2id(endtime)
    if (not keyword_set(starttime)) then begin
        jstart = jend - 90.0d    
        starttime = (strsplit(jd2str(jstart),'T',/EXTRACT))[0]
        nstart = mdi_date2id(starttime)
    endif

    dates = countday(starttime,endtime,nday=nday)

    if not (iflag) then fnc = ''
    if not (mflag) then fnm = ''

    for i=0l,nday-1 do begin
     
        id = mdi_date2id(dates[i])
        print,strtrim(nday-i,2)+': Loading data files for the day '+dates[i]+' ('+strtrim(id,2)+')...' 
    
        if not (iflag) then begin
            fnc_patt = 'fd_Ic_6h_01d.'+string(id,format='(i4.4)')+'.????.fits'
            fnc_i = file_search(data_dir,fnc_patt)
            if (fnc_i[0] eq '') and (DOWNLOAD) then mdi_wget,'intensity',index=id,locfname=fnc_i,output_dir=dat_dir,SILENT=SILENT
            if (fnc_i[0] eq '') then begin
                message,/INFO,'No Ic fits file found for the current day!'
                continue
            endif
            fnc = [fnc,fnc_i]
        endif
    
        if not (mflag) then begin
            fnm_patt = 'fd_M_96m_01d.'+string(id,format='(i4.4)')+'.????.fits'
            fnm_i = file_search(data_dir,fnm_patt)
            if (fnm_i[0] eq '') and (DOWNLOAD) then mdi_wget,'los_magnetic_field',index=id,locfname=fnm_i,output_dir=dat_dir,SILENT=SILENT
            if (fnm_i[0] eq '') then begin
                message,/INFO,'No M fits file found for the current day!'
                continue
            endif
            fnm = [fnm,fnm_i]
        endif
    endfor
    if (n_elements(fnc) eq 1) then begin
        message,/CONT,'Empty Ic file list!'
        return
    endif
    if (n_elements(fnm) eq 1) then begin
        message,/CONT,'Empty M file list!'
        return
    endif
    fnc = fnc[1:*]
    fnm = fnm[1:*]
endif
nfnc = n_elements(fnc)
nfnm = n_elements(fnm)
print,strtrim(nfnc,2)+' SOHO/MDI Ic fits file(s) loaded.'
print,strtrim(nfnm,2)+' SOHO/MDI M fits file(s) loaded.'

print,'Launching MDISS...'

mdi_ss_finding_observations,0,0,fnc,fnm,$
			               output_dir=out_dir,$
				           write_png=write_png, $
                           /WRITE_CSV,/WRITE_LOG, $
			               UNPROCESSED=UNPROCESSED, $
                           NO_SHOW=NO_SHOW


print,'Program has ended correctly'
print,'Elapsed time: '+strtrim(systime(/SEC)-launch_time,2)+' sec.'
print,''
END