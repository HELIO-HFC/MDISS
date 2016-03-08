PRO mdi_wget,physobs,date_obs=date_obs,index=index, $
             output_dir=output_dir,locfname=locfname, $
             tries=tries,provider=provider,SILENT=SILENT

;+
; NAME:
;		mdi_wget
;
; PURPOSE:
;		Download SOHO/MDI fits file from distant server using wget command.
;
; CATEGORY:
;		I/O
;
; GROUP:
;		None.
;
; CALLING SEQUENCE:
;		IDL> mdi_wget,physobs,date_obs=date_obs
;
; INPUTS:
;		physobs - Physical parameter for which fits file(s) must downloaded.
;                 It can be 'intensity' (fd_Ic_6h_01d files) or 'los_magnetic_field' (fd_M_96m_01d files).
;
; OPTIONAL INPUTS:
;       date_obs   - Scalar or Vector containing the date(s) of observation for which fits file(s) must
;                    be downloaded.
;       index      - Scalar or Vector containing the date index(es) of fits file(s) to download.
;                    (Ignored if date_obs is set.)
;       output_dir - Full path name to the directory where file(s) will be saved.
;       provider   - Name of the data provider. Only medoc is available.
;       tries      - Number of server connection try.
;                    Default is 1.
;
; KEYWORD PARAMETERS: 
;		/SILENT		    - Quiet mode.
;		
; OUTPUTS:
;		None.				
;
; OPTIONAL OUTPUTS:
;		locfname - Full path name of the downloaded file(s).		
;		
; COMMON BLOCKS:
;		None.		
;
; SIDE EFFECTS:
;		None.
;
; RESTRICTIONS/COMMENTS:
;		- An internet access is required to download SOHO/MDI data from MEDOC using wget command.
;
; CALL:
;		mdi_date2id
;       mdi_id2date
;
; EXAMPLE:
;		None.
;		
; MODIFICATION HISTORY:
;		Written by:		X.Bonnin, 10-JAN-2012. 
;
;-

if (n_params() lt 1) then begin
    message,'Call is:'
    print,'mdi_wget,physobs,date_obs=date_obs, $'
    print,'         index=index, locfname=locfname, $'
    print,'         output_dir=output_dir,tries=tries, $'
    print,'         provider=provider,/SILENT'
    return
endif

SILENT = keyword_set(SILENT)

if (not keyword_set(tries)) then try = '1' else try = strtrim(tries[0],2)
provider = 'medoc' ;only provider available for now
pvdr = strlowcase(strtrim(provider[0],2))

pobs = strlowcase(strtrim(physobs[0],2))

dflag = +keyword_set(date_obs)
iflag = -keyword_set(index)

if not (keyword_set(output_dir)) then cd,current=output_dir

if not (dflag) and not (iflag) then begin
    message,/CONT,'You must provide the date(s) of observation or the index(es) of file(s)!'
    return
endif

if (dflag) then nobs = n_elements(date_obs) else nobs = n_elements(index)
if not (SILENT) then print,strtrim(nobs,2)+' day(s) to download.'

case pvdr of
    'medoc':url = 'www.medoc-ias.u-psud.fr/archive/5/private/data/processed/mdi/'
    else:begin
        message,/CONT,'Unknown provider!'
        return
    end
endcase

locfname = ""
for i=0l,nobs-1l do begin
   
    if (dflag) then begin
        date_obs_i = date_obs[i]
        index_i = mdi_date2id(date_obs_i)  
    endif else begin
        index_i = index[i]
        date_obs_i = mdi_id2date(index_i)
    endelse

    print,'Downloading SOHO/MDI - '+pobs+' fits file(s) for the day '+date_obs_i+' ('+strtrim(index_i,2)+')...'

    if (pobs eq 'intensity') then begin
        dirname_i4 = 'fd_Ic_6h_01d/fd_Ic_6h_01d.'+string(index_i,format='(i4.4)')+'/' 
        dirname_i6 = 'fd_Ic_6h_01d/fd_Ic_6h_01d.'+string(index_i,format='(i6.6)')+'/'
        basename_i = 'fd_Ic_6h_01d.'+string(index_i,format='(i4.4)')
    endif else begin
        dirname_i4  = 'fd_M_6h_01d/fd_M_6h_01d.'+string(index_i,format='(i4.4)')+'/' 
        dirname_i6  = 'fd_M_6h_01d/fd_M_6h_01d.'+string(index_i,format='(i6.6)')+'/'
        basename_i = 'fd_M_6h_01d.'+string(index_i,format='(i4.4)') 
    endelse

    if (SILENT) then vb = '-q' else vb = '-v'
    url_i4 = url + dirname_i4
    url_i6 = url + dirname_i6
    cmd4 = 'wget '+vb+' -r -nc -nd -t '+try+' --no-parent -A.fits '+url_i4+' -P '+output_dir
    
    if not (SILENT) then print,cmd4
    spawn,cmd4
    
    locfname_i = file_search(output_dir + path_sep() + basename_i+'.????.fits')
    if (locfname_i[0] ne '') then begin
        locfname = [locfname,locfname_i] 
        continue
    endif

    cmd6 = 'wget '+vb+' -r -nc -nd -t '+try+' --no-parent -A.fits '+url_i6+' -P '+output_dir
    if not (SILENT) then print,cmd6
    spawn,cmd6
        
    locfname_i = file_search(output_dir + path_sep() + basename_i+'.????.fits')
    if (locfname_i[0] ne '') then begin
        locfname = [locfname,locfname_i] 
        continue
    endif else begin
        if not (SILENT) then print,'No file found!'
    endelse

endfor
nloc = n_elements(locfname)
if (nloc gt 1) then locfname = locfname[1:*]
if not (SILENT) then print,strtrim(nloc-1l,2)+ ' file(s) downloaded.'

END