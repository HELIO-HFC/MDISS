;+
; NAME:
;		mdi_ss_hfc_struct
;
; PURPOSE:
; 		Load some meta-data information concerning
;       frc, observatory, observations, pp_output, pp_info, feature
;       into structres using HFC table definitions.
;
; CATEGORY:
;		I/O
;
; GROUP:
;		SSMDI
;
; CALLING SEQUENCE:
;		mdi_ss_hfc_struct, oby_str, frc_str, obs_str, pp_info_str, pp_out_str, feat_str
;
; INPUTS:
;       None.	
;	
; OPTIONAL INPUTS:
;		None.
;
; KEYWORD PARAMETERS:
;		None.	
;
; OUTPUTS:
;		oby_str       - Structure containing information about observatory.
;		frc_str       - Structure containing information about SDO SS code.	
;		obs_str       - Structure containing information about observations.
;		pp_info_str	  - Structure containing information about Pre-Process code.
;       pp_out_str    - Structure containing information about Pre-Process obsevations.
;       feat_str      - Empty structure containing output feature parameters.	
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
;		None
;		
; CALL:
;		hfc_observations__define
;       hfc_observatory__define
;       hfc_frc_info__define
;       hfc_pp_output__define
;       hfc_pp_info__define.pro
;       hfc_sunspots__define
;
; EXAMPLE:
;		None.		
;
; MODIFICATION HISTORY:
;		Written by:		X.Bonnin, 18-OCT-2011.
;
;-


PRO mdi_ss_hfc_struct,oby_str, frc_str, obs_str, pp_info_str, pp_out_str, feat_str


;Observatory information
oby_str = {hfc_observatory}
oby_str = replicate(oby_str,2)
oby_str.id_observatory = [1,2]
oby_str.observat=['SOHO','SOHO']
oby_str.telescop=['MDI','MDI']
oby_str.instrume=['continuum','magnetogram']
oby_str.units=['Counts','Gauss']
oby_str.wavemin=['676.780','676.780']
oby_str.wavemax=['676.780','676.780']
oby_str.waveunit=['nm','nm']
oby_str.wavename=['Ni I','Ni I']
oby_str.spectral_name=['visible','line_of_sight magnetic field']
oby_str.obs_type=['Remote-sensing','Remote-sensing']
oby_str.comment=['','']

;Observations information
obs_str = {hfc_observations}
obs_str = replicate(obs_str,2)
obs_str.id_observations = [1,2]
obs_str.observatory_id = [1,2]

;Feature Recognition Code information
references = {names:'Publication',$
			  links:'http://adsabs.harvard.edu/abs/2005SoPh..228..361Z/',$
			  types:'html'} 

frc_str = {hfc_frc_info}
frc_str.id_frc_info = 1
frc_str.feature_name='SUNSPOTS'
frc_str.code='MDISS'
frc_str.version=''
frc_str.institut='MSSL'
frc_stc.person = 'Sergei Zharkov'
frc_str.contact='sz2@mssl.ucl.ac.uk'
frc_str.enc_met='RASTERSCAN/CHAINCODE'
frc_str.reference=references.links
			
;Pre-Processed code information
pp_info_str = {hfc_pp_info}
pp_info_str.id_pp_info = 1
pp_info_str.institut='OBSPM'
pp_info_str.code='MDISS'
pp_info_str.version=''
pp_info_str.person='Xavier Bonnin'
pp_info_str.contact='xavier.bonnin@obspm.fr'
pp_info_str.reference=references.links

;Pre-Processed observations information
pp_out_str = {hfc_pp_output}
pp_out_str = replicate(pp_out_str,2)
pp_out_str.id_pp_output = [1,2]
pp_out_str.pp_info_id = [1,1]
pp_out_str.observations_id = [1,2]
pp_out_str.url=['NULL','NULL']
pp_out_str.efit=[0,0]
pp_out_str.lineclean=[0,0]

;Load feature structure
feat_str = {hfc_sunspots}
feat_str.frc_info_id = 1
feat_str.observations_id_ic = 1
feat_str.observations_id_m = 2


error = 0
END
