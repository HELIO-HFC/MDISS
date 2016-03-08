pro ss_xml, sp, fn


n=sp->n_spots()

if n eq 0 then return
fn=strarr(n)

data=sp->dataarray()
for i=0, n-1 do begin
s1=sp->getspot(i)
s1->getbndrct, arc_arr=ab
rscan=s1->getrasterscan()
cc=s1->getchaincode()
; Create IDL structure for a Sunspot event
event = struct4event('ss')
; Populate the structure with required values
event.required.OBS_Observatory = 'SOHO'
event.required.OBS_Instrument = 'MDI'
event.required.OBS_ChannelID = 'SOHO/MDI'
event.required.OBS_MeanWavel = '6767.8'
event.required.OBS_WavelUnit = 'Angstroms'
event.required.FRM_Name = 'EGSO_SFC'
event.required.FRM_Identifier = 'EGSO_SFC'
event.required.FRM_Institute ='LMSAL'
event.required.FRM_HumanFlag = 'no'
event.required.FRM_ParamSet = 'n/a'
event.required.FRM_DateRun = anytim(!stime, /ccsds); '2007/01/03 12:00:00'
event.required.FRM_Contact = 's.zharkov at sheffield dot ac dot uk'
event.required.FRM_URL = 'n/a'
event.required.Event_StartTime = sp->getdate()
event.required.Event_EndTime = sp->getdate()
event.required.Event_CoordSys = 'UTC-HPC-TOPO'
event.required.Event_CoordUnit = 'arcsec'
event.required.Event_Coord1 = string(data(i, 0), form='(F9.4)')
event.required.Event_Coord2 =  string(data(i, 1), form='(F9.4)')
event.required.Event_C1Error = '0'
event.required.Event_C2Error = '0'
event.required.BoundBox_C1LL =  string(min(cc(0, *)), form='(F9.4)')
;string(ab(0), form='(F9.4)') ;Coordinates of lower-left
event.required.BoundBox_C2LL =  string(min(cc(1, *)), form='(F9.4)') 
;string(ab(1), form='(F9.4)') ;Corner of bounding box
event.required.BoundBox_C1UR =   string(max(cc(0, *)), form='(F9.4)') 
;string(ab(2), form='(F9.4)') ;Coordinates of upper-right    
event.required.BoundBox_C2UR =   string(max(cc(1, *)), form='(F9.4)') 
;string(ab(3), form='(F9.4)') ;Corner of bounding box 


event.optional.rasterscan=rscan
event.optional.rasterscantype='EGSO_SFC'
ccstring= string(cc(*), form='(F8.3, ",")')
k=n_elements(ccstring)
ccstr2=''
for j=0, k-1 do ccstr2=ccstr2+ccstring(j)+' '
k=strlen(ccstr2)
ccstring2=strmid(ccstr2, 0, k-2)
event.optional.Bound_ChainCode=ccstring2
event.optional.Bound_CCStartC1=''
event.optional.Bound_CCStartC2=''

date=sp->getdate()
tt=strpos(date, ':')
while tt ne -1 do begin
    date=strmid(date, 0, tt)+'.'+strmid(date, tt+1, strlen(date))
    tt=strpos(date, ':')
    endwhile
tt=strpos(date, '-')
while tt ne -1 do begin
    date=strmid(date, 0, tt)+'.'+strmid(date, tt+1, strlen(date))
    tt=strpos(date, '-')
endwhile
;stop
;If you want, add a description
event.description="My first sunspot"

;If you want, add up to ten references
;You must provide a name, link and type for each reference
;Must choose between "html", "image" and "movie" for reference_types.

event.reference_names[0] = "Publication" 
event.reference_links[0] = "http://adsabs.harvard.edu/abs/2005SoPh..228..361Z/"
event.reference_types[0] = "html"


;stop
;Now export the IDL structure to an XML file.

if i ne 0 then $
  form_string='(I'+string(fix(alog10(i))+1, form='(I1)')+')' $
  else form_string='(I1)'
i0=string(i, form=form_string) 
export_event, event, /write, suff=date+'_'+i0, filenameout=fn0


fn(i)=fn0
endfor
;, outfil="sunspot_example.xml"
end
