
pro local2earth, xypos_pix=xypos_pix, index=index, xypos_hc_earth=xypos_hc_earth, $
  clonlat=clonlat, clon_local=clon_local, array_clonlat=array_clonlat

;+
; Purpose:
;   Calculate Earth view xcen and ycen coordinates from SECCHI xcen coordinates
; Inputs:
;   xypos_pix - Input pair in pixel coordinates (relative to lower
;               left corner)
;   index structure for associated SECCHI image
; Outputs:
;   xypos_hc_earth - Two element [x,y] vector of coordinates in earth
;                 view heliocentric coordinates 
;   clonlat - Two element vector of Carrington longitude and latitude
;             corresponding to input coordinate pair
;-

n=n_elements(xypos_pix)
array_clonlat=dblarr(2, n/2)
;if gt_tagval(index, /instrume, missing='') eq 'SECCHI' then $
;        box_message, 'Applying SECCHI->EARTH transform'


for i=0, n/2-1 do begin

    xpos_asec = comp_fits_crval(index.xcen,index.cdelt1,index.naxis1,xypos_pix(0, i))
    ypos_asec = comp_fits_crval(index.ycen,index.cdelt2,index.naxis2,xypos_pix(1, i))

    if gt_tagval(index, /instrume, missing='') eq 'SECCHI' then begin
     ;   box_message, 'Applying SECCHI->EARTH transform'
        if not tag_exist(index,'CROTA') then index = secchi_index2crota(index)
        index_last = last_nelem(index)
        index_last.crpix1 = 0
        index_last.crpix2 = 0
        index_last.xcen = xpos_asec
        index_last.ycen = ypos_asec
    endif

; Create normalized (rsun=1) Cartesian corrdinates:
    yn = index_last.xcen/index_last.rsun
    zn = index_last.ycen/index_last.rsun
    temp = yn*yn + zn*zn
    xn = sqrt(1 - temp)

; Convert to local view spherical coordinates using local b0, p angles:
    rtp0 = c2s([xn,yn,zn], b0=-index_last.hglt_obs, roll=-index_last.crota)

; Local view lat and lon:
    lat = rtp0[1]
    lon = rtp0[2]

; Earth-view lat and lon:
    lat_earth = lat
    lon_earth = lon + index_last.hgln_obs

    clon_local = index_last.crln_obs
    clon_earth = (tim2carr(anytim(index_last.date_obs,/yoh)) + lon_earth) mod 360
    if clon_earth gt 180. then clon_earth = clon_earth - 360.
    clonlat = [clon_earth, lat_earth]
;    clonlat = [clon_local, lat_earth]

; Earth-view heliocentric:
    xypos_hc_earth = hel2arcmin(lat_earth, lon_earth, date=anytim(index_last.date_obs,/yoh)) * 60.
    array_clonlat(*, i)=clonlat

endfor

end
