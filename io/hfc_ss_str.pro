; ******************************************

pro  hfc_ss_str, ss_str = ss_str, obs_str=obs_str
  
; Defines structure containing observation data

JD = {jul,int:0L, frac:0.0d0}


 ;       id_observations, observatory_id, date_obs, date_end, jdint, jdfrac, exp_time, c_rotation,
 ;       bscale, bzero, bitpix, naxis1, naxis2, r_sun, center_x, center_y, cdelt1, cdelt2, quality
 ;       filename, date_obs_string, date_end_string, comment, loc_filename, id_ascii




obs_str = {hfc_obs,    $
      id:0l, $
      observatory_id:0, $
      date_obs:'', $
      date_end:'', $
      jdint:0l, $
      jdfrac:0.d0, $
      exp_time:0.d0, $
      c_rotation:0, $
      bscale:0.d0, $
      bzero:0.d0, $
      bitpix:0, $
      nx:0, ny:0, $
      r:0.d0, $
      xc:0.d0, yc:0.d0, cdx:0.d0, cdy:0.0, $
      quality:'', filename:'', date_obs_string:'', $
      date_end_string:'', comment:'', fn:'', id_ascii:0l$
   }


ss_str = {hfc_sunspot, $
    id_spot:0l, frc_info_id:0, procd_obs_id:0, run_date:'', date_obs:'', $
    feat_x_arcsec:0.d0, feat_y_arcsec:0.d0, feat_x_pix:0.d0, feat_y_pix:0.d0, $
    feat_hg_long_deg:0.d0, feat_hg_lat_deg:0.d0, feat_carr_long_deg:0.d0, feat_carr_lat_deg:0.d0,$
    feat_area_pix:0l,feat_area_deg2:0.d0, feat_diameter:0.d0, feat_mean2qsun:0.d0, $
    feat_min_int:0.d0, feat_max_int:0.d0, $
    flux_total:0.d0, flux_abs_total:0.d0, flux_min:0.d0, flux_max:0.d0, $
    umbra_number:0, umbra_area_pix:0l, $
    umbra_mean2qsun_int:0.d0, umbra_min_int:0.d0, umbra_max_int:0.d0, umbra_tot_flux:0.d0, $
    umbra_tot_absflux:0.d0, umbra_min_flux:0.d0, umbra_max_flux:0.d0, $
    brx0:0, bry0:0, brx1:0, bry1:0, abx0:0.d0, aby0:0.d0, abx1:0.d0, aby1:0.d0, $
    fn:'', rs_length:0l, rscan:''}
return
end
