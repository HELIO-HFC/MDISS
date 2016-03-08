pro ss_read_output_ascii, fno, fnss, obs, ss

hfc_ss_str, ss_str = ss_str, obs_str=ssobs

obs_string=rd_tfile(fno(0), delim=';', 25)
in=size(obs_string)
nn=in(2)

obs=replicate(ssobs, nn)
obs.id=reform(long(obs_string(0, *)))

       obs.observatory_id=reform(fix(obs_string(1, *)))
       obs.date_obs=reform((obs_string(2, *)))
       obs.date_end=reform((obs_string(3, *)))
       obs.jdint=reform(long(obs_string(4, *)))
      obs.jdfrac=reform(double(obs_string(5, *)))
      obs.exp_time=reform(double(obs_string(6, *)))
      obs.c_rotation=reform(fix(obs_string(7, *)))
      obs.bscale=reform(double(obs_string(8, *)))
      obs.bzero=reform(double(obs_string(9, *)))
      obs.bitpix=reform(fix(obs_string(10, *))) & j=11
      obs.nx=reform(fix(obs_string(j++, *)))
      obs.ny=reform(fix(obs_string(j++, *)))
      obs.r=reform(double(obs_string(j++, *)))
      obs.xc=reform(double(obs_string(j++, *))) 
      obs.yc=reform(double(obs_string(j++, *)))
      obs.cdx=reform(double(obs_string(j++, *)))
      obs.cdy=reform(double(obs_string(j++, *)))
      obs.quality=reform((obs_string(j++, *))) 
      obs.filename=reform((obs_string(j++, *))) 
      obs.date_obs_string=reform((obs_string(j++, *)))
      obs.date_end_string=reform((obs_string(j++, *))) 
      obs.comment = reform((obs_string(j++, *))) 
      obs.fn=reform((obs_string(j++, *)))
      ;obs.id_ascii=reform(long(obs_string(j++, *)))
                       
    
    ss_string=rd_tfile(fnss(0), delim=';', 43)
    in=size(ss_string)
    nn=in(2)

    ss=replicate(ss_str, nn)
    
    j=0
    ss.id_spot =reform(long(ss_string(j++, *)))  
    ss.frc_info_id =reform(fix(ss_string(j++, *))) 
    ss.procd_obs_id =reform(fix(ss_string(j++, *)))  
    ss.run_date =reform((ss_string(j++, *)))   
    ss.date_obs =reform((ss_string(j++, *)))   
    ss.feat_x_arcsec =reform(double(ss_string(j++, *)))  
    ss.feat_y_arcsec =reform(double(ss_string(j++, *)))  
    ss.feat_x_pix =reform(double(ss_string(j++, *)))  
    ss.feat_y_pix =reform(double(ss_string(j++, *)))  
    ss.feat_hg_long_deg =reform(double(ss_string(j++, *))) 
    
    ss.feat_hg_lat_deg =reform(double(ss_string(j++, *))) 
    ss.feat_carr_long_deg =reform(double(ss_string(j++, *)))
    ss.feat_carr_lat_deg =reform(double(ss_string(j++, *))) 
    ss.feat_area_pix =reform(long(ss_string(j++, *)))
    ss.feat_area_deg2 =reform(double(ss_string(j++, *)))  
    ss.feat_diameter =reform(double(ss_string(j++, *))) 
    ss.feat_mean2qsun =reform(double(ss_string(j++, *)))  
    ss.feat_min_int =reform(double(ss_string(j++, *))) 
    ss.feat_max_int =reform(double(ss_string(j++, *))) 
    ss.flux_total =reform(double(ss_string(j++, *)))
    ss.flux_abs_total =reform(double(ss_string(j++, *)))
    ss.flux_min =reform(double(ss_string(j++, *))) 
    ss.flux_max =reform(double(ss_string(j++, *))) 
    ss.umbra_number =reform(fix(ss_string(j++, *)))
    ss.umbra_area_pix =reform(long(ss_string(j++, *)))  
    ss.umbra_mean2qsun_int =reform(double(ss_string(j++, *)))  
    ss.umbra_min_int =reform(double(ss_string(j++, *))) 
    ss.umbra_max_int =reform(double(ss_string(j++, *)))  
    ss.umbra_tot_flux =reform(double(ss_string(j++, *))) 
    ss.umbra_tot_absflux =reform(double(ss_string(j++, *)))  
    ss.umbra_min_flux =reform(double(ss_string(j++, *)))  
    ss.umbra_max_flux =reform(double(ss_string(j++, *)))  
    ss.brx0 =reform(fix(ss_string(j++, *))) 
    ss.bry0 =reform(fix(ss_string(j++, *))) 
    ss.brx1 =reform(fix(ss_string(j++, *))) 
    ss.bry1 =reform(fix(ss_string(j++, *))) 
    ss.abx0 =reform(double(ss_string(j++, *)))  
    ss.aby0 =reform(double(ss_string(j++, *))) 
    ss.abx1 =reform(double(ss_string(j++, *)))  
    ss.aby1 =reform(double(ss_string(j++, *))) 
    ss.fn= reform((ss_string(j++, *)))
    ss.rs_length=reform(long(ss_string(j++, *)))  
    ss.rscan=reform((ss_string(j++, *)))
    
    jd=anytim2jd(ss.date_obs) & jd=jd.int+jd.frac
    plot, jd, ss.feat_area_deg2, xtickf='label_date', /yst, /xst
    stop
 ;  endwhile             

end
