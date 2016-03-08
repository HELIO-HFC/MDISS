pro check_detected_sunspots, ss, obs, cwait=cwait, istart=istart

n=n_elements(obs)

if not keyword_set(istart) then istart=0l
for i=istart, n-1 do begin
  locs=where(ss.date_obs eq obs(i).date_obs)
  if locs(0) ne -1 then begin
   ; stop
    mreadfits, ss(locs(0)).fn, ind, da
    da=rot(da, ind.solar_p, 1, ind.x0, ind.y0)
    imb=da & imb(*)=0 & imb0=imb
    da1=da
    
    for j=0, n_elements(locs)-1 do begin
      imb0(*)=0
      s0=ss(locs(j))
      nx=s0.brx1-s0.brx0+1
      xx=efr_raster2image(s0.rscan, nx)
      imb0(s0.brx0:s0.brx1, s0.bry0:s0.bry1)=XX
      imb=imb+imb0;+3*j
     end
     
     xlocs=where(imb mod 3 eq 1) 
     ylocs=where(imb mod 3 eq 2)
     
     if xlocs(0) ne -1 then da1(xlocs)=0
     if ylocs(0) ne -1 then da1(ylocs)=max(da1)
     
    
     x0=min(ss(locs).brx0)-10
     x1=max(ss(locs).brx1)+10
     
     y0=min(ss(locs).bry0)-10
     y1=max(ss(locs).bry1)+10

	!p.multi=[0, 1, 2]
	tvframe, da(x0:x1, y0:y1)>0, /asp, tit=ss(locs(0)).date_obs
     tvframe, da1(x0:x1, y0:y1)>0, /asp, tit=ss(locs(0)).fn+string(i)
	!p.multi=0
     ;contour, imb(x0:x1, y0:y1) mod 3, /ov, lev=[1, 2]
    if not keyword_set(cwait) then stop else wait, cwait
  endif
endfor


end
