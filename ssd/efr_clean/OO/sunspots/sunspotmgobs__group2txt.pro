PRO SunspotMgObs::group2txt, gr_ind, ng=ng

spa=*self.spotarr
image=*self.image
mgim=*self.mgdata
eol='999.9999'
n=n_elements(spa)


fn='c:/gl2003-10-28_062333_ascii.txt'
get_lun, lun
openw, lun, fn
print, n

gr_ind=intarr(n)
group_number=0
for i=0, n-1 do begin

		sp=spa[i]
		locs=sp->getlocs()
		ulocs=sp->getumbra()
		im=image
		tmp=bytarr(self.nx, self.ny)
		im[locs]=2*im[locs]
		tmp[locs]=1
		if ulocs[0] ne -1 then begin
				tmp[locs[ulocs]]=2
			endif

if not keyword_set(ng) then begin
		window, 1, xs=self.nx, ys=self.ny, xp=300
		tvscl, tmp

		window, 2, xs=self.nx, ys=self.ny, xp=600
		tvscl, im

	read, group_number
  endif else	group_number=0



		gr_ind[i]=group_number
	;	stop




		; get gravity center coordinates (helio & arcseconds)
		sp->getgc, gc_arcx, gc_arcy, gc_carlon, gc_carlat, diameter=diameter
	; get #pixels in the sunspot, umbra, #of umbras & helio size
		samplecount=sp->getscounts(n_umb=n_umb, ucount=upixsize, area=area)

	; get intensity parameters - means and quiet sun coefficients
		feat_mean_int=sp->getphysparams( image, max_value=feat_max_int, $
				min_value=feat_min_int, m2QS=mean_int_ratio, /direct_image)

	; get magnetic parameters
		locs=sp->getlocs()
		tot_flux=total(mgim[locs])
		abs_flux=total(abs(mgim[locs]))
		max_flux=max(mgim[locs])
		min_flux=min(mgim[locs])


	if n_umb ne 0 then begin
			ulocs=sp->getumbra()
			umblocs=locs[ulocs]
			utot_flux=total(mgim[umblocs])
			uabs_flux=total(abs(mgim[umblocs]))
			umax_flux=max(mgim[umblocs])
			umin_flux=min(mgim[umblocs])
		endif else begin
			utot_flux=0
			uabs_flux=0
			umax_flux=0
			umin_flux=0

		end
			fmt='(2(I6, ";"), 6(F15.8, ";"), 5(I8, ";"), 10(F20.8, ";"), A15)'

			printf, lun, $
			format=fmt,$
			i, gr_ind[i], gc_arcx, gc_arcy, gc_carlon, gc_carlat, $
				area, diameter, samplecount, n_umb, upixsize, $
				feat_max_int, feat_min_int, feat_mean_int, self.QuietSunInt, $
				max_flux, min_flux, tot_flux, abs_flux, $
				umax_flux, umin_flux, utot_flux, uabs_flux, eol
	endfor

close, lun
free_lun, lun
end
