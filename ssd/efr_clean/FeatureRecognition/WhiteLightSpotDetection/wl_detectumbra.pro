function WL_DetectUmbra, image, binary, imdi, bimdi, nx, ny, qsint


LRIM=label_region(binary)
hlrm=histogram(LRIM, reverse_indices=r)

Rnum=n_elements(hlrm)-1
;print, 'Number of Sunspots on the Image:', Rnum


imageUmbra=bytarr(nx, ny)
for i=1, Rnum do begin

	plocs=r(r[i]:r[i+1]-1)


	area=n_elements(plocs)


	xpos=plocs mod nx
	ypos=plocs / nx

	xmin=min(xpos)
	xmax=max(xpos)
	ymin=min(ypos)
	ymax=max(ypos)
	nxc=xmax-xmin+1
	nyc=ymax-ymin+1


	sc=10
	imdis1=intarr(nx, ny)
	imdis1[*, *]=-1
	imdis1[plocs]=image[plocs]
	im_dis=imdis1[xmin:xmax, ymin:ymax]

;	window, 0, xs=nxc*sc, ys=nyc*sc, xp=0, yp=0
;	tvscl, congrid(im_dis, nxc*sc, nyc*sc)
	if area gt 5 then begin
		bin_umbra=bytarr(nxc, nyc)
		mom=moment(image[plocs], mdev=mdev)
;		print, 'quiet sun int:', qsint
;		print, 'mean:', mom[0]
;		print, 'st dev:', sqrt(mom[1])
;		print, 'mean abs dev:', mdev


		if min(image[plocs]) lt .7*Qsint then begin

				Threshold=mom[0]-sqrt(mom[1]) > .55*Qsint
;				print, QSint, .7*Qsint, Threshold
;				print, ''
				xl=where(im_dis ge 0 and im_dis lt Threshold)
		;		bin_umbra[xl]=1

		;		window, 2, xs=nxc*sc, ys=nyc*sc, xp=nxc*sc, yp=0
		;		tvscl, congrid(bin_umbra, nxc*sc, nyc*sc)
				if xl[0] ne -1 then	imageumbra[(xl mod nxc)+xmin, (xl / nxc) + ymin]=1
				endif ;else print, 'No Umbra in this spot'
	endif else begin
			xl=where(image[plocs] lt .6*qsint)
			if xl[0] ne -1 then begin
				imageumbra[plocs[xl]]=1
			;	print, 'Small Sunspot, umbra detected'
			endif
		endelse
;	stop, im_dis
end

;

return, imageumbra
end