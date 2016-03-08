FUNCTION Sunspot::createRecord, wl

wlst=wl->getstructure()
spst=self->getstructure()
if comparefdstr(wlst, spst) eq 2 then begin

		wl->getheader, header
		wl->getfilename, filename=filename, path=path, ext=ext

		ln=strlen(path)
		k=strpos(path, 'fd_Ic')
		folder=strmid(path, k-16, ln-k+16)

		record=obj_new('SunspotDB')
		recST=record->getStruct()

		recST.instrument=fxpar(header, 'INSTRUME')
		recST.telescop=fxpar(header, 'TELESCOP')

		recST.StartTime=fxpar(header, 'T_START')
		recST.StopTime=fxpar(header, 'T_STOP')

		recST.DataBitPix=fxpar(header, 'BITPIX')

		recST.WaveLen='NOT SET'
		recST.orgnObsFilename=folder+filename+ext
		recST.Source=fxpar(header, 'ORIGIN')
		recST.FeatureName='Sunspot'
		recST.dateObs=self.date
		recST.ObsQuietSun=self.QuietSunInt
		recST.SampleCount=n_elements(*self.locs)
		recST.n_umbra=self.n_umbras

		nU=n_elements(*self.umbra)
		if nU ne 1 then	recST.UmbralArea=nU else $
			if *self.umbra eq -1 then recST.UmbralArea=0 $
				else recST.UmbralArea=nU
		recST.irradiance=self.meanInt/double(self.QuietSunInt)
		recST.helArea=self.helArea
		recST.barycArc=[self.gcArcx, self.gcArcy]

		recST.OBSxc=self.xc
		recST.OBSyc=self.yc
		recST.OBSsolarR=self.R
		recST.OBScdelt1=self.cdx
		recST.OBScdelt2=self.cdy
		recST.OBSnaxis1=self.nx
		recST.OBSnaxis2=self.ny



		clon=tim2carr(self.date)
		case self.origin of
			'soho':		GChel=arcmin2hel(self.GCarcx/60, self.GCarcy/60, date=self.date, /soho)
			'ground':	GChel=arcmin2hel(self.GCarcx/60, self.GCarcy/60, date=self.date)
		endcase
		GChel[1]=GChel[1]+clon[0]
		recST.barychel=GCHel

		self->getbndrct, x0, y0, x1, y1

		recST.boundRpix[0, 0]=x0
		recST.boundRpix[0, 1]=y0
		recST.boundRpix[1, 0]=x1
		recST.boundRpix[1, 1]=y1

		recST.boundRarc[*, 0]=(recST.boundRpix[*, 0]-self.xc)*self.cdx
		recST.boundRarc[*, 1]=(recST.boundRpix[*, 1]-self.yc)*self.cdy
		recST.CRnum=tim2carr(self.date, /dc)

		if abs(self.gcArcx) gt 800. then res=0.01 else res=0.1
		helArray=self->pix2hel(res=res)
		recST.diameter=calcdiam(imb=helArray)*res
return, recST
endif else begin
		print, 'Incorrect Observation'
		return, -1
endelse

END

