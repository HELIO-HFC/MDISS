;=======================================================

PRO Feature::GetBndRct, x0, y0, x1, y1, arc_arr=arc_arr

locs=*self.locs
xp=locs mod self.nx
yp=locs / self.nx

x0=min(xp)
y0=min(yp)

x1=max(xp)
y1=max(yp)

arc_arr=dblarr(4)
		br_arcx0=(x0-self.xc)*self.cdx
		br_arcx1=(x1-self.xc)*self.cdx
		br_arcy0=(y0-self.yc)*self.cdy
		br_arcy1=(y1-self.yc)*self.cdy

	arc_arr[0]=br_arcx0
	arc_arr[1]=br_arcy0
	arc_arr[2]=br_arcx1
	arc_arr[3]=br_arcy1

END
