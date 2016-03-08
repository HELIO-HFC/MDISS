;=======================================================================

function frclean_gcentre, im

;	FUNCTION
;		calculates the gravity centre of the im
;	RETURN:
;		gravity centre coordinates
;	Restrictions:
;		Grayscale Input Image (use of size function)

info=size(im)
nx=info[1]
ny=info[2]
cc=dblarr(2)

sum=0.D0
xmom=0.D0
ymom=0.D0

for i=0,nx-1,1 do $
	for j=0, ny-1,1 do begin
		sum=sum+im[i,j]
		xmom=xmom+double(im[i,j])*double(i)
		ymom=ymom+double(im[i,j])*double(j)
	end
cc[0]=xmom/sum
cc[1]=ymom/sum

return, cc
end

