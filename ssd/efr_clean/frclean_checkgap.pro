;=======================================================================

function frclean_checkgap, binary, xc, yc

;	PROCEDURE
;
;		checking the angular gap between the non zero points on
;		a binary image
;
;	INPUT
;		binary	input binary image
;		xc, yc	center in pixel coordinates
;
;	RETURN
;		gap

info=size(binary)
nx=info[1]
locs=where(binary)
n=n_elements(locs)

theta_array=dblarr(n)
dist_array=dblarr(n)
df_theta=dblarr(n)

;	Generate Angle Array

for k=0, n-1, 1 do begin
	j=fix(locs[k]/nx)
	i=locs[k] mod nx
	polar=cv_coord(FROM_RECT=[i-xc, j-yc], /TO_POLAR, /DEGREES, /DOUBLE)
	dist_array[k]=polar[1]
	if polar[0] ge 0 then theta_array[k]=polar[0] else theta_array[k]=360+polar[0]
endfor
ind=sort(theta_array)


;	Calculate Difference Array

df_theta[0]=theta_array[ind[0]]
for i=1, n-1, 1 do df_theta[i]=theta_array[ind[i]]-theta_array[ind[i-1]]

gap=max(df_theta)

return, gap
end
