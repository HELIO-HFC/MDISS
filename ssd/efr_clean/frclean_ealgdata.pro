;==================================================================================

pro frclean_ealgdata, xc, yc, v, xc1, yc1, R1, R2, theta, noprint=noprint

;		PROCEDURE
;			calculating ellipse from algebraic coefficients
;		INPUT:
;			xc, yc		coordinate centre pixel coordinates
;			v			6 element vector containing algebraic coefficients
;
;		OUTPUT:
;			xc1, yc1	ellipse centre pixel coordinates
;			R1, R2  	ellipse axis
;			theta		angle
;
;		KEYWORD:
;			noprint		if not set, prints the resulting parameters on the screen


v2=dblarr(6)

det=v[1]*v[1]-4*v[0]*v[2]

if det eq 0 then return

xc1=xc+(2*v[2]*v[3]-v[1]*v[4])/det
yc1=yc+(2*v[0]*v[4]-v[1]*v[3])/det


;Radius etc

root=(v[0]-v[2])*(v[0]-v[2])+v[1]*v[1]
r1=.5*(v[0]+v[2]+sqrt(root))
r2=.5*(v[0]+v[2]-sqrt(root))

R01=2.0*(r1-v[0])/v[1]
u1=1./sqrt(1.0+R01*R01)
v01=R01*u1
t1=atan(v01, u1)


R02=2.0*(r2-v[0])/v[1]
u2=1./sqrt(1.0+R02*R02)
v02=R02*u2
t2=atan(v02,u2)



a1=1/sqrt(abs(r1))
a2=1/sqrt(abs(r2))





theta=t1
R1=a1
R2=a2

if not keyword_set(noprint) then begin

print, 'Ellipse Center:',  xc1, yc1
print, 'Ellipse Radius:', R1, R2
print, 'theta:',  (t1)*180/!pi ;,(theta)*180/!pi,  ;, (theta)*180/!pi

endif

theta=t1
R1=a1
R2=a2

end

