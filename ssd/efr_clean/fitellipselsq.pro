;=============================================================================

function fitellipselsq, binary, xc, yc

;	Fitting Ellipse
;	RESULT:	vector of 6 coefficients, double
;	PARAMETERS:
;		binary:		binary image containing data
;					to be fitted
;		xc,yc:		coordinate center wrt to image coordinates


locs=where(binary ne 0)
in=size(locs)
N=in[1]
;	print, 'N =', N
d=dblarr(6,N)
S_inv=dblarr(6,6)
Q=dblarr(6,6)
C=dblarr(6,6)
C[2,0]=2
C[1,1]=-1
C[0,2]=2
;xc=state.x_center
;yc=state.y_center
im=binary
info=size(im)
nx=info[1]
ny=info[2]
for i=0l,N-1,1 do $
	begin
		yi=fix(locs[i]/nx)
		xi=locs[i]-long(yi)*nx
		d[0,i]=(xi-xc)*(xi-xc)
		d[1,i]=(xi-xc)*(yi-yc)
		d[2,i]=(yi-yc)*(yi-yc)
		d[3,i]=xi-xc
		d[4,i]=yi-yc
		d[5,i]=1
	endfor

Dr=dblarr(6,N)
Dr[*,0:N-1]=d[*,0:N-1]


Dr_T=transpose(Dr)
S=dblarr(6, 6)
S=Matrix_Multiply(Dr, Dr_T)
;	print, 'Scatter Matrix:'
;	print, S
A=S
Trired, A, D1, E,/double
;	print, 'Diagonal:'
;	print, D
Triql, D1,E,A,/double

;vec=A[*,0]
;vec1=vec
;for i=0,5,1 do vec1[i]=D[0]*vec[i]
;print, S # vec - vec1

;	print, 'Scatter eigenvalues:',D
if min(abs(D1)) lt 1d-06 then $
		for i=0,5,1 do begin
			if abs(D1[i]) lt 1e-06 then begin
				vectorx=A[*,i]
				;print, i
				;print, vectorx
				if (vectorx[1]^2-4.0d00*vectorx[0]*vectorx[2]) le 0 then begin
					vector=vectorx
					;print, vector
					endif else print, (vectorx[1]^2-4.0d00*vectorx[0]*vectorx[2])
			endif
			endfor else begin
		for i=0,5,1 do begin
			if D1[i] gt 0 then S_inv[i,i]=1./D1[i] else S_inv[i,i]=0
			S_inv1=Matrix_Multiply(A,Matrix_Multiply(S_inv,A,/btranspose))
			;	print, S_inv1#S
			B=S_inv1
			Trired, B, D2, E,/double
			Triql, D2, E, B,/double

			;	print, 'X eigenvalues:',D1
			j=where(D2 eq max(D2))
			;eigen=1./max(D1)
			;for i=0,5,1 do if abs(D1[i]) gt 1d-10 then begin
			;	print, i
			;	print, D[j]
			eigen=1.d00/D1[j]
			;	print, 'Egenvalue:',eigen
			e_vec=B[*,j]
			vector=double(e_vec)
			;end
		endfor
		end
for i=0,5,1 do vector[i]=vector[i]/vector[5]

;print, 'Coefficients:', vector

Return, vector
end
