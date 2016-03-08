;===========================================================

function FitEllipseSVD, binary, xc, yc

;	FUNCTION
;			fitting ellipse to the binary image
;
;	INPUT
;		binary		binary image
;		xc, yc		coordinate centre ini pixel coordinates
;
;	OUTPUT
;		6-vector of ellipse algebraic coefficients

locs=where(binary ne 0)
N=n_elements(locs)
vector=dblarr(6)

if N lt 2 then return, vector
;	print, 'N =', N
d=dblarr(6,N)
S_inv=dblarr(6,6)
Q=dblarr(6,6)

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

SVDC, S, W, U, V, /double

;recover the least value

;print, W
least=W[0]
m=0
for i=1, 5, 1 do $
	if W[i] lt least then begin
		least=W[i]
		m=i
	endif

for i=0, 5, 1 do vector[i]=V[m, i];/V[m, 5]

Return, vector
end