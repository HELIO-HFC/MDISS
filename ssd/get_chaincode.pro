pro get_chaincode, image, code_locations=code_locations, imc=imc

; PROCEDURE
;    deducing chaincode from binary image
;
;   INPUT
;         image    2D binary image with region mask


locs=where(image ne 0)
in=size(image)

nx=in(1) & ny=in(2)
xp=locs mod nx
yp=locs / nx

x0=min(xp) & x1=max(xp)
y0=min(yp) & y1=max(yp)

;if x0 eq 0 then x0=1
;if y0 eq 0 then y0=1
imc=image
imc(*)=0
for i=x0, x1 do for j=y0, y1 do begin
    if not(image(i, j)) then continue
    da=image(i-1:i+1, j) & da=da(*)
    da=[da, image(i, j+1), image(i, j-1)]
    if (min(da) ne max(da))  then $
      imc(i, j)=1
endfor

code_locations=where(imc ne 0)
;stop
end
