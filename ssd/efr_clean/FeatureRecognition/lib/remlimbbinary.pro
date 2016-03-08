;-------------------------------------------------------

function RemLimbBinary, binary, xc, yc, R, nx

; Function to remove limb edge from binary image

binary2=binary
locs=where(binary ne 0)
N=n_elements(locs)
ypos=locs / nx
xpos=locs mod nx
R2=(R-5)*(R-5)
for i=0l, N-1 do $
	if (xpos[i]-xc)*(xpos[i]-xc)+(ypos[i]-yc)*(ypos[i]-yc) gt R2 then binary2[locs[i]]=0

return, binary2
end
