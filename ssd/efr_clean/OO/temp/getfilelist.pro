;=======================================================

function GetFileList, ext, path

; FUNCTION
;		to return the list of files with
; 		given extension located inside the path
;		[checking level1 subdirectories]

l=strlen(path)
last_smb = strmid(path, l-1, 1)
if last_smb ne '/' and last_smb ne '\' then $
		path=path+'\'
print, path

fn=findfile(path+'*', count=k)

;if fn eq '' then stop
;print, k

fcnt=15
list=strarr(k*fcnt)
wlfn=strarr(k*fcnt)
n0=0

for i0=2, k-1 do begin

	fn2=findfile(fn[i0]+ext, count=l)
		if l ne 0 then begin
	;	help, fn2
		list[n0:n0+l-1]=fn2
	;	WLfn=strarr(l)

		n0=n0+l
	end
endfor

if n0 eq 0 then return, -1

list=temporary(list[0:n0-1])
ind=sort(list)
return, list[ind]
end

