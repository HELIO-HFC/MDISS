
;======================================================================

function Image2Raster, image, code=code

n=n_elements(image)
;print, 'Straight String Length:', n

str=''

;	straight forward raster

if not keyword_set(code) then $
	for i=0l, n-1 do $
	;	if image[i] eq 0 then str=str+'0' else $
			str=str+trim(string(image[i], /print)) $
	else begin

		sprtr='.'
		abs_count=0l
		while abs_count lt n do begin

			count=0
			flag=0
			value=image[abs_count]
			while (abs_count+count) lt n and flag eq 0 do $
				if image[abs_count+count] eq value then count=count+1 else flag=1
			abs_count=abs_count+count

			str=str+trim(string(value, /print))+$
					trim(string(count, /print))+ sprtr
		endwhile
endelse


;stop

cnt=0


return, str
end
