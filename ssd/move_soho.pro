pro move_soho, fnc, fnm

n=n_elements(fnc)

for i=0, n-1 do begin
    imc=readfits(fnc(i), hdc)
    imm=readfits(fnm(i), hdm)

    lc=strlen(fnc(i))
    lm=strlen(fnm(i))

    writefits, 'SOHO/'+strmid(fnc[i], 20, /rev), imc, hdc
    writefits, 'SOHO/'+strmid(fnm[i], 20, /rev), imm, hdm
;    stop

endfor


end
