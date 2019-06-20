;+
; NAME:
;   hogg_name2radec
; PURPOSE:
;   Takes properly-formatted iau name in "J" format and converts it into ra,dec
;   USE AT YOUR OWN RISK; NO IMPLIED WARRANTY
; INPUTS:
;   name       - name
; OUTPUTS:
;   ra         - J2000 deg
;   dec        - J2000 deg
; OPTIONAL INPUTS:
; COMMENTS:
; BUGS:
;   USE AT YOUR OWN RISK; NO IMPLIED WARRANTY
;   Requires J2000 "J".
;   Ultra-brittle.
; REVISION HISTORY:
;   2004-10-22  started by you know who
;-
pro hogg_name2radec, name,ra,dec
ra= -1
dec= -1
splog, name
rastr= strcompress(name,/remove_all)
jindx= strpos(rastr,'J')
if (jindx EQ -1) then begin
    splog, 'ERROR: no "J" found'
    return
endif
rastr= strmid(rastr,jindx+1)
splog, rastr
signindx= stregex(rastr,'[+-]')
if (signindx EQ -1) then begin
    splog, 'ERROR: no sign found'
    return
endif
decstr= (stregex(strmid(rastr,signindx),'([+-0123456789\.]+)',/subexpr,/extract))[0]
rastr= strmid(rastr,0,signindx)
splog, rastr
splog, decstr
rahstr= strmid(rastr,0,2)
ramstr= strmid(rastr,2,2)
rasstr= strmid(rastr,4)
dedstr= strmid(decstr,0,3)
demstr= strmid(decstr,3,2)
desstr= strmid(decstr,5)
help, rahstr,ramstr,rasstr,dedstr,demstr,desstr
string2radec, rahstr,ramstr,rasstr,dedstr,demstr,desstr,ra,dec
return
end
