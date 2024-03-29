;+
; NAME:
;	  logsum
; PURPOSE: (one line)
;   sum quantities when you have them as logs (preserving dynamic range)
; DESCRIPTION:
;
; CATEGORY:
;       Numerical
; CALLING SEQUENCE:
;   res= logsum(logs [,/double])
; INPUTS:
; OPTIONAL INPUT PARAMETERS:
; KEYWORD PARAMETERS:
;   /double - assume double precision input (otherwise assumes float)
; OUTPUTS:
; COMMON BLOCKS:
; SIDE EFFECTS:
; RESTRICTIONS:
;   seems to not have expected precision
; PROCEDURE:
; MODIFICATION HISTORY:
;	  Blanton and Roweis 2003-02-18 
;-
function logsum, logs, double=double, const=const

maxlog=max(logs)
logxmax=alog((machar(double=double)).xmax)
const=logxmax-(2.D)*(alog(double(n_elements(logs)))-maxlog)
logsum=alog(total(exp(logs+const),double=double))-const
return,logsum

end
