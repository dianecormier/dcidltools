function mmt_transform, signal, nscales
;+
; NAME:
;  mmt_transform()
;
; PURPOSE:
;  Compute the multiresolution median transform of a vector. See
;  Starck et al. (199, A&AS 134, 135) for the algorithm
;
; CATEGORY:
;  Maths
;
; CALLING SEQUENCE:
;  t = mmt_transform(signal, nscales)
;
; INPUTS:
;  signal  : dblarr - a vector
;  nscales : int    - number of scales
;
; OPTIONAL INPUTS:
;  None
;
; KEYWORD PARAMETERS:
;  None
;
; OUTPUTS:
;  returns a dblarr (nscales+1, n_elements(signal)) 
;
; OPTIONAL OUTPUTS:
;  None
;
; COMMON BLOCKS:
;  None
;
; SIDE EFFECTS:
;  None
;
; RESTRICTIONS:
;  None
;
; PROCEDURE:
;  See Starck et al. (1999)
;
; EXAMPLE:
;  v = randomn(seed, 1000) 
;  t = mmt_transform(v, 6)
;
; MODIFICATION HISTORY:
;  A long time ago, JLS
;  More recently, HA
;-

nb = n_elements(signal)

output = dblarr(nscales+1, nb)

j = 0l
l = 1l
cj = signal

while j LT nscales do begin
    win = 2*l + 1
    cjp1  = median(signal, win)
    output[j, *] = cj - cjp1
    cj = cjp1
    l = 2 * l
    j = j + 1
endwhile

output[j, *] = cj

return, output

end
