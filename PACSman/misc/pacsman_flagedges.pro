function pacsman_flagedges, imag, h, arr = arr

  extast, h, astr
  
  flux = reform(imag(*, *, 0))
  flag = flux*0.
  ;help, flux, flag
  s = size(flux, /dim)
  ;help, s
  
  i = 0
  arr = intarr(1000, 2)
  for x = 0, s(0) - 1 do for y = 0, s(1) - 1 do begin
    
    aroundfinite = 1
    if x ge 1 then if not finite(flux(x-1, y)) then aroundfinite = 0
    if y ge 1 then if not finite(flux(x, y-1)) then aroundfinite = 0
    if x lt s(0)-1 then if not finite(flux(x+1, y)) then aroundfinite = 0
    if y lt s(1)-1 then if not finite(flux(x, y+1)) then aroundfinite = 0
    
    if finite(flux(x, y)) and ((aroundfinite eq 0) or (x eq 0) or (x eq s(0)-1) or (y eq 0) or (y eq s(1)-1)) then begin
      arr(i, *) = [x, y] 
      flag(x,y) = 1.
      i += 1
    endif
   
  endfor
  ;print, flag(*,100)
  ;print, x, y 
  npoints = i
  arr = arr(0:npoints - 1, *)
  
  ind_edges = where(flag eq 1., c)
  
  ;writefits, '~/Donnees/PACS/flag_test.fits', flag, h
  
  return, ind_edges
  
end
