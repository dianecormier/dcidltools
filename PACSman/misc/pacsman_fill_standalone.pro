
;--------------------------------------------------------------------------
;fill zeros within maps
function PACSman_fill, imag, h

  nx = n_elements(imag[*, 0])
  ny = n_elements(imag[0, *])
  ;print, nx, ny

  ind_zeros = where( finite(imag, /nan), count )
  
  if count eq 0. then return, imag
  
  for boxsize = 35, 7, -2 do begin
  ;boxsize = 25                    ;we check n x n patches only, n needs to be odd
  midboxsize = (boxsize-1) / 2
  for x = midboxsize, nx - midboxsize - 1 do for y = midboxsize, ny - midboxsize - 1 do begin
    box = imag[x-midboxsize:x+midboxsize, y-midboxsize:y+midboxsize]
    ind = where( finite(box, /nan), c )
    if c eq 0 or c eq boxsize^2 then continue ;nothing to fill or all nans
    ;print, x, y, c
    out = [ reform(box[0, 0:boxsize-1]), reform(box[boxsize-1, 0:boxsize-1]), reform(box[0:boxsize-1, 0]), reform(box[0:boxsize-1, boxsize-1])]
    ind = where( finite(out, /nan), c )
    if c gt 0 then continue ;at least one nan in the outer edge of the box
    print, boxsize, x, y
    newbox1 = box
    for x2 = 0, n_elements(box[*, 0]) - 1 do begin
       ind = where( finite(box[x2, *], /nan), c, complement = noind )
       if c gt 0 and n_elements(noind) gt 3 then newbox1[x2, ind] = interpol(box[x2, noind], noind, ind)
    endfor
    newbox2 = box 
    for y2 = 0, n_elements(box[0, *]) - 1 do begin
       ind = where( finite(box[*, y2], /nan), c, complement = noind )
       if c gt 0 and n_elements(noind) gt 3 then newbox2[ind, y2] = interpol(box[noind, y2], noind, ind)
    endfor
    newbox = 0.5 * ( newbox1 + newbox2 )
    imag[x-midboxsize:x+midboxsize, y-midboxsize:y+midboxsize] = newbox
    ;horizontal = imag(*, y)
    ;ind = where( finite(horizontal) )
    ;val1 = interpol(horizontal(ind), ind, x)
    ;vertical = imag(x, *)
    ;ind = where( finite(vertical) )
    ;val2 = interpol(vertical(ind), ind, y)
    ;imag(x, y) = avg([val1, val2])
  endfor
  endfor
  
  
  return, imag
  
end
