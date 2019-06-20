pro pacsman_makereg, file, plot = plot, color = color

  if not keyword_set(color) then color = 'black'
  
  im = readfits(file, h, /sil)
  extast, h, astr
  
  flux = reform(im(*, *, 0))
  s = size(flux, /dim)
  node = [s/2, s/2]
  
  arr = [ [[indgen(s(0)), replicate(s(0)-1, s(1)), reverse(indgen(s(0))), replicate(0, s(1))]], $
    [[replicate(0, s(0)), indgen(s(1)), replicate(s(1)-1, s(0)), reverse(indgen(s(1)))]] ]
  arr = reform(arr)
  
  if keyword_set(plot) then begin
    loadct, 0, /sil
    window, 0, xsize = 400, ysize = 400
    !p.position = [0.1, 0.1, 0.95, 0.95]
    plot, indgen(10), /nodata, xran = [0, s(0)], yran = [0, s(1)], xsty = 1, ysty = 1
  endif
  
  text = ['# Region file format: DS9 version 4.1', $
    '# Filename: '+file, $
    'global color=green dashlist=8 3 width=1 font="helvetica 10 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1', $
    'fk5' ]
    
  st = 'polygon(
  
  for i = 0, n_elements(arr(*, 0)) - 1 do begin
  
    endpoint = reform(arr(i, *))
    
    nPoints = ABS(endpoint(0)-node(0)+1) > ABS(endpoint(1)-node(1)+1)
    nPoints *= 5
    xloc = s(0) - 1 - (node(0) + (endpoint(0) - node(0)) * Findgen(nPoints) / (nPoints - 1))
    yloc = node(1) + (endpoint(1) - node(1)) * Findgen(nPoints) / (nPoints - 1)
    
    
    cut = im(xloc, yloc)
    
    ind = where( finite(cut) eq 0, c )
    if c eq 0 then ind = n_elements(cut)-1 else ind = ind(0)
    
    if keyword_set(plot) then begin
      loadct, 0, /sil
      oplot, xloc, yloc, ps = -4, colo = 150
      loadct, 38, /sil
      plots, xloc(ind), yloc(ind), ps = 4, symsize = 2, color = 150
    endif
    
    xy2ad, xloc(ind), yloc(ind), astr, a, d
    
    if i gt 0 then st += ','
    st += a2str(a)+','+a2str(d)
    
  endfor
  
  st += '# color='+color+' width=2 font="helvetica 14 normal roman" text={}'
  text = [text, st]
  
  openw, lun, strtrans(file, '\.fits', '\.reg'), /get_lun, error = err
  if err ne 0 then print, !ERROR_STATE.MSG else begin
    for i = 0, n_elements(text) - 1 do printf, lun, text(i)
    close, /all
  endelse
  
end
