pro pacsman_makereg2, file, color = color

  if not keyword_set(color) then color = 'black'
  
  imag = readfits(file, h, /sil)

  
  header = [ '# Region file format: DS9 version 4.1', $
    '# Filename: '+file, $
    'global color=green dashlist=8 3 width=1 font="helvetica 10 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1', $
    'fk5']
  
              
  pretext = 'polygon(' 
  posttext = ') # color='+color+' width=2 font="helvetica 14 normal roman" text={}'
  
  edges = pacsman_flagedges(imag, h, arr=arr)
  
  ;help, arr
  arrsort = arr
  n = 0
  i = 0
  npoints = n_elements(arr(*, 0))
  remaining = intarr(npoints) + 1
  while n lt npoints do begin
    ;print, i, n
    ind = where( remaining eq 1 )
    ;help, ind
    ind2 = sort( sqrt( (arr(i, 0)-arr(ind, 0))^2. + (arr(i, 1)-arr(ind, 1))^2. ) )
    ind = ind(ind2)
    ;print, ind(0:5)
    arrsort(n, *) = arr(ind(0), *)
    remaining(ind(0)) = 0
    n += 1
    i = ind(0)
  endwhile
  
                                ;create ra dec array
  npoints = n_elements(arr(*, 0))
  adarr = dblarr(npoints, 2)
  extast, h, astr
  for i = 0, npoints - 1 do begin
    xy2ad, arrsort(i, 0), arrsort(i, 1), astr, a, d
    adarr[i, *] = [a, d]
 endfor
                                ;sort to have continuous region
  adarr2 = dblarr(npoints, 2)
  adarr2[0, *] = adarr[0, *]
  n = 0
  distarr = dblarr(npoints)
  while n lt npoints-1 do begin
     gcirc, 2, adarr2[n, 0], adarr2[n, 1], adarr[*, 0], adarr[*, 1], dist
     ind = sort(dist)
     newind = ind[0]
     distarr[n] = dist[newind]
     adarr2[n+1, *] = [adarr[newind, 0], adarr[newind, 1]]
     adarr[newind, *] = !values.f_nan
     n += 1
  endwhile

                                ;split disconnected regions
                                ;look for subsequent pixels with
                                ;distance more than ~2 pixels
  openw, lun, strtrans(file, '\.fits', '\.reg'), /get_lun, error = err  
  if err ne 0 then begin
     print, !ERROR_STATE.MSG
     retall
  endif 
  printf, lun, header

  pxscale = 0.5 * total(astr.cdelt) * 3600.
  print, pxscale
  ind = where( distarr gt 2.*pxscale, c )
  print, c
  
  if c eq 0 then begin 
     text = pretext
     for i = 0, npoints-1 do begin
        text += a2str(adarr2[i, 0])+','+a2str(adarr2[i, 1])
        if i ge 0 and i lt npoints-1 then text += ','
     endfor
     text += posttext
     printf, lun, text
  endif else begin
     startpoint = [0, ind+1]
     endpoint = [ind, npoints-1]
     for j = 0, c do begin        
                                ;make string
        text = ''
        for i = startpoint[j], endpoint[j] do begin
           text += a2str(adarr2[i, 0])+','+a2str(adarr2[i, 1])
           if i ge startpoint[j] and i lt endpoint[j] then text += ','
        endfor
        text = pretext + text + posttext 
        printf, lun, text
     endfor
  endelse

  close, /all
  
end
