pro merge_cubes, cubelist = cubelist, line = line

  if not keyword_set(cubelist) and not keyword_set(line) then begin
    print, 'Calling sequence.'
    print, ''
    print, " merge_cubes, [line = 'CII157'], [cubelist = 'dir/cube.sav', 'otherdir/cube.sav']]"
    retall
  endif

  print, '___________________________________________________'
  if keyword_set(line) then begin
    cubelist = file_search('./', 'cube*.sav')
    print, 'All cubes: '
    print, cubelist
    ind = where(strpos(cubelist, line) gt -1 and strpos(cubelist, '3x3') eq -1 and strpos(cubelist, '5x5') eq -1 and strpos(cubelist, 'cube_merged') eq -1 and strpos(cubelist, 'before_subtracting') eq -1)
    cubelist = cubelist(ind)
    print, 'Requested line: ', line
  endif

  if keyword_set(cubelist) and ~keyword_set(line) then begin
    labels = ['CII157','OI63','OIII88','OI145','NII122','NII205','NIII57']
    foundline=0
    for j=0,n_elements(labels)-1 do begin
      ind = where(strpos(cubelist, labels[j]) gt -1,ctlab)
      if ctlab gt 0 then begin
       line = labels[j]
       foundline=1
       print, 'Requested line: ', line
      endif
    endfor
    if foundline eq 0 then begin
      print, 'No line found from the list.'
      retall
    endif
  endif
  
  print, 'Files found: ', cubelist
  print, '___________________________________________________'
  
  n = n_elements(cubelist)
  print, 'Reading cube : ', cubelist(0)
  restore, cubelist(0)
  merged_cube_lres = cube_lres
  merged_cube_fitparams = cube_fitparams
  merged_cube_spectra = cube_spectra
  merged_cube_params = cube_params
  
  for i = 1, n - 1 do begin
     ;if strpos(cubelist(i), 'cube_merged') gt -1 then continue
    print, 'Reading cube : ', cubelist(i)
    restore, cubelist(i)
    
    ;increment number of rasters
    print, ''
    print, 'Number of rasters before adding file: ', merged_cube_params.n_rasters
    ind = uniq(merged_cube_lres.raster, sort(merged_cube_lres.raster))
    ;print, merged_cube_lres[ind].raster
    
    
    ;increment number of rasters
    cube_lres.raster += total(merged_cube_params.n_rasters)
    cube_spectra.raster += total(merged_cube_params.n_rasters)
    
    ;merged_cube_params.n_rasters += total(cube_params.n_rasters)
    
    ;possibility to use struct_replace_field.pro
    
    ;need to merge cube_params as well to get the obsid for instance
    merged_cube_lres = [merged_cube_lres, cube_lres]
    merged_cube_fitparams = [merged_cube_fitparams, cube_fitparams]
    merged_cube_params = [merged_cube_params, cube_params]
    merged_cube_spectra = [merged_cube_spectra, cube_spectra]
  endfor
  
   ;increment number of rasters
    print, ''
    print, 'Number of rasters after adding file: ', merged_cube_params.n_rasters
    ind = uniq(merged_cube_lres.raster, sort(merged_cube_lres.raster))
    ;print, merged_cube_lres[ind].raster

  
  
  cube_lres = merged_cube_lres
  cube_fitparams = merged_cube_fitparams
  cube_spectra = merged_cube_spectra
  cube_params = merged_cube_params
    

    if file_test('merged',/directory) eq 0 then file_mkdir, 'merged'
    if file_test('merged/'+line,/directory) eq 0 then file_mkdir, 'merged/'+line
    dir = './merged/'+line+'/'

  
  save, filename = dir+'/cube_merged.sav', cube_lres, cube_fitparams, cube_spectra, cube_params, cube_header
  
end
