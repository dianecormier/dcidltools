;______________________________________________________________________________
;fringes + baseline + gaussian models
function SpectrumModel, x, p, poly_degree = poly_degree, ncomps = ncomps
  
  result = 0.
  
  result += p[0] * exp(-(x-p[1])^2./(2.*p[2]^2.)) 

  if ncomps gt 0 then for nci = 0, ncomps - 1 do result += p[3+3*nci] * exp(-(x-p[3+3*nci+1])^2./(2.*p[3+3*nci+2]^2.))
  
  result += p[3+3*ncomps] * sin( (2.*!pi/p[3+3*ncomps+1])*x + p[3+3*ncomps+2] )
  
  for i = 0, poly_degree do result +=  p[i+3+3*ncomps+2+1]*x^float(i)
  
  return, result
end
;______________________________________________________________________________
;gaussian model
function LineModel, x, p
  return, p[0] * exp(-(x-p[1])^2./(2.*p[2]^2.))
end
;-----------------------------------___________________________________________
;polynomial model
function ContinuumModel, x, p, poly_degree = poly_degree
  result = 0.
  result += p[0] * sin( (2.*!pi/p[1])*x + p[2] )
  for i = 0, poly_degree do result +=  p[i+3]*x^float(i)
  return, result
end
;______________________________________________________________________________
pro pacsman_convert, lines = lines, datacloud = datacloud, dsm = dsm

  if not keyword_set(lines) then begin
    spawn, 'find . -type d -maxdepth 1 -mindepth 1', lines
    for line_i = 0, n_elements(lines) - 1 do lines(line_i) = strtrans(lines(line_i), './', '')
  endif
  
  print, lines
  
  for line_i = 0, n_elements(lines) - 1 do begin
  
    line = lines(line_i)
    
    spawn, '\rm '+line+'/cubeAcubeB*.dat'
    
    file = line+"/cubeAcubeB.sav"
    if not file_test(file) then continue
    
    ;if tag_exist(cube_params, 'pacsman_version') then begin
    ;  if cube_params.pacsman_version ge 3.55 then ncomps = n_elements(cube_lres(0).components_flux)
    ;endif else ncomps = 0
    ;help, ncomps
    
    spawn, '\rm -f ' + strtrans(file, '.sav', '_spectra.dat.gz')
    restore, file, /verbose
    openw, lun, strtrans(file, '.sav', '_spectra.dat'), width = 800, /get_lun
    printf, lun, '# line ' + cube_params.line
    printf, lun, '# band ' + cube_params.band
    printf, lun, '# order ' + a2str(cube_params.order)
    printf, lun, '# reference fwhm ' + a2str(cube_params.reference_fwhm)
    printf, lun, '# lambda_rest ' + a2str(cube_params.lambda_rest)
    printf, lun, '# redshift ' + a2str(cube_params.redshift)
    printf, lun, '# fringes ' + a2str(cube_params.fringes)
    printf, lun, '# poly_degree ' + a2str(cube_params.poly_degree)
    printf, lun, '# range ' + a2str(cube_params.range)
    if tag_exist(cube_params, 'pacsman_version') then if cube_params.pacsman_version ge 3.55 then printf, lun, '# addcomponents ' + a2str(cube_params.addcomponents)
    printf, lun, '# OBS_ID ' + a2str(cube_params.obsid)
    printf, lun, '# version ' + a2str(cube_params.hipe_version)
    printf, lun, '# COLUMNS: '
    printf, lun, '# raster x y lambda flux error norm_flux cont_sub_flux lambda_obs fit'
    for x = 1, 5 do for y = 1, 5 do for raster = 1, max(cube_spectra.raster) do begin
      ind = where( cube_spectra.spaxel_X eq x and cube_spectra.spaxel_Y eq y and cube_spectra.raster eq raster, c )
      specfit = SpectrumModel(cube_spectra(ind).lambda - cube_params.lambda_rest, cube_fitparams(ind).params, poly_degree = cube_params.poly_degree, ncomps = 0)
      specfit *= ( 1.e-6 / (9.4*!pi/180./3600.)^2. ) 
      for i = 0, n_elements(cube_spectra(ind).lambda) - 1 do begin
        if not finite(cube_spectra(ind).lambda(i)) then continue
        printf, lun, raster, x, y, cube_spectra(ind).lambda(i), cube_spectra(ind).flux(i), cube_spectra(ind).error(i), $
          cube_spectra(ind).normalized_flux(i), cube_spectra(ind).cont_subtracted_flux(i), cube_spectra(ind).lambda_obs(i), specfit(i)
      endfor
    endfor
    close, /all
    ;spawn, '\gzip -f ' + strtrans(file, '.sav', '_spectra.dat')
    ;spawn, '\rm -f ' + strtrans(file, '.sav', '_spectra.dat')
    
    
    
    
    openw, lun, strtrans(file, '.sav', '_results.dat'), width = 800, /get_lun
    printf, lun, '# line ' + cube_params.line
    printf, lun, '# band ' + cube_params.band
    printf, lun, '# order ' + a2str(cube_params.order)
    printf, lun, '# reference fwhm ' + a2str(cube_params.reference_fwhm)
    printf, lun, '# lambda_rest ' + a2str(cube_params.lambda_rest)
    printf, lun, '# redshift ' + a2str(cube_params.redshift)
    printf, lun, '# fringes ' + a2str(cube_params.fringes)
    printf, lun, '# poly_degree ' + a2str(cube_params.poly_degree)
    printf, lun, '# range ' + a2str(cube_params.range)
    ;printf, lun, '# addlocal ' + a2str(cube_params.addlocal)
    if tag_exist(cube_params, 'pacsman_version') then if cube_params.pacsman_version ge 3.55 then printf, lun, '# addcomponents ' + a2str(cube_params.addcomponents)
    printf, lun, '# COLUMNS: '
    printf, lun, '# raster x y flux error int_flux int_error cont cont_error velo velo_error fwhm fwhm_error ra dec chisq'
    for x = 1, 5 do for y = 1, 5 do for raster = 1, max(cube_spectra.raster) do begin
      ind = where( cube_lres.spaxel_X eq x and cube_lres.spaxel_Y eq y and cube_lres.raster eq raster, c )
      printf, lun, raster, x, y, cube_lres(ind).flux, cube_lres(ind).error, cube_lres(ind).integrated_flux, cube_lres(ind).integrated_noise, $
        cube_lres(ind).continuum, cube_lres(ind).continuum_error, cube_lres(ind).velocity, cube_lres(ind).velocity_error, cube_lres(ind).fwhm, $
        cube_lres(ind).fwhm_error, cube_lres(ind).ra, cube_lres(ind).dec, cube_lres(ind).chisq
    endfor
    close, /all
    spawn, '\gzip -f ' + strtrans(file, '.sav', '_results.dat')
    spawn, '\rm -f ' + strtrans(file, '.sav', '_results.dat')
     
     
     
     
    if keyword_set(datacloud) then begin
    
      for raster = 1, max(cube_spectra.raster) do begin
        file = line+"/cubeAcubeB_raster"+a2str(raster)+"_datacloud.sav"
        restore, file, /verbose
        openw, lun, strtrans(file, '.sav', '.dat'), width = 800, /get_lun
        printf, lun, '# line ' + cube_params.line
        printf, lun, '# band ' + cube_params.band
        printf, lun, '# order ' + a2str(cube_params.order)
        printf, lun, '# reference fwhm ' + a2str(cube_params.reference_fwhm)
        printf, lun, '# lambda_rest ' + a2str(cube_params.lambda_rest)
        printf, lun, '# redshift ' + a2str(cube_params.redshift)
        printf, lun, '# fringes ' + a2str(cube_params.fringes)
        printf, lun, '# poly_degree ' + a2str(cube_params.poly_degree)
        printf, lun, '# range ' + a2str(cube_params.range)
        ;printf, lun, '# addlocal ' + a2str(cube_params.addlocal)
        if tag_exist(cube_params, 'pacsman_version') then if cube_params.pacsman_version ge 3.55 then printf, lun, '# addcomponents ' + a2str(cube_params.addcomponents)
        printf, lun, '# COLUMNS: '
        printf, lun, '# x y lambda flux error lambda_obs fit mask ra dec'
        for x = 1, 5 do for y = 1, 5 do begin
          ind = where( datacloud.spaxel_X eq x and datacloud.spaxel_Y eq y )
          specfit = SpectrumModel(datacloud(ind).lambda - cube_params.lambda_rest, cube_fitparams(ind).params, poly_degree = cube_params.poly_degree, ncomps = 0)
          specfit *= ( 1.e-6 / (9.4*!pi/180./3600.)^2. ) 
          for i = long(0), n_elements(datacloud(ind).lambda) - 1 do begin
            if not finite(datacloud(ind).lambda(i)) then continue
            printf, lun, x, y, datacloud(ind).lambda(i), datacloud(ind).flux(i), datacloud(ind).error(i), $
              datacloud(ind).lambda_obs(i), specfit(i), datacloud(ind).mask(i), datacloud(ind).ra(i), datacloud(ind).dec(i)
          endfor
        endfor
        close, /all
        spawn, '\gzip -f ' + strtrans(file, '.sav', '.dat')
        spawn, '\rm -f ' + strtrans(file, '.sav', '.dat')
      endfor
      
    endif
    
  endfor
  
  
end