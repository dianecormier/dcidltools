pro measure_rms

  dirs = ['./']
  
  ;for dirs_i = 0, n_elements(dirs) - 1 do begin
  ;
  ;    spawn, '\find ' + dirs(dirs_i) + ' -name "*_Flux.fits"', files
  ;    for i = 0, n_elements(files) - 1 do begin
  ;        print, files(i)
  ;        pacsman_makereg, files(i)
  ;    endfor
  ;
  ;endfor
  
  openw, lun, 'rms.txt', /get_lun
  printf, lun, '# OBS SIGMA(MJy/sr) SIGMA(Jy)'
  
  for dirs_i = 0, n_elements(dirs) - 1 do begin
  
    spawn, '\find ' + dirs(dirs_i) + ' -name "*_spectra.fits"', files
    for i = 0, n_elements(files) - 1 do begin
    
      print, files(i)
      
      st = strsplit(files(i), '/', /extract)
      pos = strpos(files(i), st(n_elements(st)-1))
      
      cubefile = strmid(files(i), 0, pos) + 'cube*.sav'
      spawn, '\ls ' + cubefile, cubefile
      cubefile = cubefile(where( strlen(cubefile) eq min(strlen(cubefile)) ) )
      print, cubefile
      if not file_test(cubefile) then continue
      restore, cubefile
      if tag_exist(cube_params, 'pacsman_version') then version = cube_params.pacsman_version else version = cube_params.version
      print, version
      
      f = readfits(files(i), hf, /sil)
      wobs = readfits(files(i), hw, exten_no=1, /sil)
      f = reform(f, n_elements(f))
      wobs = reform(wobs, n_elements(wobs))
      warr = wobs/(1.+cube_params.redshift) - cube_params.lambda_rest
      
      if tag_exist(cube_params, 'instrumental_fwhm') then ref_fwhm = cube_params.instrumental_fwhm else ref_fwhm = cube_params.reference_fwhm
      
      linewidth = sqrt( (ref_fwhm)^2. + (cube_params.broadening(1))^2. ) ;in km s-1
      ;linewidth *= cube_params.lambda_rest*(1.+cube_params.redshift) / 299792.458 ;in um
      
      cube_params.constraints_continuum /= sqrt( cube_params.broadening(1)^2. + ref_fwhm^2. )/ref_fwhm
      
      ;print,cube_params.constraints_continuum(0)*linewidth 
      ind = where( abs(warr) gt cube_params.constraints_continuum(0)*linewidth )
      
      
      signal = f(ind) - median(f(ind), 13)
      
      ;plot, warr(ind), signal
      
      if version ge 3.52 then begin
        sig1 = robust_sigma(signal)
        sig2 = sig1 / (1.e-6 / (9.4*!pi/180./3600.)^2.)
      endif else begin
        sig2 = robust_sigma(signal)
        sig1 = sig2 * (1.e-6 / (9.4*!pi/180./3600.)^2.)
      endelse
      
      print, sig1, sig2
      
      ;correct for continuum?
      ;st = strsplit(files(i), '/', /extract)
      ;file = strtrans(files(i), st(n_elements(st)-1), 'cube*.sav')
      ;spawn, '\ls ' + file, file
      ;file = file(0)
      ;restore, file
      
      printf, lun, cube_params.line, sig1, sig2
      
    endfor
  endfor
  
  close, /all
  
end
