;---------------------------------------------------------------------------------
FUNCTION str_round, input
  x = input
  ;stupid function to round the number in the string...
  if strpos(x, 'NaN') gt -1 or strpos(x, 'Inf') gt -1 or finite(x) eq 0 then return, 'NaN' else begin
    if strpos(x, 'e') gt -1 then begin
      expo = strmid(x, strpos(x, 'e'), strlen(x)-strpos(x, 'e')+1)
      x = strmid(x, 0, strpos(x, 'e'))
    endif else expo = ''
    if strpos(x, '.') eq -1 then x = x + '.'
    x = x + '0000000000000000000'
    if strpos(x, 'e') ge 0 then x = '0.0'
    sx = strsplit(x, '.', /extract)
    ;sx[1] = sx[1] + '000000000000000'
    ;r = nint(strmid(sx[1], 0, 3)/10.)
    ;if (r lt 10) then rst = '0'+sm_a2str(r) else rst = sm_a2str(r)
    r = a2str(long(strmid(sx[0]+sx[1], 0, strlen(sx[0])+3))/1000.)
    if strpos(r, '.') eq -1 then r = r + '.'
    r = r + '0000000000000000000'
    sr = strsplit(r, '.', /extract)
    rst = strmid(sr[1], 0, 3)
    tmp = sr[0] + '.' + rst
    if tmp eq '-0.00' then tmp = '0.00'
    return, tmp + expo
  endelse
END
;______________________________________________________________________________
;print on screen and on file
pro write_message, textarr, lun
  for nn = 0, n_elements(textarr) - 1 do begin
    printf, lun, a2str(textarr[nn])
    print, textarr[nn]
  endfor
end
;############################################################################################################
pro pacsman_footprint, lines = lines, cube = cube, selraster = selraster, redshift = redshift, center = center, sumcomponents = sumcomponents, upperlimits = upperlimits, nsrc = nsrc, background = background, integrated = integrated

  if not keyword_set(cube) then cube = "cubeAcubeB"
  if not keyword_set(nsrc) then nsrc = 1

  if not keyword_set(lines) then begin
    spawn, 'find . -maxdepth 1 -mindepth 1 -type d ', lines
    for line_i = 0, n_elements(lines) - 1 do lines(line_i) = strtrans(lines(line_i), './', '')
  endif
  
  print, lines
  for line_i = 0, n_elements(lines) - 1 do begin
  
    if strpos(lines(line_i), 'ANALYSIS') gt -1 then continue
    line = lines(line_i)
    if not file_test(line+"/"+cube+".sav") then begin
       print, "Not found: " + line+"/"+cube+".sav"
       continue
    endif
    ;individual spectra
    restore, line+"/"+cube+".sav"
    
    val = max(cube_lres.raster)
    if ~keyword_set(selraster) then begin 
      selraster = 1
      ;;-- DC modify 09-10-13
      ;;-- selects the raster with brightest peak (likely better centered?) for multi-raster observations
      ;;-- added keyword 'whichraster' to pacsman_fitpsf as well
      if val ne 1 then begin
        indmax = where(cube_lres.flux eq max(cube_lres.flux,/nan) and cube_lres.flux/cube_lres.error gt 3., ctmax)
        if ctmax gt 0 then selraster = cube_lres(indmax[0]).raster
        print, 'Multi-raster observation: selected brightest raster #',selraster,ctmax
    endif else print, '-- It is a one-raster observation'
    endif else print, 'Manually selected raster: #', selraster
    cube_lres = cube_lres(where(cube_lres.raster eq selraster))

    case nsrc of
       1: begin
          pacsman_fitpsf, line+"/"+cube+".sav", /plot, res = res, status = statusres, whichraster=selraster, background = background, integrated = integrated
          pacsman_fitpsf, line+"/"+cube+".sav", /plot, /cont, res = rescont, status = statusrescont, whichraster=selraster, background = background, integrated = integrated
       end
       2: begin
          pacsman_fittwopsf, line+"/"+cube+".sav", /plot, res = res, status = statusres, whichraster=selraster, background = background, integrated = integrated
          pacsman_fittwopsf, line+"/"+cube+".sav", /plot, /cont, res = rescont, status = statusrescont, whichraster=selraster, background = background, integrated = integrated    
       end
       else: retall
    endcase

    ;pacsman_fitpsf, line+"/"+cube+".sav", /plot, res = res, status = statusres
    ;pacsman_fitpsf, line+"/"+cube+".sav", /plot, /cont, res = rescont, status = statusrescont
    
    ;lambda = 0.
    ;if strpos(line, 'CII157') gt -1 then lambda = 157.7409
    ;if strpos(line, 'OI63') gt -1 then lambda = 63.183705
    ;if strpos(line, 'OI145') gt -1 then lambda = 145.525439
    ;if strpos(line, 'NII122') gt -1 then lambda = 121.89806
    ;if strpos(line, 'OIII52') gt -1 then lambda = 52.
    ;if strpos(line, 'NIII57') gt -1 then lambda = 57.33
    ;if strpos(line, 'OIII88') gt -1 then lambda = 88.356
    ;if strpos(line, 'OH79') gt -1 then lambda = 79.
    ;if lambda eq 0. then stop
    lambda = cube_params.lambda_rest
    if not keyword_set(redshift) then redshift = cube_params.redshift
    lambda *= (1.+ redshift)
    
    PACSman_dir = find_with_def('pacsman_map.pro', !path)
    tmp = strsplit(PACSman_dir, '/', /extract)
    PACSman_dir = '/'
    for i = 0, n_elements(tmp) - 2 do PACSman_dir = PACSman_dir + tmp[i] + '/'
    if PACSman_dir eq '/' then PACSman_dir ='./'
    print, 'Reading corrections in: ', PACSman_dir+"/calib/PACS_PointSourceCorrection.dat"
    readcol, PACSman_dir+"/calib/PACS_PointSourceCorrection.dat", pls_lambda, pls_frac_in1, tmp, pls_frac_in3x3, pls_frac_in5x5, /silent
    
    pls_corr = 1./interpol(pls_frac_in1, pls_lambda, lambda)
    pls_corr = pls_corr[0]
    
    pls3by3_corr = 1./interpol(pls_frac_in3x3, pls_lambda, lambda)
    pls3by3_corr = pls3by3_corr[0]
    
    pls5by5_corr = 1./interpol(pls_frac_in5x5, pls_lambda, lambda)
    pls5by5_corr = pls5by5_corr[0]
    
    ;#######################################################################
    
    ;signal and SNR
    signal_int = cube_lres.integrated_flux
    noise_int = cube_lres.integrated_noise
    
    ;fit
    signal = cube_lres.flux
    noise = cube_lres.error
    if keyword_set(sumcomponents) and tag_exist(cube_lres, 'components_flux') then begin
      for i = 0, 24 do begin
        signal(i) += total(cube_lres(i).components_flux(*), /nan)
        ind = where( finite(cube_lres(i).components_error(*)), c )
        if c gt 0 then noise(i) += total(cube_lres(i).components_error(ind))
      endfor
    endif
    snr = signal / noise
    
    noise_norm = noise / max(signal, /nan)
    signal_norm = signal / max(signal, /nan)
    
    ;continuum
    continuum = cube_lres.continuum
    ;ind = where( continuum lt 0., c )
    ;if c gt 0 then continuum(ind) = 0.
    ;if min(continuum, /nan) lt 0. then continuum += -min(continuum, /nan)
    continuum_noise = cube_lres.continuum_error
    snrcont = continuum / continuum_noise
    
    continuum_noise_norm = continuum_noise / max(continuum, /nan)
    continuum_norm = continuum / max(continuum, /nan)
    
    if keyword_set(center) then begin
      print, 'Using input reference spaxel. '
      indbrightest = where( abs(cube_lres.spaxel_X-center(0)) eq 0 and abs(cube_lres.spaxel_Y-center(1)) eq 0, c )
    endif else begin
      ;brightest spaxel
      if max(snr, /nan) le 1. then begin
        print, 'No detection in any spaxel. Assuming the source is in the central spaxel. '
        indbrightest = where( abs(cube_lres.spaxel_X-3) eq 0 and abs(cube_lres.spaxel_Y-3) eq 0, c )
      endif else begin
        indbrightest = where( signal eq max(signal, /nan), count )
      endelse
    endelse
    
    ;#######################################################################
    ;First method= brightest spaxel and PLS correction
    flux1 = (signal(indbrightest) * pls_corr)(0)
    error1 = (noise(indbrightest) * pls_corr)(0)
    flux1int = (signal_int(indbrightest) * pls_corr)(0)
    error1int = (noise_int(indbrightest) * pls_corr)(0)
    
    flux1_cont = (continuum(indbrightest) * pls_corr)(0)
    error1_cont = (continuum_noise(indbrightest) * pls_corr)(0)
    
    ;#######################################################################
    ;Second method= 3x3 around the brightest
    ind3by3 = where( abs(cube_lres.spaxel_X-cube_lres(indbrightest).spaxel_X) le 1 and abs(cube_lres.spaxel_Y-cube_lres(indbrightest).spaxel_Y) le 1, c )
    
    ;for each quantity, we take at least the uncertainty
    if keyword_set(upperlimits) then begin
      ind = where( signal(ind3by3) lt noise(ind3by3), c )
      if c gt 0 then begin
        signal(ind3by3(ind)) = noise(ind3by3(ind))
      endif
    endif
    flux9 = total( signal(ind3by3), /nan ) * pls3by3_corr
    error9 = total( noise(ind3by3), /nan ) * pls3by3_corr
    
    if keyword_set(upperlimits) then begin
      ind = where( signal_int(ind3by3) lt noise_int(ind3by3), c )
      if c gt 0 then begin
        signal_int(ind3by3(ind)) = noise_int(ind3by3(ind))
      endif
    endif
    flux9int = total( signal_int(ind3by3), /nan ) * pls3by3_corr
    error9int = total( noise_int(ind3by3), /nan ) * pls3by3_corr
    
    if keyword_set(upperlimits) then begin
      ind = where( continuum(ind3by3) lt continuum_noise(ind3by3), c )
      if c gt 0 then begin
        continuum(ind3by3(ind)) = continuum_noise(ind3by3(ind))
      endif
    endif
    flux9_cont = total( continuum(ind3by3), /nan ) * pls3by3_corr
    error9_cont = total( continuum_noise(ind3by3), /nan ) * pls3by3_corr
    
    ind = where( signal gt 2.*noise, ndetected )
    if ndetected eq 0 then begin
      fluxdetected = !values.f_nan
      errordetected = !values.f_nan
      fluxdetectedint = !values.f_nan
      errordetectedint = !values.f_nan
      fluxdetected_cont = !values.f_nan
      errordetected_cont = !values.f_nan
    endif else begin
      fluxdetected = total( signal(ind), /nan )
      errordetected = total( noise(ind), /nan )
      fluxdetectedint = total( signal_int(ind), /nan )
      errordetectedint = total( noise_int(ind), /nan )
      fluxdetected_cont = total( continuum(ind), /nan )
      errordetected_cont = total( continuum_noise(ind), /nan )
    endelse
    
    ;#######################################################################
    ;Third method= 5x5
    ;for each quantity, we take at least the uncertainty
    ;what was done for the 3 by 3 is already corrected
    if keyword_set(upperlimits) then begin
      ind = where( signal lt noise, c )
      if c gt 0 then begin
        signal(ind) = noise(ind)
      endif
    endif
    flux25 = total( signal, /nan ) * pls5by5_corr
    error25 = total( noise, /nan ) * pls5by5_corr
    
    if keyword_set(upperlimits) then begin
      ind = where( signal_int lt noise, c )
      if c gt 0 then begin
        signal_int(ind) = noise_int(ind)
      endif
    endif
    flux25int = total( signal_int, /nan ) * pls5by5_corr
    error25int = total( noise_int, /nan ) * pls5by5_corr
    
    if keyword_set(upperlimits) then begin
      ind = where( continuum lt noise, c )
      if c gt 0 then begin
        continuum(ind) = continuum_noise(ind)
      endif
    endif
    flux25_cont = total( continuum, /nan )* pls5by5_corr
    error25_cont = total( continuum_noise, /nan )* pls5by5_corr

;######################    
    ;combined spectrum
    if file_test(line+"/"+cube+"_3x3.sav") then begin
      restore, line+"/"+cube+"_3x3.sav"
      flux9combint = total(cube_lres(selraster-1).integrated_flux, /nan )
      error9combint = total(cube_lres(selraster-1).integrated_noise, /nan )
      flux9comb = total(cube_lres(selraster-1).flux, /nan)
      error9comb = total(cube_lres(selraster-1).error)
      if keyword_set(sumcomponents) and tag_exist(cube_lres, 'components_flux') then begin
        ;there is only one element for the combined spectra
        flux9comb += total(cube_lres(selraster-1).components_flux(*), /nan)
        ind = where( finite(cube_lres(selraster-1).components_error(*)), c )
        if c gt 0 then error9comb += total(cube_lres(selraster-1).components_error(ind))
      endif
      flux9comb_cont = total(cube_lres(selraster-1).continuum, /nan )
      error9comb_cont = total(cube_lres(selraster-1).continuum_error, /nan )
    endif else begin
      flux9comb = !values.f_nan
      error9comb = !values.f_nan
      flux9combint = !values.f_nan
      error9combint = !values.f_nan
      flux9comb_cont = !values.f_nan
      error9comb_cont = !values.f_nan
    endelse

    ;combined spectrum
    if file_test(line+"/"+cube+"_5x5.sav") then begin
      restore, line+"/"+cube+"_5x5.sav"
      flux25combint = total(cube_lres(selraster-1).integrated_flux, /nan )
      error25combint = total(cube_lres(selraster-1).integrated_noise, /nan )
      flux25comb = total(cube_lres(selraster-1).flux, /nan)
      error25comb = total(cube_lres(selraster-1).error)
      if keyword_set(sumcomponents) and tag_exist(cube_lres, 'components_flux') then begin
        ;there is only one element for the combined spectra
        flux25comb += total(cube_lres(selraster-1).components_flux(*), /nan)
        ind = where( finite(cube_lres(selraster-1).components_error(*)), c )
        if c gt 0 then error25comb += total(cube_lres(selraster-1).components_error(ind))
      endif
      flux25comb_cont = total(cube_lres(selraster-1).continuum, /nan )
      error25comb_cont = total(cube_lres(selraster-1).continuum_error, /nan )
    endif else begin
      flux25comb = !values.f_nan
      error25comb = !values.f_nan
      flux25combint = !values.f_nan
      error25combint = !values.f_nan
      flux25comb_cont = !values.f_nan
      error25comb_cont = !values.f_nan
    endelse


    ;#######################################################################
    ;Results
    
    ;nsigma = str_round(res.flux/res.error)
    ;print, '  - PLS (Optimal)      : ' + str_round(res.flux) + ' +/- ' + str_round(res.error) + ' W m-2 (' + nsigma + ' sigma)'
    
    ;nsigma = str_round(flux1/error1)
    ;print, '  - PLS (fit)          : ' + str_round(flux1) + ' +/- ' + str_round(error1) + ' W m-2 (' + nsigma + ' sigma) (PLS correction factor: '  + str_round(pls_corr) + ')'
    ;nsigma = str_round(flux1int/error1int)
    ;print, '  - PLS (int)          : ' + str_round(flux1int) + ' +/- ' + str_round(error1int) + ' W m-2 (' + nsigma + ' sigma) (PLS correction factor: '  + str_round(pls_corr) + ')'
    
    ;nsigma = str_round(flux9/error9)
    ;print, '  - 3x3 (fit+add)      : ' + str_round(flux9) + ' +/- ' + str_round(error9) + ' W m-2 (' + nsigma + ' sigma) (PLS correction factor: '  + str_round(pls3by3_corr) + ')'
    ;nsigma = str_round(flux9int/error9int)
    ;print, '  - 3x3 (int+add)      : ' + str_round(flux9int) + ' +/- ' + str_round(error9int) + ' W m-2 (' + nsigma + ' sigma) (PLS correction factor: '  + str_round(pls3by3_corr) + ')'
    
    ;nsigma = str_round(fluxdetected/errordetected)
    ;print, '  - detected (fit+add) : ' + str_round(fluxdetected) + ' +/- ' + str_round(errordetected) + ' W m-2 (' + nsigma + ' sigma) (>3-sigma: '+a2str(ndetected)+')'
    ;nsigma = str_round(fluxdetectedint/errordetectedint)
    ;print, '  - detected (int+add) : ' + str_round(fluxdetectedint) + ' +/- ' + str_round(errordetectedint) + ' W m-2 (' + nsigma + ' sigma) (>3-sigma: '+a2str(ndetected)+')'
    
    ;nsigma = str_round(flux9comb/error9comb)
    ;print, '  - 3x3 (comb+fit)     : ' + str_round(flux9comb) + ' +/- ' + str_round(error9comb) + ' W m-2 (' + nsigma + ' sigma) (PLS correction factor: '  + str_round(pls3by3_corr) + ')'
    ;nsigma = str_round(flux9combint/error9combint)
    ;print, '  - 3x3 (comb+int)     : ' + str_round(flux9combint) + ' +/- ' + str_round(error9combint) + ' W m-2 (' + nsigma + ' sigma) (PLS correction factor: '  + str_round(pls3by3_corr) + ')'
    
    ;nsigma = str_round(flux25/error25)
    ;print, '  - 5x5 (fit+add)      : ' + str_round(flux25) + ' +/- ' + str_round(error25) + ' W m-2 (' + nsigma +' sigma)'
    ;nsigma = str_round(flux25int/error25int)
    ;print, '  - 5x5 (int+add)      : ' + str_round(flux25int) + ' +/- ' + str_round(error25int) + ' W m-2 (' + nsigma +' sigma)'
    
    ;nsigma = str_round(flux25comb/error25comb)
    ;print, '  - 5x5 (comb+fit)     : ' + str_round(flux25comb) + ' +/- ' + str_round(error25comb) + ' W m-2 (' + nsigma + ' sigma) (PLS correction factor: '  + str_round(pls5by5_corr) + ')'
    ;nsigma = str_round(flux25combint/error25combint)
    ;print, '  - 5x5 (comb+int)     : ' + str_round(flux25combint) + ' +/- ' + str_round(error25combint) + ' W m-2 (' + nsigma + ' sigma) (PLS correction factor: '  + str_round(pls5by5_corr) + ')'

    
    ;; write results in text file
    openw, lun_tmp, line+'/footprint_results.txt', width = 400, /get_lun
    
    write_message, ['', 'Continuum', ''], lun_tmp
    
    write_message,    '  - RA/DEC               : ' + str_round(double(rescont.position[0])) + '/' + str_round(double(rescont.position[1])), lun_tmp
    nsigma = str_round(rescont.flux[0]/rescont.error[0])
    write_message,    '  - F_opt (Optimal)      : ' + str_round(rescont.flux[0]) + ' +/- ' + str_round(rescont.error[0]) + ' Jy (' + nsigma + ' sigma)', lun_tmp
    
    if n_elements(rescont.flux) eq 2 then begin
       write_message, '  - RA/DEC               : ' + str_round(double(rescont.position[2])) + '/' + str_round(double(rescont.position[3])), lun_tmp
       nsigma2 = str_round(rescont.flux[1]/rescont.error[1])
       write_message, '  - F_opt (Optimal)      : ' + str_round(rescont.flux[1]) + ' +/- ' + str_round(rescont.error[1]) + ' Jy (' + nsigma2 + ' sigma)', lun_tmp
    endif

    nsigma = str_round(flux1_cont/error1_cont)
    write_message,    '  - F_1 (central spaxel) : ' + str_round(flux1_cont) + ' +/- ' + str_round(error1_cont) + ' Jy (' + nsigma + ' sigma) (PLS correction factor: '  + str_round(pls_corr) + ')', lun_tmp
    
    nsigma = str_round(flux9_cont/error9_cont)
    write_message,    '  - F_3x3                : ' + str_round(flux9_cont) + ' +/- ' + str_round(error9_cont) + ' Jy (' + nsigma + ' sigma) (PLS correction factor: '  + str_round(pls3by3_corr) + ')', lun_tmp
    
    nsigma = str_round(fluxdetected_cont/errordetected_cont)
    write_message,    "  - F'_3x3 (detected)    : " + str_round(fluxdetected_cont) + ' +/- ' + str_round(errordetected_cont) + ' Jy (' + nsigma + ' sigma) (>3-sigma: '+a2str(ndetected)+')', lun_tmp
    
    nsigma = str_round(flux9comb_cont/error9comb_cont)
    write_message,    '  - F_3x3 (comb)         : ' + str_round(flux9comb_cont) + ' +/- ' + str_round(error9comb_cont) + ' Jy (' + nsigma + ' sigma) (PLS correction factor: '  + str_round(pls3by3_corr) + ')', lun_tmp
    
    nsigma = str_round(flux25_cont/error25_cont)
    write_message,    '  - F_5x5                : ' + str_round(flux25_cont) + ' +/- ' + str_round(error25_cont) + ' Jy (' + nsigma +' sigma)', lun_tmp
    
   
    write_message, ['', 'Line', ''], lun_tmp

    write_message, line, lun_tmp
    write_message, 'redshift = ' + a2str(redshift), lun_tmp
    
    write_message,    '  - RA/DEC               : ' + str_round(double(res.position[0])) + '/' + str_round(double(res.position[1])), lun_tmp
    nsigma = str_round(res.flux[0]/res.error[0])
    write_message,    '  - F_opt (Optimal)      : ' + str_round(res.flux[0]) + ' +/- ' + str_round(res.error[0]) + ' W m-2 (' + nsigma + ' sigma)', lun_tmp
    
    if n_elements(res.flux) eq 2 then begin
       write_message, '  - RA/DEC               : ' + str_round(double(res.position[2])) + '/' + str_round(double(res.position[3])), lun_tmp
       nsigma2 = str_round(res.flux[1]/res.error[1])
       write_message, '  - F_opt (Optimal)      : ' + str_round(res.flux[1]) + ' +/- ' + str_round(res.error[1]) + ' W m-2 (' + nsigma2 + ' sigma)', lun_tmp
    endif

    nsigma = str_round(flux1/error1)
    write_message,    '  - F_1 (central spaxel) : ' + str_round(flux1) + ' +/- ' + str_round(error1) + ' W m-2 (' + nsigma + ' sigma) (PLS correction factor: '  + str_round(pls_corr) + ')', lun_tmp
    nsigma = str_round(flux1int/error1int)
    write_message,    '  - F_1 (integrated)     : ' + str_round(flux1int) + ' +/- ' + str_round(error1int) + ' W m-2 (' + nsigma + ' sigma) (PLS correction factor: '  + str_round(pls_corr) + ')', lun_tmp
    
    nsigma = str_round(flux9/error9)
    write_message,    '  - F_3x3 (fit+add)      : ' + str_round(flux9) + ' +/- ' + str_round(error9) + ' W m-2 (' + nsigma + ' sigma) (PLS correction factor: '  + str_round(pls3by3_corr) + ')', lun_tmp
    nsigma = str_round(flux9int/error9int)
    write_message,    '  - F_3x3 (int+add)      : ' + str_round(flux9int) + ' +/- ' + str_round(error9int) + ' W m-2 (' + nsigma + ' sigma) (PLS correction factor: '  + str_round(pls3by3_corr) + ')', lun_tmp
    
    nsigma = str_round(fluxdetected/errordetected)
    write_message,    "  - F'_3x3 (fit+add)     : " + str_round(fluxdetected) + ' +/- ' + str_round(errordetected) + ' W m-2 (' + nsigma + ' sigma) (>3-sigma: '+a2str(ndetected)+')', lun_tmp
    nsigma = str_round(fluxdetectedint/errordetectedint)
    write_message,    "  - F'_3x3 (int+add)     : " + str_round(fluxdetectedint) + ' +/- ' + str_round(errordetectedint) + ' W m-2 (' + nsigma + ' sigma) (>3-sigma: '+a2str(ndetected)+')', lun_tmp
    
    nsigma = str_round(flux9comb/error9comb)
    write_message,    '  - F_3x3 (comb+fit)     : ' + str_round(flux9comb) + ' +/- ' + str_round(error9comb) + ' W m-2 (' + nsigma + ' sigma) (PLS correction factor: '  + str_round(pls3by3_corr) + ')', lun_tmp
    nsigma = str_round(flux9combint/error9combint)
    write_message,    '  - F_3x3 (comb+int)     : ' + str_round(flux9combint) + ' +/- ' + str_round(error9combint) + ' W m-2 (' + nsigma + ' sigma) (PLS correction factor: '  + str_round(pls3by3_corr) + ')', lun_tmp
    
    nsigma = str_round(flux25/error25)
    write_message,    '  - F_5x5 (fit+add)      : ' + str_round(flux25) + ' +/- ' + str_round(error25) + ' W m-2 (' + nsigma +' sigma)', lun_tmp
    nsigma = str_round(flux25int/error25int)
    write_message,    '  - F_5x5 (int+add)      : ' + str_round(flux25int) + ' +/- ' + str_round(error25int) + ' W m-2 (' + nsigma +' sigma)', lun_tmp
 
    nsigma = str_round(flux25comb/error25comb)
    write_message,    '  - F_5x5 (comb+fit)     : ' + str_round(flux25comb) + ' +/- ' + str_round(error25comb) + ' W m-2 (' + nsigma + ' sigma) (PLS correction factor: '  + str_round(pls5by5_corr) + ')', lun_tmp
    nsigma = str_round(flux25combint/error25combint)
    write_message,    '  - F_5x5 (comb+int)     : ' + str_round(flux25combint) + ' +/- ' + str_round(error25combint) + ' W m-2 (' + nsigma + ' sigma) (PLS correction factor: '  + str_round(pls5by5_corr) + ')', lun_tmp
    
    close, lun_tmp


    result = { line: line, integrated: keyword_set(integrated), $
      brightest: [nint(cube_lres(indbrightest).spaxel_X), nint(cube_lres(indbrightest).spaxel_Y)], $
      ;pls_proba: pls_proba, pls_proba_acc: pls_proba_acc, $
      ;plscont_proba: plscont_proba, plscont_proba_acc: plscont_proba_acc, $
      ;shifted: shifted, shifted_cont: shifted_cont, $
      flux1: double(flux1), error1: double(error1), $
      flux9: double(flux9), error9: double(error9), $
      fluxdetected: double(fluxdetected), errordetected: double(errordetected), $
      ndetected: ndetected, $
      flux9comb: double(flux9comb), error9comb: double(error9comb), $
      flux25: double(flux25), error25: double(error25), $
      flux25comb: double(flux25comb), error25comb: double(error25comb), $
      flux1int: double(flux1int), error1int: double(error1int), $
      flux9int: double(flux9int), error9int: double(error9int), $
      fluxdetectedint: double(fluxdetectedint), errordetectedint: double(errordetectedint), $
      flux9combint: double(flux9combint), error9combint: double(error9combint), $
      flux25int: double(flux25int), error25int: double(error25int), $
      flux25combint: double(flux25combint), error25combint: double(error25combint), $
      flux1_cont: double(flux1_cont), error1_cont: double(error1_cont), $
      flux9_cont: double(flux9_cont), error9_cont: double(error9_cont), $
      fluxdetected_cont: double(fluxdetected_cont), errordetected_cont: double(errordetected_cont), $
      flux9comb_cont: double(flux9comb_cont), error9comb_cont: double(error9comb_cont), $
      flux25_cont: double(flux25_cont), error25_cont: double(error25_cont), $
      flux25comb_cont: double(flux25comb_cont), error25comb_cont: double(error25comb_cont), $
      fluxpsf: double(res.flux), errorpsf: double(res.error), chi2psf: res.redchi2, $
      fluxpsf_cont: double(rescont.flux), errorpsf_cont: double(rescont.error), chi2psf_cont: rescont.redchi2, $
      psfbrightest: res.brightest, psfbrightest_cont: rescont.brightest, $
      psfposition: res.position $
      }
      
    save, filename = line+'/footprint_results.sav', result
    spawn, '\rm -f '+line+'/extent_test.*'
    

 print, ''
    print, 'Results written in ' + line+'/footprint_results.sav'
 


  endfor
  
  
end
