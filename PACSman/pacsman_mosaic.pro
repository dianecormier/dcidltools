;______________________________________________________________________________
;used for theoretical line profile (changing sigma for instance)
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
;fringes + baseline + gaussian model
;function LineModel, x, p
;  return, p[0] * exp(-(x-p[1])^2./(2.*p[2]^2.))
;end
;;-----------------------------------___________________________________________
;;polynomial model
;function ContinuumModel, x, p, poly_degree = poly_degree
;  result = 0.
;  result += p[0] * sin( (2.*!pi/p[1])*x + p[2] )
;  for i = 0, poly_degree do result +=  p[i+3]*x^float(i)
;  return, result
;end
;______________________________________________________________________________
;----------------------------------------
pro pacsman_mosaic, lines = lines, cube = cube, dsm = dsm, showcontinuum = showcontinuum

  if not keyword_set(cube) then cube = "cubeAcubeB"
  
  if not keyword_set(lines) then begin
    spawn, 'find . -maxdepth 1 -mindepth 1 -type d ', lines
    for line_i = 0, n_elements(lines) - 1 do lines(line_i) = strtrans(lines(line_i), './', '')
  endif
  
  for line_i = 0, n_elements(lines) - 1 do begin
  
    if strpos(lines(line_i), 'ANALYSIS') gt -1 then continue
    line = lines(line_i)
    print, line
    if not file_test(line+"/"+cube+".sav") then continue
    print, 'building mosaic...'
  
    set_plot, 'PS'
    ;set_plot, 'z'
    ;ERASE
    ;device, set_pixel_depth = 24, decomposed = 0, set_resolution = [800, 800]
    spawn, '\rm - f '+line+'/Mosaic_*.png'
    device, filename = line+'/Mosaic_'+line+'.eps', /encaps, /color, xsize = 20, ysize = 20
    !p.position = [0.1, 0.9, 0.1, 0.9]
    !p.multi = [0, 5, 5, 0, 0]
    !x.margin = 0
    !y.margin = 0
    
    restore, line+"/"+cube+".sav", /verbose
   
    lambda = cube_spectra.lambda - cube_params.lambda_rest  
    flux = cube_spectra.flux / ( 1.e-6 / (9.4*!pi/180./3600.)^2. ) ;MJy sr-1 -> JY
    if cube_params.pacsman_version gt 3.5 then fwhm = cube_params.reference_fwhm else fwhm = cube_params.instrumental_fwhm
    c = cube_params.constraints_continuum*fwhm
    
    deltalambda = cube_params.lambda_rest * sqrt(fwhm^2.+cube_params.broadening(1)^2.) / 3.e5
        
    for x = 1, 5 do for y = 1, 5 do begin
      ;ind = where( lambda gt -c(1) and lambda lt c(2) )
      ind = where( abs(lambda) lt 2.*deltalambda )
      yran = minmax(flux(ind), /nan)
    endfor
    yran(0) -= 0.1
    yran(1) *= 1.5
    
    xran = 5.*deltalambda * [-1., 1.]

    for x = 1, 5 do for y = 1, 5 do begin
    
      if x eq 3 and y eq 3 then tit = line else tit = '''
      
      ind = where( cube_spectra.spaxel_X eq x and cube_spectra.spaxel_Y eq y and cube_spectra.raster eq 1, count )

      lambda = cube_spectra(ind).lambda - cube_params.lambda_rest
      flux = cube_spectra(ind).flux / ( 1.e-6 / (9.4*!pi/180./3600.)^2. ) ;MJy sr-1 -> JY
    
      ind2 = where( abs(lambda) lt 5.*deltalambda )
      
      sc = 0.18
      sh = 0.
      if x eq 2 then sh = 0.1
      !p.position = [ (x-1)*sc+0.05, (y-1)*sc+sh, (x-1)*sc + sc+0.05, (y-1)*sc + sc + sh]
      plot, lambda, flux, title = tit, charsize = 1.5, /nodata, yran = yran, ysty = 1, xsty = 1, xran = xran
      
      loadct, 0, /sil
      plots, [0, 0], yran, line = 1, color = 200
      plots, !x.crange, [0, 0], line = 1, color = 200
      if ~keyword_set(showcontinuum) then begin
        addcomponents = 0.
        indparams_contfringes = 3 + indgen(3+cube_params.poly_degree+1)
        params_cont = cube_fitparams[ind].params[indparams_contfringes]
        oplot, lambda, flux-ContinuumModel(lambda, params_cont, poly_degree = cube_params.poly_degree)
        loadct, 38, /sil
        oplot, lambda[ind2], LineModel(lambda[ind2], cube_fitparams[ind].params[0:2]), color = 200, thick = 4
      endif else begin
        oplot, lambda, flux
        loadct, 38, /sil
        oplot, lambda[ind2], SpectrumModel(lambda[ind2], cube_fitparams[ind].params, poly_degree = cube_params.poly_degree, ncomps = 0), color = 200, thick = 4
      endelse
      ;print, cube_fitparams(ind).params
      legend, /top, /left, '('+a2str(x)+','+a2str(y)+')', box = 0
    
    endfor
    
    ;write_png, line+'/Mosaic_'+line+'.png', tvrd(/true)
    
  endfor
  
  spawn,'\rm -f ' + line+'/Mosaic_'+line+'.png'
  device, /close
  set_plot, 'X'
  !p.multi = 0
  
end