;##############################################################################
;##############################################################################
pro pacsman_fittwopsf_scale, templates, footprint_obs = footprint_obs, footprint_errobs = footprint_errobs, allplot = allplot, res = res, $
    status = status, background = background
    
  bestres = 1.e18
  if keyword_set(background) then degree = 0 else degree = -1
  
  for ti = 0, n_elements(templates) - 1 do begin        
     
     if abs(templates[ti].offsetx) gt 24. or abs(templates[ti].offsety) gt 24. then continue

     for tj = ti, n_elements(templates) - 1 do begin



     if abs(templates[tj].offsetx) gt 24. or abs(templates[tj].offsety) gt 24. then continue
    ;if templates[ti].offsetx ne 8. or templates[ti].offsety ne 4. then continue
    footprint1 = reform(templates[ti].footprint)
    footprint2 = reform(templates[tj].footprint)
    footprint = footprint1 + footprint2

    if keyword_set(allplot) then begin
    endif
    
    obs = reform(footprint_obs, 25)
    err = reform(footprint_errobs, 25)
    modelpsf = reform([[[reform(footprint1, 25)]],[[reform(footprint2, 25)]]])
    model = dblarr(n_elements(obs), 2+degree+1) ;matrix with 1 psf and background coefs
    
    w = fltarr(25)
    for i = 0, n_elements(obs) - 1 do begin
      w[i] = 1./err[i] ;err should be err^2 but pb with float, just take sig^2 as the final error at the end
      if (finite(obs[i]) eq 0) or (finite(err[i]) eq 0) then begin
        w[i] = double(0.)
        obs[i] = double(0.) ;mregress doesn't like nans
      endif
      model(i, 0) = double( modelpsf(i, 0) )
      model(i, 1) = double( modelpsf(i, 1) )
      if keyword_set(background) then model[i, 2] = double( 1. ) ;that's the background
    endfor
    coeffs = sm_mregress(model, obs, weights = w, status = stat, sigma = sig)
    
    a = coeffs[0:1] ;PSF scaling factors
    e = sig[0:1]^2. ;errors
    
    if keyword_set(allplot) then begin
   endif

    if keyword_set(background) then bg = coeffs(2) else bg = 0.
    model2 = a[0]*model[*, 0]+a[1]*model[*, 1]+bg

    chi2 = total( (abs(obs-median(obs))-model2)^2. / err^2. ) ;/ 1. 1 ;degree of freedom
    residual = chi2
    ;residual = total( abs( abs(obs-median(obs)) - model*a ) * err, /nan ) / total(err, /nan)
    ;residual = total( abs( abs(obs) - model*a ) * err, /nan ) / total(err, /nan)
    ;print, residual
    
    if residual lt bestres and a[0] gt 0 and a[1] gt 0 then begin
       if keyword_set(allplot) then begin
          set_plot, 'X'
          wset, 0
          loadct, 0, /sil
          !p.position = [0.1, 0.5, 0.45, 0.95]
          pacsman_disp, footprint_obs, /sq, title = 'Observed Footprint' ;, /log
          for x = 1, 4 do plots, x*[1, 1], [0, 5], color = 100
          for y = 1, 4 do plots, [0, 5], y*[1, 1], color = 100
          for x = 0, 4 do for y = 0, 4 do if footprint_obs(x, y) lt 2.*footprint_errobs(x, y)  then begin
                                ;plots, x+.5, y+.5, ps = 4
             plots, [x, x+1], [y, y+1], color = 100
             plots, [x+1, x], [y, y+1], color = 100
          endif
                                ;!p.position = [0.5, 0.5, 0.95, 0.95]
                                ;pacsman_disp, footprint, /sq, /noerase, title = 'Template Footprint';, /log
                                ;for x = 1, 4 do plots, x*[1, 1], [0, 5], color = 100
                                ;for y = 1, 4 do plots, [0, 5], y*[1, 1], color = 100
          
          loadct, 0, /sil
          !p.position = [0.5, 0.5, 0.95, 0.95]
          pacsman_disp, reform(model2, 5, 5), /sq, /noerase, title = 'Template Footprint' ;, /log
          for x = 1, 4 do plots, x*[1, 1], [0, 5], color = 100
          for y = 1, 4 do plots, [0, 5], y*[1, 1], color = 100
          
          !p.position = [0.17, 0.15, 0.95, 0.45]
          loadct, 0, /sil
          ploterror, obs, err, ps = -4,thick = 1, xtit = 'Array index', ytit = 'Flux scaling factor', /noerase
                                ;ind_good = where( obs gt 2.*err )
          loadct, 38, /sil
          oploterror, obs, err, ps = 4, color = 200, errcolor = 200, symsize = 2
          plots, !x.crange, [0, 0], line = 0
          oplot, model2, color = 100, ps = -4, thick = 2
          plots, !x.crange, replicate(bg, 2), line = 1
          wait, 0.1
       endif
      bestmodel = ti
      bestres = residual
      scalingfactor = a
      scalingfactorerror = e
      bestpos = [templates(ti).offsetx, templates(ti).offsety, templates(tj).offsetx, templates(tj).offsety]
      ;print, bestres
      ;print, bestpos
    endif
    
  endfor
  endfor
  
  if bestres eq 1.e18 then begin
    print, 'unable to converge'
    status = -1
    return
  endif
  
  status = 0
  
  ;order sources by distance to center
  dist0 = sqrt(total(bestpos[0:1]^2.))
  dist1 = sqrt(total(bestpos[2:3]^2.))
  
  if dist1 le dist0 then begin
     scalingfactor = reverse(scalingfactor)
     scalingfactorerror = reverse(scalingfactorerror)
     bestpos = [bestpos[2:3], bestpos[0:1]]
  endif

  print, '--------------------------------'
  print, 'Best model   : ', bestmodel
  print, 'Reduced Chi2 : ', bestres
  print, 'sqrt(Chi2)   : ', sqrt(bestres)
  print, '--------------------------------'
  print, 'Background   : ', bg
  print, '--------------------------------'
  print, 'Offset       : ', bestpos
  print, 'Flux         : ', scalingfactor
  print, 'Error        : ', scalingfactorerror
  print, '--------------------------------'
  
  ;print, max(scalingfactor*bestmodel)
  
  res = { redchi2: sqrt(bestres), offset: bestpos, flux: scalingfactor, error: scalingfactorerror, brightest: [0, 0], position: dblarr(4), background: bg }
  
end
;----------------------------------------------------------
;----------------------------------------------------------
;----------------------------------------------------------
pro pacsman_fittwopsf, file, removebg = removebg, plot = plot, allplot = allplot, res = res, cont = cont, status = status, background = background, whichraster = whichraster, integrated = integrated

  resnan = { redchi2: !values.f_nan, offset: !values.f_nan, offset2nd: !values.f_nan, flux: !values.f_nan, error: !values.f_nan, brightest: [!values.f_nan, !values.f_nan], position: [!values.f_nan, !values.f_nan], background: !values.f_nan }
  
  if ~keyword_set(whichraster) then whichraster = 1
  
  ;----------------------------------------------------------
  ;prepare observation
  if not keyword_set(file) then file = 'cubeAcubeB.sav'
  restore, file
  if keyword_set(whichraster) then cube_lres = cube_lres(where(cube_lres.raster eq whichraster)) ;;- DC modify 09-10-13
  lambda = cube_params.lambda_rest * (1.+cube_params.redshift)
  footprint_obs = fltarr(5, 5)
  footprint_errobs = footprint_obs
  for x = 1, 5 do for y = 1, 5 do begin
    ind = where( cube_lres.spaxel_X eq x and cube_lres.spaxel_Y eq y )
    if keyword_set(cont) then begin
      footprint_obs(x-1, y-1) = cube_lres(ind).continuum
      footprint_errobs(x-1, y-1) = cube_lres(ind).continuum_error
    endif else begin
       if keyword_set(integrated) then begin
          footprint_obs(x-1, y-1) = cube_lres(ind).integrated_flux
          footprint_errobs(x-1, y-1) = cube_lres(ind).integrated_noise
       endif else begin
          footprint_obs(x-1, y-1) = cube_lres(ind).flux
          footprint_errobs(x-1, y-1) = cube_lres(ind).error
       endelse
    endelse
    if footprint_errobs(x-1, y-1) le 0. then begin
      print, 'Error is <= 0, replacing with median error'
      footprint_errobs(x-1, y-1) = median( cube_lres.error )
    endif
  endfor
  
  ind_good = where( footprint_obs gt 2.*footprint_errobs, c_good, complement = ind_bad )
  print, 'Detected spaxels (>2 sigma) : ', c_good
  if keyword_set(removebg) then begin
    footprint_obs -= median(footprint_obs)
  endif
  
  if keyword_set(allplot) then begin
    set_plot, 'X'
    window, 0, xsize = 900, ysize = 800
  endif
  
  psffile = filepath(root_dir = ProgramRootDir(), subdir = 'calib/', 'PSF_5x5_templates.sav')
  restore, psffile
  
  ;----------------------------------------------------------
  ;prepare ans scale PSF
  ;lambda = 160.
  l = templates[UNIQ(templates.lambda, SORT(abs(templates.lambda-lambda)))].lambda
  print, lambda, l(0), l(1)
  w = 1./abs(lambda - l(0:1))
  ind = where( finite(w, /inf), c, complement = noind )
  if c gt 0 then begin
    w(ind) = 1.
    w(noind) = 0.
  endif
  
  print, '--------> ' + a2str(l(0))
  ind0 = where( templates.lambda eq l(0) )
  pacsman_fittwopsf_scale, templates(ind0), footprint_obs = footprint_obs, footprint_errobs = footprint_errobs, allplot = allplot, res = res, status = status, background = background
  if status eq -1 then begin
    res = resnan
    return
  endif
  res0 = res
  
  print, '--------> ' + a2str(l(1))
  ind1 = where( templates.lambda eq l(1) )
  pacsman_fittwopsf_scale, templates(ind1), footprint_obs = footprint_obs, footprint_errobs = footprint_errobs, allplot = allplot, res = res, status = status, background = background
  if status eq -1 then begin
    res = resnan
    return
  endif
  res1 = res  

  print, '--------> Final'
  fluxfinal = ( w(0)*res0.flux + w(1)*res1.flux ) / total(w)
  errorfinal = ( w(0)*res0.error + w(1)*res1.error ) / total(w)
  bgfinal = ( w(0)*res0.background + w(1)*res1.background ) / total(w)
  print, 'Source flux      : ', fluxfinal[0], errorfinal[0]
  print, 'Source flux (#2) : ', fluxfinal[1], errorfinal[1]

  res = res0
  res.flux = fluxfinal
  res.error = errorfinal
  res.background = bgfinal
  res.redchi2 = min([res0.redchi2, res1.redchi2])
  ;res.offset = [avg([res0.offset(0), res1.offset(0)]), avg([res0.offset(1), res1.offset(1)]), avg([res0.offset(2), res1.offset(2)]), avg([res0.offset(3), res1.offset(3)])]
  res.offset = res0.offset
  
  pa = -cube_params.position_angle * !pi/180.
  rot_matrix = [ [cos(pa), -sin(pa)], [sin(pa), cos(pa)] ]
  print, res.offset
  ;source 1
  res.position[0:1] = (res.offset[0:1]/3600.) # rot_matrix
  indcs = where( cube_lres.spaxel_X eq 3 and cube_lres.spaxel_Y eq 3 )
  res.position(0) /= cos((res.position(1)+cube_lres(indcs).dec)*!pi/180.)
  res.position[0:1] += [cube_lres(indcs).ra, cube_lres(indcs).dec]
  print, 'Source centroid : ', res.position[0:1]
  ;source 2
  res.position[2:3] = (res.offset[2:3]/3600.) # rot_matrix
  indcs = where( cube_lres.spaxel_X eq 3 and cube_lres.spaxel_Y eq 3 )
  res.position(2) /= cos((res.position(3)+cube_lres(indcs).dec)*!pi/180.)
  res.position[2:3] += [cube_lres(indcs).ra, cube_lres(indcs).dec]
  print, 'Source centroid : ', res.position[2:3]
  
  ;----------------------------------------------------------
  ;plots
  if keyword_set(plot) or keyword_set(allplot) then begin
  
    set_plot, 'PS'
    ncolors = 256
    if keyword_set(cont) then plotfile = strtrans(file, '.sav', '_Opt_cont.eps') else plotfile = strtrans(file, '.sav', '_Opt.eps')
    device, filename = plotfile, /encaps, xsize = 20, ysize = 17, /color
    !p.font = 0
    !p.charsize = 1.2
    ;set_plot, 'Z'
    ;ERASE
    ;device, set_resolution=[600, 500], set_pixel_depth = 24, decomposed = 0
    
    ;test
    indf0 = where( templates(ind0).offsetx eq res0.offset(0) and templates(ind0).offsety eq res0.offset(1), c0 )
    indf1 = where( templates(ind1).offsetx eq res1.offset(0) and templates(ind1).offsety eq res1.offset(1), c1 )
    if c0*c1 ne 1 then begin
      print, 'more than 1 minimum?!'
      print, c
      retall
    endif
    indg0 = where( templates(ind0).offsetx eq res0.offset(2) and templates(ind0).offsety eq res0.offset(3), c0 )
    indg1 = where( templates(ind1).offsetx eq res1.offset(2) and templates(ind1).offsety eq res1.offset(3), c1 )
    if c0*c1 ne 1 then begin
      print, 'more than 1 minimum?!'
      print, c
      retall
    endif
    
    ;plot footprints
    loadct, 0, /sil
    !p.position = [0.1, 0.5, 0.45, 0.95]
    pacsman_disp, footprint_obs, /sq, title = 'Observed Footprint';, /log
    for x = 1, 4 do plots, x*[1, 1], [0, 5], color = 200
    for y = 1, 4 do plots, [0, 5], y*[1, 1], color = 200
    
    ;footprint_model = (w(0)*templates(ind0[indf0]).footprint+w(1)*templates(ind1[indf1]).footprint)/total(w)
    footprint_model = ( w(0) * (res0.flux[0]*templates(ind0[indf0]).footprint+res1.flux[1]*templates(ind0[indg0]).footprint) + $
                       w(1) * (res0.flux[0]*templates(ind1[indf1]).footprint+res1.flux[1]*templates(ind1[indg1]).footprint) ) / total(w) + res.background

    writefits, strtrans(plotfile, '.eps', '.fits'), [[[footprint_obs]], [[footprint_model]], [[footprint_model-footprint_obs]]]

    for x = 0, 4 do for y = 0, 4 do begin
      if footprint_obs(x, y) lt 2.*footprint_errobs(x, y)  then begin
        ;plots, x+.5, y+.5, ps = 4
        plots, [x, x+1], [y, y+1], color = 200
        plots, [x+1, x], [y, y+1], color = 200
      endif
      if footprint_model(x, y) eq max(footprint_model) then res.brightest = [x+1, y+1]
    endfor
    !p.position = [0.5, 0.5, 0.95, 0.95]
    pacsman_disp, footprint_model, /sq, /noerase, title = 'Template Footprint';, /log
    for x = 1, 4 do plots, x*[1, 1], [0, 5], color = 200
    for y = 1, 4 do plots, [0, 5], y*[1, 1], color = 200
    
    ;plot curves
    obs = reform(footprint_obs, 25)
    err = reform(footprint_errobs, 25)
    model0 = reform(res0.flux[0]*templates(ind0(indf0)).footprint, 25) + reform(res0.flux[1]*templates(ind0(indg0)).footprint, 25)
    model1 = reform(res1.flux[0]*templates(ind1(indf1)).footprint, 25) + reform(res1.flux[1]*templates(ind1(indg1)).footprint, 25) 
    arr = reverse(sort(model0))
    obs = obs(arr)
    err = err(arr)
    model0 = model0(arr)
    model1 = model1(arr)
    ind_good = where( obs gt 2.*err, c_good, complement = ind_bad )
    model = model0 * 0.
    for i = 0, n_elements(model0) - 1 do model(i) = ( w(0)*model0(i) + w(1)*model1(i) ) / total(w) + res.background
    
    !p.position = [0.17, 0.10, 0.87, 0.45]
    xarr = 1 + findgen(25)
    xran = [0, 25.5]
    ind = where( obs gt 0. ,ct0 )
    
    ;if ct0 gt 0 then yran = minmax([obs(ind)+err(ind), model*(fluxfinal+errorfinal)], /nan) else $
    ;                 yran = [0, max(model*(fluxfinal+errorfinal), /nan)]  
    yran = minmax([obs(ind)+err(ind)], /nan)
    
    if keyword_set(cont) then ylog = 0 else ylog = 1
    
    plot, indgen(10), ps = -4,thick = 1, xtit = 'Array index', ytit = 'Flux scaling factor', $
      /nodata, /noerase, xran = xran, xsty = 1, yran = yran, ylog = ylog, ysty = 1+8
    axis, yaxis = 1, ysty = 1, yran = yran / (9.4*!pi/180./3600.)^2., color = 0, ytitle = '[W m'+textoidl('^{-2}')+' sr'+textoidl('^{-1}')+']'
    
    polyfill, [xran(0), xran(1), xran(1), xran(0), xran(0)], [yran(0), yran(0), res.background, res.background, yran(0)], color = 240
    
    plot, indgen(10), ps = -4,thick = 1, xtit = 'Array index', ytit = 'Flux scaling factor', $
      /nodata, /noerase, xran = xran, xsty = 1, yran = yran, ylog = ylog, ysty = 1+8
    axis, yaxis = 1, ysty = 1, yran = yran / (9.4*!pi/180./3600.)^2., color = 0, ytitle = '[W m'+textoidl('^{-2}')+' sr'+textoidl('^{-1}')+']'
    
    plots, !x.crange, [0, 0]
    loadct, 0, /sil
    plots, xran, [0, 0], line = 0, col = 50
    plots, [1, 1], yran, line = 0, col = 50
    plot, xarr, [obs+err, model*(fluxfinal+errorfinal) + res.background], ps = -4,thick = 1, xtit = 'Array index', ytit = 'Flux scaling factor', $
      /nodata, /noerase, xran = xran, xsty = 1, yran = yran, ylog = ylog, ysty = 1+8
    ;for i = 0, n_elements(model) - 2 do begin
    ;  x = 1+[i, i+1, i+1, i, i]
    ;  y_u = model * (fluxfinal[0] + errorfinal[0]) + res.background
    ;  y_d = model * (fluxfinal[0] - errorfinal[0]) + res.background
    ;  polyfill, x, [y_d(i), y_d(i+1), y_u(i+1), y_u(i), y_d(i)], color = 100, noclip = 0
    ;endfor
    loadct, 38, /sil
    oplot, xarr, model, color = 0, ps = 4, thick = 2, symsize = 2.5
    
    loadct, 38, /sil
    oploterror, xarr, obs, err, ps = 3,thick = 1, col = 200, errcol = 200
    if c_good gt 0 then for i = 0, n_elements(ind_good) do oploterror, xarr(ind_good), obs(ind_good), err(ind_good), ps = sym(4), color = 200, errcolor = 200, symsize = 2
    
    al_legend, ["F = " + a2str(res.flux) + ' +/- ' + a2str(res.error), "R_Chi2 = " + str_round(res.redchi2)], box = 0, /top, /right, charsize = 0.8
    
    al_legend, ["Offset = " + str_round(res.offset(0)) + ', ' + str_round(res.offset(1)), 'Source centroid: ' + a2str(res.position(0)) + ', ' + a2str(res.position(1)), "Offset = " + str_round(res.offset(2)) + ', ' + str_round(res.offset(3)), 'Source centroid: ' + a2str(res.position(2)) + ', ' + a2str(res.position(3))], box = 0, /bottom, /left, charsize = 0.8
    
    ;device, /close
    if keyword_set(cont) then filename = strtrans(file, '.sav', '_Opt_continuum.png') else filename = strtrans(file, '.sav', '_Opt.png')
    ;write_png, filename, tvrd(/true)
    
    device, /close
    set_plot, 'X'
    
  endif
  
  
end
