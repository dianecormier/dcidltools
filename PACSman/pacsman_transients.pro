
;################################################################################
function dipper_model, x, p
  return, p[0] - p[1] * exp(-p[2]*x) + p[3] * x
;continuum, depth, curvature
;           50., 1.
end
;################################################################################
function fader_model, x, p
  return, p[0] - ( p[1] * exp(-p[2]*x^p[3]) ) ;- p[4] * x
end
;################################################################################
function fader_model_discontinuity, x, p, center = center
  f = p[0]
  f += (x gt center+1) * ( p[1] * exp(-p[2]*x^p[3]) ) - p[4] * x
  return, f
end
;################################################################################
FUNCTION PACSman_BINARY, number, COLOR=color, SEPARATE=separate, NOREVERSE=noreverse

  ON_ERROR, 2
  
  ; What kind of number is this?
  thisType = SIZE(number, /Type)
  
  CASE thisType OF
  
    1: BEGIN ; Byte value
    
      bin = STRARR(8)
      FOR j=0,7 DO BEGIN
        powerOfTwo = 2L^j
        IF (LONG(number) AND powerOfTwo) EQ powerOfTwo THEN $
          bin(j) = '1' ELSE bin(j) = '0'
      ENDFOR
      IF Keyword_Set(color) THEN bin = [bin, STRARR(16)+'0']
    ENDCASE
    
    2: BEGIN ; Integer value.
    
      bin = STRARR(16)
      FOR j=0,15 DO BEGIN
        powerOfTwo = 2L^j
        IF (LONG(number) AND powerOfTwo) EQ powerOfTwo THEN $
          bin(j) = '1' ELSE bin(j) = '0'
      ENDFOR
      IF Keyword_Set(color) THEN bin = [bin, STRARR(8)+'0']
    ENDCASE
    
    3: BEGIN ; Long integer value.
    
      number = LONG(number)
      bin = STRARR(32)
      FOR j=0,31 DO BEGIN
        powerOfTwo = 2L^j
        IF (LONG(number) AND powerOfTwo) EQ powerOfTwo THEN $
          bin(j) = '1' ELSE bin(j) = '0'
      ENDFOR
      IF Keyword_Set(color) THEN bin = bin[0:23]
      
    ENDCASE
    
    ELSE: Message, 'Only BYTE, INTEGER, and LONG values allowed.'
    
  ENDCASE
  
  ; Do we need to separate in groups of 8?
  IF Keyword_Set(separate) THEN BEGIN
    CASE N_Elements(bin) OF
      8:
      16: bin = [bin[0:7], '  ', bin[8:15]]
      24: bin = [bin[0:7], '  ', bin[8:15], '  ', bin[16:23]]
      32: bin = [bin[0:7], '  ', bin[8:15], '  ', bin[16:23], '  ', bin[24:31]]
    ENDCASE
  ENDIF
  
  ; Reverse the array, so highest bits are on left and lowest bits are on right.
  if keyword_set(noreverse) then RETURN, bin else RETURN, Reverse(bin)
  
END
;################################################################################
;;converts binary mask into a boolean mask
;function convert2bool, maskin, readouts
;  maskout = long(intarr(readouts, 25, 18))
;  for xx = 0, 24 do for yy = 0, 17 do begin
;    arr = long(intarr(1))
;    for i = long(0), n_elements(maskin(*, xx, yy)) - 1 do arr = [arr, nint(PACSman_binary(maskin(i, xx, yy), /norev))]
;    ;maskout(x, y, *) = arr[1:n_elements(arr)-1]
;    ;help, arr, readouts, maskin
;    maskout(*, xx, yy) = arr[1:readouts]
;  endfor
;  arr = 0
;  return, maskout
;end
;______________________________________________________________________________
;converts binary mask into a boolean mask
function convert2bool, maskin
  nn = n_elements(maskin(0, 0, *)) * 32
  maskout = intarr(5, 5, nn)
  for x = 0, 4 do for y = 0, 4 do begin
    arr = intarr(1)
    for i = long(0), n_elements(maskin(x, y, *)) - 1 do arr = [arr, nint(PACSman_binary(maskin(x, y, i), /noreverse))] ;somehow the mask extensions are the transposed matrix (5x5)
    maskout(x, y, *) = arr[1:nn]
  endfor
  return, maskout
end
;################################################################################
;################################################################################
;################################################################################
;################################################################################
;################################################################################
pro PACSman_transients, frame = frame, dippers = dippers, verbose = verbose, save_plots = save_plots, faders = faders, denoise = denoise, $
  redomask = redomask, reflag = reflag;, brightline = brightline

  ;cd, current = rootdir
  ;rootdir += '/'

  if file_test('dir_out') then begin
    readcol, 'dir_out', dir, format = '(a)'
    dir = dir(0) + '/'
  endif else dir = './'
  if keyword_set(frame) then file = frame else file = dir + 'frame.fits'
  if file_test(file) eq 0 then begin
    print, 'File not found: ' + file
    retall
  endif
  
  ;###########################################
  ;READ CUBE'S EXTENSIONS
  if keyword_set(verbose) then begin
    print, '###########################################'
    print, ' Reading extensions...'
  endif
  !ERROR_STATE.CODE=0
  extensions = strarr(25)
  n = 0
  while !ERROR_STATE.CODE eq 0 do begin
    test = readfits(file, h, exten_no=n, /silent)
    test = 0
    extensions(n) = strlowcase(strtrim(sxpar(h, 'EXTNAME')))
    n += 1
  endwhile
  extensions = extensions(0:n-1)
  
  signal = readfits(file, h, exten_no=(where(extensions eq 'signal'))[0], /silent)
  ;;  noise = readfits(file, h, exten_no=(where(extensions eq 'noise'))[0], /silent)
  ind = where(extensions eq 'noise', c)
  if c gt 0 then noise = readfits(file, h, exten_no=(where(extensions eq 'noise'))[0], /silent) else noise = 0.1*signal
  
  ;###########################################
  ;FLAG PIXELS
  if keyword_set(verbose) then begin
    print, '###########################################'
    print, ' flagging bad pixels...'
  endif
  readouts = n_elements(signal(*, 0, 0))
  if file_test(strtrans(file, '.fits', '_flag.fits')) and not keyword_set(redomask) then mask_cube = readfits(strtrans(file, '.fits', '_flag.fits')) else begin
    mask_cube = 1 + intarr(readouts, 25, 18)
    mask_names = ['saturation', 'rawsaturation', 'noisypixels', 'badpixels', 'glitch', 'outliers', 'blindpixels'];, 'uncleanchop', 'gratmove'
    for mm = 0, n_elements(mask_names) - 1 do begin
      ind = where(extensions eq strtrim(mask_names(mm)), c)
      if c gt 0 then begin
        ind = ind[0]
        ;mask_i_cube = long(readfits(file, h, exten_no=ind, /silent))
        mask_i_cube = readfits(file, h, exten_no=ind, /silent) ;longint when reading
        ;help, mask_i_cube
        mask_i_cube = convert2bool(mask_i_cube);, readouts)
        ind2 = where( mask_i_cube, c )
        str = [strtrim(mask_names(mm))+': '+a2str(c)+'/'+a2str(readouts*25*18)]
        if keyword_set(verbose) then print, str
        if c gt 0 then mask_cube(ind2) = 0
      endif
    endfor
    mask_i_cube = 0
    writefits, strtrans(file, '\.fits', '_flag.fits'), mask_cube
  endelse
  
  spawn, '\rm -rf ' + dir + '/plots/'
  spawn, '\mkdir -p ' + dir + '/plots/'
  spawn, '\mkdir -p ' + dir + '/plots/dippers/'
  
  ;###########################################
  ;MAIN LOOP
  ;###########################################
  set_plot, 'PS'
  device, filename = dir+'/plots/faders.ps', /color
  
  
  if keyword_set(faders) then begin
  
    if keyword_set(verbose) then begin
      print, '###########################################'
      print, ' removing transients...'
      print, ' /!\ this is supposed to work only for sources with negligible continuum emission!!'
      print, ''
    endif
    
    for spaxel = 0, 24 do begin
    
      for spectral_pixel = 1, 16 do begin
      
        ;if (spaxel ne 12) or (spectral_pixel ne 8) then continue
        ;if (spaxel ne 15) or (spectral_pixel ne 5) then continue
      
        mask = reform(mask_cube(*, spaxel, spectral_pixel))
        sig = reform(signal(*, spaxel, spectral_pixel))
        ;sig_orig = sig
        ;;      err = reform(noise(*, spaxel, spectral_pixel))
        
        ind_bad = where( mask eq 0, c_bad, complement = ind_good )
        if n_elements(ind_good)/float(n_elements(sig)) lt 0.5 or n_elements(where(finite(sig) eq 1)) lt 10 then begin
          print, spaxel, spectral_pixel
          print, ' most of the pixels are bad! Skipping...'
          print, a2str(n_elements(ind_good))+'/'+a2str(n_elements(sig))
          print, ' spaxel : ', spaxel
          print, ' spectral pixel : ', spectral_pixel
          continue
        endif
        ind_good_orig = ind_good
        
        loadct, 0, /sil
        plot, ind_good, sig(ind_good), xtitle = 'Readouts', ytitle = 'Flux density [Jy]', /nodata, title = 'Signal vs. readouts, spaxel '+a2str(spaxel+1)+', spectralpixel '+a2str(spectral_pixel+1), ysty = 1
        if c_bad gt 0 then oplot, ind_bad, sig(ind_bad), color = 100, ps = 4
        oplot, ind_good, sig(ind_good)
        loadct, 38, /sil
        
        ;if ind_bad(0) ne -1 then oplot, ind_bad, sig(ind_bad), psym = 4, color = 200
        
        ;for k = long(0), n_elements(sig), 7500 do plots, [k+1500, k+1500], !y.crange, linestyle = 2, color = 100
        ;for k = long(0), 76499, 7500 do   plots, [k+4500, k+4500], !y.crange, linestyle = 2, color = 100
        ;for k = long(0), 76499, 7500 do   plots, [k+7500, k+7500], !y.crange, linestyle = 2, color = 100
        
        ;      arr = fltarr(10000)
        ;      k = 0
        ;      step = 100
        ;      for ii = long(0), n_elements(sig), step do begin
        ;        indd = where( ind_good gt ii and ind_good lt ii + step)
        ;        tmp = reform(sig(ind_good(indd)))
        ;        tmpsmooth = median(tmp, 25)
        ;        sdev = median( abs(tmp - tmpsmooth) )
        ;        ;m = moment(sig(ind_good(indd)), sdev = sdev)
        ;        ;print, sdev
        ;        plots, [ii+0.5*step, ii+0.5*step], [sdev, sdev]*median(sig(ind_good))/10., ps = 4, color = 100
        ;        arr(k) = sdev
        ;        k += 1
        ;      endfor
        ;      arr = arr(0:k-1)
        ;      plots, !x.crange, [median(arr), median(arr)]*median(sig(ind_good))/10., linestyle = 1
        
        
        nscales = 8
        
        ;;-- DC modified 2013-02-15; expand signal size to correct transient after calibration block
        sig2 = fltarr(readouts+400)
        sig2(0:199) = sig(0)
        sig2(200:readouts+199) = sig
        sig2(readouts+200:readouts+399) = sig(readouts-1)
        t2 = mmt_transform(sig2, nscales)
        t = t2(*,200:readouts+199)
        ;t = mmt_transform(sig, nscales)
        
        ind = where(t(nscales, *) eq sig, c, complement = noind)
        
        sig(noind) -= t(nscales, noind) ;transients
        
        ;print, t(nscales, *)
        if keyword_set(denoise) then sig -= t(0, *) + t(1, *) ;denoise
        
        arr = reform(t(nscales, *))
        arr = arr[2^(nscales-1)+10:n_elements(arr)-2^(nscales-1)-10-1] ;skip the edges
        
        ;  if (min(arr) ge 0.) then ind = where(arr lt 1.1*min(arr)) $
        ;else ind = where(arr lt 0.8*min(arr)) ;guesses the real level
        
        
          ;;-- DC modified 2012-11-06; min(arr) was negative on N44 frame
          if (min(arr, /nan) ge 0.) then indlvl = where(arr lt 1.1*min(arr, /nan)) $
            else indlvl = where(arr lt 0.8*min(arr, /nan)) ;guesses the real level
          sig += median(arr(indlvl)) ;re-shift up
        
        ;indlvl = where(arr lt 1.1*min(arr))
        ;sig(noind) += median(arr(indlvl)) ;re-shift up
        
        oplot, t(nscales, *), color = 200
        t = 0
        arr = 0
        
        ;arr = fltarr(10000)
        ;k = 0
        ;for ii = long(0), 75000, 400 do begin
        ;  indd = where( ind_good gt ii and ind_good lt ii + 400)
        ;  tmp = reform(sig(ind_good(indd)))
        ;  tmpsmooth = median(tmp, 25)
        ;  sdev = median( abs(tmp - tmpsmooth) )
        ;  ;m = moment(sig(ind_good(indd)), sdev = sdev)
        ;  ;print, sdev
        ;  plots, [ii+200, ii+200], [sdev, sdev]*median(sig(ind_good))/10., ps = 4
        ;  arr(k) = sdev
        ;  k += 1
        ;endfor
        ;arr = arr(0:k-1)
        ;plots, !x.crange, [median(arr), median(arr)]*median(sig(ind_good))/10., linestyle = 1
        
        loadct, 38, /silent
        oplot, ind_good, sig(ind_good), color = 150
        
        
        signal(*, spaxel, spectral_pixel) = temporary(sig)
        
      endfor
      
    endfor
    
  endif
  device, /close
  
  ;###########################################
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  ;###########################################
  
  if keyword_set(dippers) then begin
  
    ;###########################################
    ;FIND GLITCHES
    if keyword_set(verbose) then begin
      print, '###########################################'
      print, ' finding dippers...'
      print, ''
    endif
    for spaxel = 0, 24 do begin
    
      for spectral_pixel = 1, 16 do begin
      
      
        mask = reform(mask_cube(*, spaxel, spectral_pixel))
        sig = reform(signal(*, spaxel, spectral_pixel))
        
        if keyword_set(reflag) then begin
          ;reflag some outliers?
          nscales = 2
          t = mmt_transform(sig, nscales)
          sdev = abs( sig-t(nscales, *) )
          ind_bad = where( sdev gt 10.*median(sdev), c_bad )
          if c_bad gt 0 then mask(ind_bad) = 0
        endif
        
        sig_orig = sig
        ind_good = where( mask eq 1, complement = ind_bad )
        ind_good_orig = ind_good
        
        pw = 25 ;edges
        k = 0
        
        ;cycling through bad pixels only
        for ii = n_elements(ind_bad) - 1, 0, -1 do begin
        
          i = ind_bad(ii)
          if (i le pw) or (i ge n_elements(sig)-1-pw) then continue
          
          x = findgen(pw)
          sig2 = sig[i:i+pw-1]
          mask2 = mask[i:i+pw-1]
          ind_bad2 = where( mask2 eq 0, c_bad2, complement = ind_good2 )
          if c_bad2 eq 0 then stop ;that wouldn't be normal
          if n_elements(ind_good2) le 1 then continue
          
          ind = where( finite(sig2(ind_good2)), c)
          if c lt 10 then continue
          
          coeffs = bspline_iterfit(x(ind_good2), sig2(ind_good2), bkspace=15, yfit=yfit, max_iter=10, nord=3)
          sdev = median(abs(sig2(ind_good2) - temporary(yfit)))
          
          ;keep positive glitches only, remask all the <0 bad pixels
          for j = 0, n_elements(ind_bad2) - 1 do begin
            base = median(sig2(max([ind_bad2(j)-5, 0]):min([ind_bad2(j)+5, n_elements(sig2)-1])))
            if sig2(ind_bad2(j)) lt base then mask2(ind_bad2(j)) = 1
          endfor
          ;reindex glitches
          ind_bad2 = where( mask2 eq 0, c_bad2, complement = ind_good2 )
          
          sig3 = sig2(sort(sig2))
          cont = median(sig3[floor(pw/2.):pw-1])
          
          w = 0.*sig2 + 1.
          
          parinfo = replicate({relstep: 0., value:0.D, fixed:0, limited:[0, 0], limits:[0.D, 0]}, 4)
          parinfo(3).fixed = 1
          ;parinfo[3].limited[0] = 1
          ;parinfo[3].limits[0]  = 0.D
          parinfo[*].value = [cont, 50., 1., 0.]
          params = mpfitfun('dipper_model', x(ind_good2), sig2(ind_good2), weights = w(ind_good2), parinfo = parinfo, /nan, /quiet, yfit = yfit, STATUS=status, $
            ERRMSG=errmsg, dof = dof, bestnorm = bestnorm, perror = perror)
          if status LE 0 then message, errmsg
          chisq = sqrt(bestnorm/dof)
          
          sig2 = 0
          sig3 = 0
          mask2 = 0
          
          ;first test
          ;    return, p[0] - p[1] * exp(-p[2]*x) + p[3] * x
          ;  depth                        depth             x-scale                 x-scale                  slope                 chi2
          if (params[1] gt 10.*sdev) and (params[1] gt 0.) and (params[2] gt 0.05) and (params[2] lt 0.7) and (params(3) ge 0.) and (chisq lt 5.) then begin
          
            ;probably a dipper, continue with 2nd fit
            pw2 = 100
            
            sig2 = sig[i:min([i+pw2, n_elements(sig)-1])]
            x = findgen(n_elements(sig2))
            mask2 = mask[i:min([i+pw2, n_elements(sig)-1])]
            ind_good2 = where( mask2 eq 1 )
            ind_bad2 = where( mask2 eq 0 )
            ;keep negative glitches only
            if ind_bad2(0) ne -1 then begin
              for j = 0, n_elements(ind_bad2) - 1 do begin
                base = median(sig2(max([ind_bad2(j)-5, 0]):min([ind_bad2(j)+5, n_elements(sig2)-1])))
                if sig2(ind_bad2(j)) lt base then mask2(ind_bad2(j)) = 1
              endfor
              ind_bad2 = where( mask2 eq 0, complement = ind_good2 )
            endif
            
            sig3 = sig2(sort(sig2))
            cont = median(sig3[floor(n_elements(sig3)/1.5):n_elements(sig3)-1])
            
            w = 0.*sig2 + 1.
            
            parinfo = replicate({relstep: 0., value:0.D, fixed:0, limited:[0, 0], limits:[0.D, 0]}, 4)
            parinfo(3).fixed = 1
            ;parinfo[3].limited[0] = 1
            ;parinfo[3].limits[0]  = 0.D
            parinfo[2].limited[0] = 1
            parinfo[2].limits[0]  = 0.2D
            parinfo[2].limited[1] = 1
            parinfo[2].limits[1]  = 2.D
            parinfo[*].value = [cont, 50., 1., 0.]
            params = mpfitfun('dipper_model', x(ind_good2), sig2(ind_good2), weights = w(ind_good2), parinfo = parinfo, /nan, /quiet, yfit = yfit, STATUS=status, $
              ERRMSG=errmsg, dof = dof, bestnorm = bestnorm, perror = perror)
            if status LE 0 then message, errmsg
            chisq = sqrt(bestnorm/dof)
            
            ;2nd test
            if (params[1] lt 10.*sdev) or (params[1] lt 0.) or (params[2] lt 0.05) or (params[2] gt 0.7) or (params(3) lt 0.) or (chisq gt 5.)then continue
            
            if keyword_set(verbose) then print, 'Dipper found!!'
            k += 1
            
            if keyword_set(verbose) then print, i, params, chisq
            
            if keyword_set(save_plots) then begin
              if k eq 1 then device, filename = dir+'/plots/dippers/spaxel'+a2str(spaxel+1)+'_spectralpixel'+a2str(spectral_pixel+1)+'_dippersfits.ps', /color
              ;xarr = max([i-50, 0])+findgen(min([i+50+pw-1, n_elements(sig)-1])-max([i-50, 0])+1)
              loadct, 0, /sil
              plot, sig[max([i-50, 0]):min([i+pw2, n_elements(sig)-1])], xtitle = 'Readouts', ysty = 1, ytitle = 'Flux density [Jy]', $
                title = 'Dip #'+a2str(k)+', readout='+a2str(i)+', depth='+a2str(params(1)), subtitle = 'curve='+a2str(params(2))+', slope='+a2str(params(3))+$
                ', sdev='+a2str(sdev)+', chi='+a2str(chisq)
              loadct, 38, /sil
              if ind_bad2(0) ne -1 then oplot, x(ind_bad2)+50, sig2(ind_bad2), color = 100, psym = 4
              oplot, x+50, dipper_model(x, params), col = 50
            endif
            
            w1 = reverse(findgen(n_elements(sig2)))
            w2 = findgen(n_elements(sig2))
            ;sig[i:min([i+pw2, n_elements(sig)-1])] *= (w1*params(0)/dipper_model(x, params) + w2) / (w1+w2)
            sig[i:min([i+pw2, n_elements(sig)-1])] += (params(0)-dipper_model(x, params))
            
            ;sig[i:min([i+pw2, n_elements(sig)-1])] *= params(0)/glitch_model(x, params)
            sig[i-1:i] = params(0)
            
            sig2 = 0
            sig3 = 0
            mask2 = 0
            
            if keyword_set(save_plots) then begin
              oplot, sig[max([i-50, 0]):min([i+pw2, n_elements(sig)-1])], col = 100, thick = 2
            endif
            
          endif
          
        endfor
        
        if keyword_set(save_plots) and k gt 0 then begin
          loadct, 0, /silent
          ;device, filename = dir+'/plots/dippers/spaxel'+a2str(spaxel+1)+'_spectralpixel'+a2str(spectral_pixel+1)+'_minusdippers.ps', /color
          plot, ind_good, sig(ind_good), xtitle = 'Readouts', ytitle = 'Flux density [Jy]', /nodata, title = 'Signal vs. readouts (dippers removed)', ysty = 1
          oplot, ind_good_orig, sig_orig(ind_good_orig)
          if ind_bad(0) ne -1 then oplot, ind_bad, sig(ind_bad), psym = 4, color = 200
          loadct, 38, /silent
          
          oplot, ind_good, sig(ind_good), color = 100
          device, /close
        endif
        ;stop
        
        
        signal(*, spaxel, spectral_pixel) = temporary(sig)
        
        mask = 0
        x = 0
        w = 0
        w1 = 0
        w2 = 0
        
        writefits, dir + 'signal_corrected.fits', signal
        
      endfor ;spectral pixels
      
    endfor ;spaxels
  endif
  
end
