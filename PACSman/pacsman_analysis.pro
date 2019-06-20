;+
; NAME:
;   PACSman_analysis
;
; AUTHOR:
;   Vianney Lebouteiller, Diane Cormier, CEA/SAp, Saclay, France
;   vianney.lebouteiller@cea.fr
;
;   UPDATED VERSIONs can be found at:
;      http://www.myravian.fr/Homepage/Softwares.html
;
; PURPOSE:
;   Map analysis (cuts, flux extraction, geometry)
;
; MAJOR TOPICS:
;
;
; CALLING SEQUENCE:
;   PACSman_analysis
;
; DESCRIPTION:
;
; INPUTS:
;
;   none, has to be called from the directory where the line maps are located
;
; RETURNS:
;
;
; OPTIONAL PARAMETERS:
;
; REFERENCES:
;
; DEPENDENCIES:
;
; MPFIT : http://cow.physics.wisc.edu/~craigm/idl/idl.html
;
; IDL ASTROLIB (readfits, extast, hrot, hastrom, ...) :
; http://idlastro.gsfc.nasa.gov/contents.html
;
; MODIFICATION HISTORY:
;   Written, Aug 2010, VL
;
;
;##############################################################################
;
; LICENSE
;
;  Copyright (C) 2010 CEA/SAp
;
;  PACSman is distributed in the hope that it will be useful, but
;  WITHOUT ANY WARRANTY; without even the implied warranty of
;  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
;  General Public License for more details.
;
;##############################################################################

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


;##############################################################################
function nicecol, j, n
  res = 35/(n/7.)*j+45
  if res gt 250 then res = 0
  return, long(res)
end
;---------------------------------------------------------------------------------
FUNCTION str_round, input
  x = input
  ;stupid function to round the number in the string...
  if strpos(x, 'NaN') gt -1 or strpos(x, 'Inf') gt -1 then return, 'NaN' else begin
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
    rst = strmid(sr[1], 0, 2)
    tmp = sr[0] + '.' + rst
    if tmp eq '-0.00' then tmp = '0.00'
    return, tmp + expo
  endelse
END
;##############################################################################
;displays a nice WCS string from RA and DEC in degrees
FUNCTION radec2radec, ra, dec
  disp_ra  = ra
  disp_dec = dec
  
  neg_dec  = disp_dec LT 0
  radec, disp_ra, abs(disp_dec), ihr, imin, xsec, ideg, imn, xsc
  wcsstring = string(ihr, imin, xsec, ideg, imn, xsc, format = '(i2.2,"h",i2.2,"m",f6.3,"s   ",i2.2,":",i2.2,":",f5.2)' )
  if (strmid(wcsstring, 6, 1) EQ ' ') then $
    strput, wcsstring, '0', 6
  if (strmid(wcsstring, 22, 1) EQ ' ') then $
    strput, wcsstring, '0', 22
  IF neg_dec THEN strput, wcsstring, '-', 15
  
  return, wcsstring
END
;##############################################################################
;converts anything to string
function a2str,number,length=length,format=format,trail=trail,character=char
  on_error,2
  n = n_elements(number)
  s = strtrim(string(number,format=format),2)
  if 1-keyword_set(char) then char = '0'
  if n_elements(length) gt 0 then begin
    ilen = strlen(s)
    for i = 0,n-1 do begin
      nz = length-ilen[i]
      if nz gt 0 then begin
        for j=0,nz-1 do begin
          if keyword_set(trail) then s[i] = s[i]+char else $
            s[i] = char+s[i]
        endfor
      endif
    endfor
  endif
  return,s
end
;;##############################################################################
function gauss2d, nx, ny, r
  x = findgen(nx) # replicate(1.0, ny)  ;Create X and Y arrays
  y = replicate(1.0, nx) # findgen(ny)
  u = ((x-r(4))/r(2))^2. + ((y-r(5))/r(3))^2.  ;Create ellipse
  z = r(0) + r(1) * exp(-u/2.)   ;to gaussian
  return, z
end
;;##############################################################################
;performs 2D fit of the maps
pro gauss2dfit_image, s, isurf = isurf

  loadct, 39, /sil
  case s.radius_unit of
    'pixel': begin
      radius = s.radius
      radius_ann = s.radius_ann
    end
    'arcsec': begin
      radius = s.radius / (9.4/3.)
      radius_ann = s.radius_ann / (9.4/3.)
    end
  endcase
  
  imag = reform(s.result(*, *, s.image_index))
  left = max([0, s.image_size(0)-s.aper_center(0)-radius])
  right = min([s.image_size(0)-1, s.image_size(0)-s.aper_center(0)+radius])
  bottom = max([0, s.aper_center(1)-radius])
  top = min([s.image_size(1)-1, s.aper_center(1)+radius])
  
  subim = double(imag[left:right, bottom:top])
  err = reform(s.noise(*, *, s.image_index))
  suberr = double(err[left:right, bottom:top])
  ;res = gauss2dfit(subim, a, /tilt)
  
  parinfo = replicate({value:0.D, fixed:0, limited:[0, 0], limits:[0.D, 0]}, 7)
  parinfo[0].value = median(subim)
  parinfo[1].value = max(subim) - parinfo[0].value
  parinfo[2:3].value = replicate(0.3*radius, 2)
  parinfo[4:5].value = 0.5*size(subim, /dim)
  ind = where( finite(subim, /nan), c)
  if c gt 0 then subim(ind) = 0.
  ind = where( finite(suberr, /nan), c)
  if c gt 0 then suberr(ind) = max(suberr, /nan)
  res = mpfit2dpeak(subim, a, /tilt, parinfo = parinfo, error=suberr, nfree = nfree, bestnorm = bestnorm, perror = perror, status = status)
  
  pixel_sr = (9.4/3.)^2. * 2.3504e-11
  res *= pixel_sr
  
  extast, s.header, astr
  
  print, ''
  print, '2D gaussian fit'
  print, 'Fit status return: ', status
  
  dof = n_elements(subim) - n_elements(a)
  chisq = sqrt(bestnorm/dof)
  pcerror = 2. * perror * chisq ;1 sigma errors assuming dof is correct
  
  constant = a(0)
  scaling_factor = a(1)
  
  ;calculate error in flux
  f1 = total(gauss2d(right-left, top-bottom, a))
  a(0) -= pcerror(0)
  a(1) += pcerror(1)
  fm = total(gauss2d(right-left, top-bottom, a))
  a(0) += 2.*pcerror(0)
  a(1) -= 2.*pcerror(1)
  fp = total(gauss2d(right-left, top-bottom, a))
  ferr = 0.5*( (f1-fm) + (fp-f1) )
  ferr *= pixel_sr
  
  ferr /= float('1.'+strmid(total(res), strpos(total(res), "e-"), strlen(total(res))))
  
  ecc = abs(a(2)-a(3))/(a(2)+a(3))
  
  extast, s.header, astr
  xy2ad, s.aper_center(0), s.aper_center(1), astr, ra, dec
  
  ;if s.annulus_use eq 0 then constant_str = 'off' else constant_str = 'on'
  str = [ $
    'Line: ' + a2str(s.lines(s.image_index)), $
    ' ', $
    'constant        : '+a2str(constant)+' (W m-2 sr-1)', $
    'scaling factor  : '+a2str(scaling_factor)+' (W m-2 sr-1)', $
    'ellipse sigma   : '+a2str(a(2))+'px, '+a2str(a(3))+'px   ('+a2str(a(2)*9.4/3.)+'", '+a2str(a(3)*9.4/3.)+'")', $
    'eccentricity    : '+a2str(ecc), $
    'center          : '+a2str((a(4)+left))+'px, '+a2str(a(5)+bottom)+'px', $
    '                : '+a2str(ra)+'deg, '+a2str(dec)+'deg', $
    'angle           : '+a2str(a(6)*180./!pi)+' deg', $
    ' ', $
    ;'Integrated flux (constant subtraction: '+constant_str+')', $
    'Integrated flux', $
    '', $
    a2str(total(res))+' +/- '+a2str(ferr)+' W m-2', $
    'With sky subtraction', $
    a2str(total(res-constant*pixel_sr))+' +/- '+a2str(ferr)+' W m-2' $
    ]
  for i = 0, n_elements(str) - 1 do print, str(i)
  widget_control, s.shape_text, set_value=str
  
  openw, lun, 'ANALYSIS/Aperture/Fluxes.dat', /get_lun, width = 400, /append
  printf, lun, '____________________________________________________'
  printf, lun, 'Fluxes calculated from the 2D gaussian fit'
  for ll = 0, n_elements(str) - 1 do printf, lun, str(ll)
  printf, lun, '', ''
  close, /all
  
  xy2ad, s.image_size(0) - (a(4)+left), a(5)+bottom, astr, rafit, decfit
  
  wset, s.wset_image
  oplot, [rafit, rafit], [decfit, decfit], color = 10, psym = 4, symsize = 2
  fwhm = a(2:3)*2.*sqrt(2.*alog(2.))*9.4/3./3600.
  loadct, 0, /sil
  tvellipse, .5*fwhm(0)/cos(decfit*!pi/180.), .5*fwhm(1), rafit, decfit, a(6)*180./!pi, /data, color = 0, linestyle = 0, thick = 2
  sra = radius_ann * 9.4/3. / 3600.
  plots, [ra-sra/cos(decfit*!pi/180.), ra+sra/cos(decfit*!pi/180.)], [dec-sra, dec-sra], linestyle=0
  plots, [ra-sra/cos(decfit*!pi/180.), ra+sra/cos(decfit*!pi/180.)], [dec+sra, dec+sra], linestyle=0
  plots, [ra-sra/cos(decfit*!pi/180.), ra-sra/cos(decfit*!pi/180.)], [dec-sra, dec+sra], linestyle=0
  plots, [ra+sra/cos(decfit*!pi/180.), ra+sra/cos(decfit*!pi/180.)], [dec-sra, dec+sra], linestyle=0
  
  ;calculate PSF FWHM
  line = a2str(s.lines(s.image_index))
  i = 0
  test = 0
  while test eq 0 and i lt 100 do begin
    tmp = strmid(line, i, 1)
    test = max(tmp eq ['0', '1',  '2',  '3',  '4',  '5',  '6',  '7',  '8',  '9'])
    i += 1
  endwhile
  if i eq 100 then begin
    print, 'Line was not found, assuming the map is not that of a line. '
    lambda = 157.
  endif else lambda = strmid(line, i-1, strlen(line)-i+1)
  
  x = [62., 77., 90., 124., 154., 180.] ;lambda
  y = [9.2, 9.2, 9.2, 9.8, 11.3, 13.] ;fwhm in arcsec
  
  fwhm = interpol(y, x, lambda)
  
  xy2ad, 0, 0, astr, ra0, dec0
  ;for sig = 0., sigma, 0.01*sigma do tvellipse, sig/3600., sig/3600., ra0+sigma/3600., dec0+sigma/3600., 0., /data, color = 50
  fwhm /= 9.4/3.
  ;tvellipse, .5*fwhm/cos(decfit*!pi/180.), .5*fwhm, rafit, decfit, 0., /data, color = 150
  
  write_png, 'ANALYSIS/Cut/'+a2str(s.lines(s.image_index))+'_map.png', tvrd(/true)
  write_png, 'ANALYSIS/Aperture/'+a2str(s.lines(s.image_index))+'_map.png', tvrd(/true)
  
  loadct, 39, /sil
  
  wset, s.wset_show3_fit
  show3, res ;, ax=45,az=20;, shade=bytscl(res)
  if keyword_set(isurf) eq 1 then isurface, res
  write_png, 'ANALYSIS/Aperture/'+a2str(s.lines(s.image_index))+'_surfacefit.png', tvrd(/true)
  
  wset, s.wset_show3_data
  show3, subim;, ax=45,az=20;, shade=bytscl(subim)
  if keyword_set(isurf) eq 1 then isurface, subim
  write_png, 'ANALYSIS/Aperture/'+a2str(s.lines(s.image_index))+'_surface.png', tvrd(/true)
  
  
end

;##############################################################################
;update plots on the right panel
pro update_plots, s, xarray = xarray, yarray = yarray, farray = farray

  case s.radius_unit of
    'pixel': begin
      radius = s.radius
      radius_ann = s.radius_ann
    end
    'arcsec': begin
      radius = s.radius / (9.4/3.)
      radius_ann = s.radius_ann / (9.4/3.)
    end
  endcase
  
  if radius gt 0. then begin
  
    spawn, '\rm -f ANALYSIS/Cut/*.{dat,png}'
    spawn, '\rm -f ANALYSIS/Aperture/*.{dat,png}'
    
    tmp = readfits(s.files(0), /silent)
    nlines = n_elements(s.lines)
    
    ;horizontal line
    x0 = s.line_cut_start(0)
    y0 = s.line_cut_start(1)+1
    x1 = s.line_cut_end(0)
    y1 = s.line_cut_end(1)
    
    tmpx = [x0, x1]
    ;ind = sort(tmpx)
    ;tmpx = tmpx(ind)
    x0 = tmpx(0)
    x1 = tmpx(1)
    tmpy = [y0, y1]
    ;tmpy = tmpy(ind)
    y0 = tmpy(0)
    y1 = tmpy(1)
    
    ;center of image
    xc = s.image_size(0)-s.aper_center(0)
    yc = s.aper_center(1)
    
    
    ;------------------------------
    loadct, 38, /sil
    
    ;from http://www.dfanning.com/ip_tips/image_profile.html
    nPoints = ABS(x1-x0+1) > ABS(y1-y0+1)
    farray = fltarr(nPoints, nlines)
    earray = fltarr(nPoints, nlines)
    xloc = s.image_size(0) - 1 - (x0 + (x1 - x0) * Findgen(nPoints) / (nPoints - 1))
    yloc = y0 + (y1 - y0) * Findgen(nPoints) / (nPoints - 1)
    extast, s.header, astr
    
    openw, lun, /get_lun, 'ANALYSIS/Cut/linecut.dat', width = 800
    str = string('#X', format='(a-15)') + string('Y', format='(a-15)') + string('RA', format='(a-15)') + string('DEC', format='(a-15)')
    for line = 0, nlines - 1 do str += string('F('+s.lines(line)+')', format='(a-15)') + string('E('+s.lines(line)+')', format='(a-15)')
    printf, lun, str
    
    for j = 0, n_elements(xloc) - 1 do begin
      str = ''
      xy2ad, s.image_size(0)-xloc(j), yloc(j), astr, ra, dec
      for line = 0, nlines - 1 do begin
        tmp = reform(s.result[*, *, line])
        profile_f = interpolate(tmp, xloc, yloc)
        tmp = reform(s.noise[*, *, line])
        profile_e = interpolate(tmp, xloc, yloc)
        farray(*, line) = profile_f ;/ max(profile_f, /nan)
        earray(*, line) = profile_e/3. ;/ max(profile_f, /nan)
        str += string(a2str(profile_f(j)), format='(a-15)') + string(a2str(profile_e(j)), format='(a-15)')
      endfor
      str = string(a2str(xloc(j)), format='(a-15)') + string(a2str(yloc(j)), format='(a-15)') + string(a2str(ra), format='(a-15)') + string(a2str(dec), format='(a-15)') + str
      printf, lun, str
    endfor
    close, /all
    
    wset, s.wset_line_cut
    !p.position = [.1, .1, .85, .95]
    if max(farray) ne 0. then begin
      ind = where( farray gt 0. and earray gt 0. and farray gt 1.*earray and finite(farray) and finite(earray), c )
      if c gt 0 then begin
        loadct, 0, /sil
        plot, farray, title="Line Cut", ytitle="LOG Normalized Flux [W/m2/sr]", thick=2, yrange=minmax(farray(ind), /nan), $
          /nodata, /ylog, xsty=1, xran=[0, nPoints], xtitle = 'pixel', ysty = 1, background = 255, color = 0
        plots, replicate(0., 2), minmax(farray(ind), /nan), psym = 2, color = 0, symsize = 2
        plots, replicate(.5*nPoints, 2), minmax(farray(ind), /nan), psym = 4, color = 0, symsize = 2
        plots, replicate(nPoints, 2), minmax(farray(ind), /nan), psym = 4, color = 0, symsize = 2
        loadct, 38, /sil
        for j=0,nlines-1 do begin
          xarr = findgen(n_elements(farray(*, j)))
          ind = where( farray(*, j) gt 0. and earray(*, j) gt 0. and farray(*, j) gt 1.*earray(*, j) and finite(farray(*, j)) and finite(earray(*, j)), c )
          if c gt 0 then begin
            oploterror, xarr(ind), farray(ind, j), earray(ind, j), color=nicecol(j, nlines), $
              ps = -4, errcol = nicecol(j, nlines);, /nohat
            xyouts, nPoints, farray(max(ind), j), s.lines(j), charsize = 2, col = nicecol(j, nlines)
          endif
        endfor
      endif
    endif
    ;colarr = fltarr(nlines)
    ;for i = 0, nlines - 1 do colarr(i) = nicecol(i, nlines)
    ;al_legend,s.lines, linestyle=0, box=0, color=colarr, thick=2, /bottom, /right, textcol = 0;, position=[s.radius-2.5,0.1+0.05*j]
    
    write_png, 'ANALYSIS/Cut/linecut.png', tvrd(/true)
    
    ;FLUX GROWTH
    bkg=[-0d0,0d0]
    tmp_result = s.result
    tmp_noise = s.noise
    testarray = fltarr(nlines, 10000, 2) ;nlines, radius iterations, radius/flux
    sumarray = fltarr(nlines, 10000, 2) ;nlines, radius iterations, radius/flux
    sumerr = fltarr(nlines, 10000, 2)
    radius_step = .1
    for line=0, nlines - 1 do begin
      k = 0
      flux_in = reform(tmp_result(*, *, line))
      test_in = 0.*flux_in + 100.
      ind = where( finite(flux_in, /nan), c, complement = noc )
      if c gt 0 then begin  
        flux_in(ind) = 0.
        test_in /= n_elements(noc)
        test_in(ind) = 0.
      endif
      noise_in = reform(tmp_noise(*,*,line))^2.
      ind = where( finite(noise_in, /nan), c )
      if c gt 0 then noise_in(ind) = max(noise_in, /nan)
      for radi = radius, 1., 0.-radius_step do begin
        aper,flux_in,xc,yc, flux_in_radi, errap, sky, skyerr,1000,radi,/flux,/exact,setsky=0d0,/NaN,/silent
        sumarray(line, k, 0:1) = [0.1*nint(10.*radi), flux_in_radi]
        aper,test_in,xc,yc, flux_in_radi, errap, sky, skyerr,1000,radi,/flux,/exact,setsky=0d0,/NaN,/silent
        testarray(line, k, 0:1) = [0.1*nint(10.*radi), flux_in_radi]
        aper,noise_in,xc,yc, error_in_radi, errap, sky, skyerr,1000,radi,/flux,/exact,setsky=0d0,/NaN,/silent
        sumerr(line, k, 0:1) = [0.1*nint(10.*radi), error_in_radi]
        k += 1
      endfor
    endfor
    sumarray = sumarray[*, 0:k-1, 0:1]
    testarray = testarray[*, 0:k-1, 0:1]
    sumerr = sqrt(sumerr[*, 0:k-1, 0:1])
    
    ;FLUX
    str = strarr(nlines)
    str_ann = strarr(nlines)      
    pixel_sr = (9.4/3.)^2. * 2.3504e-11
    for line=0, nlines - 1 do begin
    
      sumarray(line, 0, 1) *= pixel_sr
      sumerr(line, 0, 1) *= pixel_sr
      ferr = sumerr(line, 0, 1) / float('1.'+strmid(sumarray(line, 0, 1), strpos(sumarray(line, 0, 1), "e-"), strlen(sumarray(line, 0, 1))))
      
      prestr = ''
      if line eq s.image_index then prestr = '---> '
      
      ;;NO SKY ANNULUS
      ;if s.annulus_use eq 0 then begin
      
      str(line) = prestr + s.lines(line)+': '+a2str(sumarray(line, 0, 1))+' +/- '+a2str(ferr)+' W m-2' + $
        ' ('+str_round(100.*sumarray(line, 0, 1)/(total(tmp_result(*,*,line), /nan)*pixel_sr))+'% of total=' + a2str(total(tmp_result(*,*,line), /nan)*pixel_sr) + ', '+str_round(testarray(line, 0, 1))+'% area)'
      
      ;endif else begin
      
      flux_in = reform(tmp_result(*,*,line))
      ind = where( finite(flux_in, /nan), c )
      if c gt 0 then flux_in(ind) = 0.
      
      ;noise_in = reform(tmp_noise(*,*,line))
      ;ind = where( finite(noise_in, /nan), c )
      ;if c gt 0 then noise_in(ind) = max(noise_in, /nan)
      
      ;SKY ANNULUS
      aper, flux_in, xc, yc, flux, errap, sky, skyerr, 1000, radius, [radius, radius_ann], /flux, /exact,/NaN,/silent
      ;;aper,sqrt(tmp_noise(*,*,line)),xc,yc, error, errap, sky, skyerr,1000,radius /flux,/exact,/NaN,/silent
      ;;ferr = error / float('1.'+strmid(flux, strpos(flux, "e-"), strlen(flux)))
      flux *= pixel_sr
      str_ann(line) = prestr + s.lines(line)+': '+a2str(flux)+' +/- '+a2str(ferr)+' W m-2'
      
    ;endelse
      
    endfor
    
    extast, s.header, astr
    xy2ad, s.aper_center(0), s.aper_center(1), astr, ra, dec
    
    str = ['Aperture properties', $
      'Radius: '+a2str(radius)+'px, ('+a2str(radius*9.4/3.)+'")', $
      'Annulus: '+a2str(radius_ann)+'px, ('+a2str(radius_ann*9.4/3.)+'")', $
      'Center: '+a2str(s.image_size(0)-s.aper_center(0))+'px, '+a2str(s.aper_center(1))+'px' , $
      '        '+a2str(ra)+'deg, '+a2str(dec)+'deg' , $
      '        '+radec2radec(ra, dec), $
      ' ', $
      'Integrated flux', $
      ' ', $
      str, $
      ' ', $
      'Integrated flux (sky subtraction)', $
      ' ', $
      str_ann ]
    for i = 0, n_elements(str) - 1 do print, str(i)
    widget_control, s.flux_text, set_value=str
    
    openw, lun, 'ANALYSIS/Aperture/Fluxes.dat', /get_lun, width = 400
    printf, lun, '____________________________________________________'
    printf, lun, 'Integrated fluxes calculated from the map data'
    for ll = 0, n_elements(str) - 1 do printf, lun, str(ll)
    printf, lun, '', ''
    close, /all
    
    wset, s.wset_sn_growth
    loadct, 0, /sil
    minf = 10.
    maxf = 0.
    nn = n_elements(sumarray(0, *, 0))
    for j = 0, nlines - 1 do begin
      for r = 0, n_elements(sumarray(0, *, 0)) - 1 do begin
        tmp = total(sumarray[j, 0:r, 1], /nan) / total(sumerr[j, 0:r, 1], /nan)
        tmp /= total(sumarray[j, *, 1], /nan) / total(sumerr[j, *, 1], /nan)
        minf = min([minf, tmp], /nan)
        maxf = max([maxf, tmp], /nan)
      endfor
    endfor
    plot, indgen(10), title="S/N Growth", xtitle="Radius (pix)", $
      ytitle="(S/N) / (S/N)_max", thick=2, xrange=[1.,radius_ann], $
      yrange=[minf, maxf], ysty=1, /nodata
    plots, [radius, radius], !y.crange, linestyle = 0, color = 100
    loadct, 39, /silent
    for j = 0, nlines - 1 do begin
      tmp = fltarr(nn)
      for r = 0, nn - 1 do begin
        tmp(r) = total(sumarray[j, 0:r, 1], /nan) / total(sumerr[j, 0:r, 1], /nan)
        tmp(r) /= total(sumarray[j, *, 1], /nan) / total(sumerr[j, *, 1], /nan)
      endfor
      oplot, sumarray[j, *, 0], tmp, color=nicecol(j, nlines)
    endfor
    
    write_png, 'ANALYSIS/Aperture/sngrowth.png', tvrd(/true)
    
    ;for line=0, nlines - 1 do begin
    ;  sumarray(line, *, 1) /= max(sumarray(line, *, 1), /nan) ;normalize with max
    ;  sumerr(line, *, 1) /= sumarray(line, 0, 1) * 1.e-18
    ;endfor
    
    wset, s.wset_flux_growth
    loadct, 0, /sil
    minf = 10.
    maxf = 0.
    nn = n_elements(sumarray(0, *, 0))
    for j = 0, nlines-1 do begin
      fluxes = (reform(sumarray[j, *, 1]) / sumarray(j, nn-1, 1)) / reform(sumarray[j, *, 0])^2.
      minf = min([minf, min(fluxes, /nan)], /nan)
      maxf = max([maxf, max(fluxes, /nan)], /nan)
    endfor
    plot, indgen(10), title="Flux Growth", xtitle="Radius (pix)", $
      ytitle="Normalized Flux / r^2", thick=2, xrange=[1., radius_ann], yrange=[minf, 2.], /nodata
    plots, [radius, radius], !y.crange, linestyle = 0, color = 100
    loadct, 39, /silent
    for j = 0, nlines-1 do begin
      fluxes = (reform(sumarray[j, *, 1]) / sumarray(j, nn-1, 1)) / reform(sumarray[j, *, 0])^2.
      oplot, reform(sumarray[j, *, 0]), fluxes, color=nicecol(j, nlines)
    endfor
    colarr = intarr(nlines)
    for i = 0, nlines - 1 do colarr(i) = nicecol(i, nlines)
    al_legend,s.lines, linestyle=0, box=0, color=colarr, thick=2, /top, /right;, position=[s.radius-2.5,0.1+0.05*j]
    
    write_png, 'ANALYSIS/Aperture/fluxgrowth.png', tvrd(/true)
    
    if xregistered('explodedview') then explode, s
    
  endif
  
end
;##############################################################################
pro explodedview_event, ev
  widget_control, ev.id, get_uvalue = s
;widget_control, s.exploded_draw, draw_xsize = ev.x
;widget_control, s.exploded_draw, draw_ysize = ev.y
end
;##############################################################################
;update plots on the right panel
pro explode, s, xarray = xarray, yarray = yarray, farray = farray

  case s.radius_unit of
    'pixel': begin
      radius = s.radius
      radius_ann = s.radius_ann
    end
    'arcsec': begin
      radius = s.radius / (9.4/3.)
      radius_ann = s.radius_ann / (9.4/3.)
    end
  endcase
  
  sc_size = get_screen_size(resolution=sc_res)
  
  nimages = float(n_elements(s.result(0, 0, *)))
  xxsize = 300.*nimages
  set_plot, 'X'
  if not xregistered('explodedview') then begin
    explodedview = widget_base(title=' ', /column, scr_xsize=sc_size(0), scr_ysize=220+xxsize/nimages, xsize=xxsize, ysize=xxsize/nimages, /scroll)
    imagebase = widget_base(explodedview, scr_xsize=sc_size(0), scr_ysize=30+xxsize/nimages, xsize=xxsize, /scroll)
    image = widget_draw(imagebase, xsize=xxsize, ysize=xxsize/nimages)
    exploded_text = widget_text(explodedview, xsize = 200, ysize = 10, /scroll)
    widget_control, explodedview, /realize
    widget_control, image, get_value=exploded_index ;for wset
    xmanager, 'explodedview', explodedview, /no_block
    s.exploded_index = exploded_index
    s.exploded_draw = image
    s.exploded_text = exploded_text
    widget_control, explodedview, set_uvalue = s
  endif
  wset, s.exploded_index
  
  x0 = s.line_cut_start(0)
  y0 = s.line_cut_start(1)
  x1 = s.line_cut_end(0)
  y1 = s.line_cut_end(1)
  
  extast, s.header, astr
  xy2ad, x0, y0, astr, ra0, dec0
  xy2ad, x1, y1, astr, ra1, dec1
  xline = [ra0, ra1]
  yline = [dec0, dec1]
  
  xy2ad, s.aper_center(0), s.aper_center(1), astr, rac, decc
  
  xy2ad, s.aper_center(0)+radius+1, s.aper_center(1), astr, rad, decd
  
  t = FINDGEN(361)/360.0 * 2.0 * 3.14
  xcirc = sin(t)
  ycirc = cos(t)
  R = (radius+1)*astr.cdelt[0]
  Rr = abs(rad-rac) ;* cos(avg(s.dec)*!pi/180.)
  
  !x.style=1
  !y.style=1
  
  for ii = 0, nimages - 1 do begin
  
    loadct, 39, /silent
    
    imag = s.result(*, *, ii)
    imag = sqrt((imag-min(imag, /nan))/max(imag, /nan))
    
    pacsman_disp, imag, xrange = [s.ra(0), s.ra(1)], xtitle = 'RA/DEC J2000 (deg)', $
      yrange = [s.dec(0), s.dec(1)], $
      aspect = 1, position = [ii/nimages, 0.1,  (ii+1.)/nimages, 1.], /noerase, /data;, xsty = 4, ysty = 4
      
    xyouts, s.ra(0), s.dec(0), s.lines(ii), /data, charsize = 2
    
    loadct, 0, /sil
    tvellipse, radius*astr.cdelt(0)/cos(decc*!pi/180.), radius*astr.cdelt(0), rac, decc, 0., /data, color = 255
    tvellipse, (radius_ann)*astr.cdelt(0)/cos(decc*!pi/180.), (radius_ann)*astr.cdelt(0), rac, decc, 0., /data, color = 255
    
    oplot, xline, yline,thick=2,linestyle=2, color = 255
    plots, xline(0), yline(0), ps = 2, color = 255
    plots, xline(1), yline(1), ps = 4, color = 255
    ;oplot, xline, yline, linestyle=0, psym=2, color= 255
    ;oplot, Rr*xcirc + rac, R*ycirc + decc,color=255,thick=2
    ;oplot, [rac, rac], [decc,decc], linestyle=0, psym=1,color=0
    
    
    arr = fltarr(3, s.image_size(0)*s.image_size(1))
    k = 0
    for x = 0, s.image_size(0) - 1 do for y = 0, s.image_size(1) - 1 do begin
      dist = sqrt( (s.image_size(0)-s.aper_center(0)-x)^2. + (s.aper_center(1)-y)^2. )
      if dist le s.radius then begin
        arr[*, k] = [x, y, s.result(x, y, ii)]
        k += 1
      endif
    endfor
    arr = arr[*, 0:k]
    ind = where(arr[2, *] eq max(arr[2, *], /nan))
    xy2ad, s.image_size(0)-(arr[0, ind]-.5), arr[1, ind]+.5, astr, ramax, decmax
    oplot, ramax, decmax, color = 0, psym = 7
    
  endfor
  
  write_png, 'ANALYSIS/Cut/exploded_view.png', tvrd(/true)
  write_png, 'ANALYSIS/Aperture/exploded_view.png', tvrd(/true)
  
;write_png, strtrans(s.files(s.image_index), '.fits', '_map.png'), tvrd(/true)
  
  
end
;##############################################################################
;updates the image and overplots
pro plotimage, s, objects = objects

  case s.radius_unit of
    'pixel': begin
      radius = s.radius
      radius_ann = s.radius_ann
    end
    'arcsec': begin
      radius = s.radius / (9.4/3.)
      radius_ann = s.radius_ann / (9.4/3.)
    end
  endcase
  
  wset, s.wset_image
  
  x0 = s.line_cut_start(0)
  y0 = s.line_cut_start(1)
  x1 = s.line_cut_end(0)
  y1 = s.line_cut_end(1)
  
  extast, s.header, astr
  xy2ad, x0, y0, astr, ra0, dec0
  xy2ad, x1, y1, astr, ra1, dec1
  xline = [ra0, ra1]
  yline = [dec0, dec1]
  
  ;if keyword_set(medline) then begin
  ;   xy2ad, p0(*,0), p0(*,1), astr, rap0,decp0
  ;   xy2ad, p1(*,0), p1(*,1), astr, rap1,decp1
  ;endif
  
  xy2ad, s.aper_center(0), s.aper_center(1), astr, rac, decc
  
  print, 'image size (px)        :', s.image_size(0), s.image_size(1)
  print, 'new center (px)        : ', s.aper_center(0), s.aper_center(1)
  print, 'RA/DEC of center (deg) : ', rac, decc
  
  
  xy2ad, s.aper_center(0)+radius+1, s.aper_center(1), astr, rad, decd
  
  t = FINDGEN(361)/360.0 * 2.0 * 3.14
  xcirc = sin(t)
  ycirc = cos(t)
  R = (radius+1)*astr.cdelt[0]
  Rr = abs(rad-rac) ;* cos(avg(s.dec)*!pi/180.)
  
  loadct, 39, /silent
  
  
  !x.style=1
  !y.style=1
  
  
  imag = s.result(*, *, s.image_index)
  imag = sqrt((imag-min(imag, /nan))/max(imag, /nan))
  
  pacsman_disp, imag, xrange = [s.ra(0), s.ra(1)], xtitle = 'RA/DEC J2000 (deg)', $
    yrange = [s.dec(0), s.dec(1)], $
    aspect = 1, position = [0.1, 0.1, 1., 1.]
    
  xyouts, s.ra(0), s.dec(0), s.lines(s.image_index), /data, charsize = 2
  
  
  loadct, 0, /sil
  tvellipse, radius*astr.cdelt(0)/cos(decc*!pi/180.), radius*astr.cdelt(0), rac, decc, 0., /data, color = 255
  tvellipse, (radius_ann)*astr.cdelt(0)/cos(decc*!pi/180.), (radius_ann)*astr.cdelt(0), rac, decc, 0., /data, color = 255
  
  oplot, xline, yline,thick=2,linestyle=2, color = 255
  plots, xline(0), yline(0), ps = 2, color = 255
  plots, .5*(xline(1)+xline(0)), .5*(yline(1)+yline(0)), ps = 4, color = 255
  plots, xline(1), yline(1), ps = 4, color = 255
  ;oplot, xline, yline, linestyle=0, psym=2, color= 255
  ;oplot, Rr*xcirc + rac, R*ycirc + decc,color=255,thick=2
  ;oplot, [rac, rac], [decc,decc], linestyle=0, psym=1,color=0
  
  
  arr = fltarr(3, s.image_size(0)*s.image_size(1))
  k = 0
  for x = 0, s.image_size(0) - 1 do for y = 0, s.image_size(1) - 1 do begin
    ;dist = sqrt( (s.image_size(0)-s.aper_center(0)-x)^2. + (s.aper_center(1)-y)^2. )
    dist = sqrt( (s.aper_center(0)-x)^2. + (s.aper_center(1)-y)^2. )
    ;if dist le s.radius then begin
    if dist le radius then begin
      arr[*, k] = [x, y, s.result(x, y, s.image_index)]
      k += 1
    endif
  endfor
  arr = arr[*, 0:k]
  ind = where(arr[2, *] eq max(arr[2, *], /nan))
  xy2ad, s.image_size(0)-(arr[0, ind]-.5), arr[1, ind]+.5, astr, ramax, decmax
  oplot, ramax, decmax, color = 0, psym = 7
  
  if keyword_set(objects) then begin
    case s.radius_unit of
      'pixel': radius = s.radius * (9.4/3.)
      'arcsec': radius = s.radius
    endcase
    extast, s.header, astr
    xy2ad, s.aper_center(0), s.aper_center(1), astr, rac, decc
    sss = 'http://simbad.harvard.edu/simbad/sim-script?submit=submit+script&script=output+script=off+console=off%0D%0A'
    sss = sss + 'set+epoch+2000%0D%0A'
    sss = sss + 'format%20object%20%22%25IDLIST%281%29%20%7c%20%25COO%28d;A%20D%29%22%0D%0A'
    sss = sss + 'query+coo+' + a2str(rac) + '+' + a2str(decc) + '+radius=' + a2str(radius/60.) + '&Radius.unit=arcmin%0D'
    query = webget(sss)
    if strpos(query.text[0], 'error') gt -1 then begin
      wmessage, text='No object found', title='Objects within aperture'
    endif else begin
      arr = query.text(0:n_elements(query.text)-2)
      ind = where(arr ne '')
      arr = arr(ind)
      wset, s.wset_image
      for i = 0, n_elements(arr) - 1 do begin
        arr(i) = a2str(i+1) + ' :' + query.text(i)
        tmp = (strsplit(query.text(i), '|', /extract))(1)
        tmp = strsplit(tmp, ' ', /extract)
        ;plots, tmp(0), tmp(1), psym=4
        xyouts, tmp(0), tmp(1), '#'+a2str(i+1)
      endfor
      wmessage, text=arr, title='Objects within aperture', ysize = min([20, n_elements(arr)]), xsize = 70
    endelse
  endif
  
;write_png, strtrans(s.files(s.image_index), '.fits', '_map.png'), tvrd(/true)
  
end

;##############################################################################
;calculates the position of the cursor in data coordinates
FUNCTION wid2pl, pos, s, x = x, y = y
  if keyword_set(x) then begin
    val = max([(pos-0.1*s.drawsize(0))*s.image_size(0)/(0.9*s.drawsize(0)), 0]) ;0.15 is the x_position given to pacsman_disp by sm_adopt_extract, 270 is the size of the image widget
    val = min([val, s.image_size(0)])
  endif
  if keyword_set(y) then begin
    val = max([(pos-0.1*s.drawsize(1))*s.image_size(1)/(0.9*s.drawsize(1)), 0]) ;0.15 is the x_position given to pacsman_disp by sm_adopt_extract, 270 is the size of the image widget
    val = min([val, s.image_size(1)])
  endif
  return, val
END

;##############################################################################
;events dispatcher
PRO PACSman_analysis_event, ev
  ;event handler
  widget_control, ev.id, get_uvalue=uvalue
  widget_control, ev.top, get_uvalue = s
  CASE uvalue OF
    'none':
    'tab_event':
    'explode': begin
      explode, s
    end
    'line_map': begin
      wset, s.wset_image
      s.image_index = ev.index
      update_plots, s
      plotimage, s
      gauss2dfit_image, s
    end
    'isurf': begin
      if ev.press eq 1 then begin
        update_plots, s
        plotimage, s
        gauss2dfit_image, s, isurf = 1
      endif
    end
    'image': begin
      ;print, ev.x, ev.y, s.fize
      x = wid2pl(ev.x, s, /x)
      y = wid2pl(ev.y, s, /y)
      widget_control, s.pos_xy_text, set_value = '('+str_round(a2str(x))+','+str_round(a2str(y))+')'
      widget_control, s.flux_text_otf, set_value= str_round(a2str(s.result(x, y, s.image_index)))
      if xregistered('explodedview') then begin
        str = ''
        for i = 0, n_elements(s.lines) - 1 do str += string(s.lines(i) + ': ' + str_round(a2str(s.result(x, y, i))), format='(a-28)')
        strline = ''
        for i = 0, n_elements(s.lines)*28 do strline = strline + '_'
        str = [str, strline]
        for i = 0, n_elements(s.lines) - 1 do begin
          str2 = ''
          for j = 0, n_elements(s.lines) - 1 do begin
            if i eq j or s.result(x, y, i) eq 0. then str2 += string('...', format='(a-28)') else $
              str2 += string(s.lines(j)+'/'+s.lines(i) + '=' + string(s.result(x, y, j)/s.result(x, y, i), format='(e-12.3)'), format='(a-28)')
          endfor
          str = [str, str2]
        endfor
        widget_control, s.exploded_text, set_value = str
      endif ;if xregistered('explodedview') then begin
      extast, s.header, astr
      xy2ad, s.image_size(0) - x, y, astr, ra, dec
      widget_control, s.pos_radec1_text, set_value = '('+str_round(a2str(ra))+', '+str_round(a2str(dec))+')'
      widget_control, s.pos_radec2_text, set_value = '('+radec2radec(ra, dec)+')'
      if ev.press eq 1 then begin
        s.aper_center(0) = s.image_size(0) - x
        s.aper_center(1) = y
        update_plots, s
        plotimage, s
        gauss2dfit_image, s
      endif
      if ev.press eq 4 then begin
        s.line_cut_start = [s.image_size(0) - x, y]
      endif
      if ev.release eq 4 then begin
        s.line_cut_end = [s.image_size(0) - x, y]
        update_plots, s
        plotimage, s
        gauss2dfit_image, s
        widget_control, s.wtab, set_tab_current = 2
      endif
    end
    'radius_change': begin
      widget_control, s.radius_btn, get_value = tmp
      s.radius = max([tmp, 0.])
      case s.radius_unit of
        'pixel': s.radius_ann = s.radius + 1.5
        'arcsec': s.radius_ann = s.radius + 1.5 * (9.4/3.)
      endcase      
      widget_control, s.radius_btn_ann, set_value = a2str(s.radius_ann)
      update_plots, s
      plotimage, s
      gauss2dfit_image, s
    end
    'radius_ann_change': begin
      widget_control, s.radius_btn_ann, get_value = tmp
      s.radius_ann = max([tmp, 0.])
      update_plots, s
      plotimage, s
      gauss2dfit_image, s
    end
    'radius_unit_pixel': begin
      s.radius_unit = 'pixel'
      update_plots, s
      plotimage, s
      gauss2dfit_image, s
    end
    'radius_unit_arcsec': begin
      s.radius_unit = 'arcsec'
      update_plots, s
      plotimage, s
      gauss2dfit_image, s
    end
    'get_objects': begin
      plotimage, s, /objects
    end
    'radius_left': begin
      s.radius -= 1.
      s.radius = max([s.radius, 0.])
      widget_control, s.radius_btn, set_value = a2str(s.radius)
      case s.radius_unit of
        'pixel': s.radius_ann = s.radius + 1.5
        'arcsec': s.radius_ann = s.radius + 1.5 * (9.4/3.)
      endcase
      widget_control, s.radius_btn_ann, set_value = a2str(s.radius_ann)
      update_plots, s
      plotimage, s
      gauss2dfit_image, s
    end
    'radius_right': begin
      s.radius += 1.
      widget_control, s.radius_btn, set_value = a2str(s.radius)
      case s.radius_unit of
        'pixel': s.radius_ann = s.radius + 1.5
        'arcsec': s.radius_ann = s.radius + 1.5 * (9.4/3.)
      endcase
      widget_control, s.radius_btn_ann, set_value = a2str(s.radius_ann)
      update_plots, s
      plotimage, s
      gauss2dfit_image, s
    end
    ;'no_annulus': begin
    ;  s.annulus_use = 0
    ;  widget_control, s.radius_panel_ann, sensitive = 0
    ;  ;print, 'Annulus subtraction is off'
    ;  update_plots, s
    ;  plotimage, s
    ;  gauss2dfit_image, s
    ;end
    ;'yes_annulus': begin
    ;  s.annulus_use = 1
    ;  widget_control, s.radius_panel_ann, sensitive = 1
    ;  ;print, 'Annulus subtraction is on'
    ;  update_plots, s
    ;  plotimage, s
    ;  gauss2dfit_image, s
    ;end
    'radius_left_ann': begin
      s.radius_ann -= 1.
      s.radius_ann = max([s.radius_ann, s.radius+1.])
      widget_control, s.radius_btn_ann, set_value = a2str(s.radius_ann)
      update_plots, s
      plotimage, s
      gauss2dfit_image, s
    end
    'radius_right_ann': begin
      s.radius_ann += 1.
      widget_control, s.radius_btn_ann, set_value = a2str(s.radius_ann)
      update_plots, s
      plotimage, s
      gauss2dfit_image, s
    end
    'radius_change': begin
      widget_control, s.radius_btn_ann, get_value = tmp
      s.radius_ann = max([tmp, s.radius+1.])
      update_plots, s
      plotimage, s
      gauss2dfit_image, s
    end
  endcase
  widget_control, ev.top, set_uvalue = s
end

;##############################################################################
;##############################################################################
;##############################################################################
;##############################################################################
;main routine
pro PACSman_analysis, object = object, save_plots = save_plots

  
  cleanplot, /sil
  set_plot, 'X'
  
  
  device,true_color=24
  device,decomposed=0
  device,retain=2
  
  if not keyword_set(object) then object = ''
  
  file_mkdir, 'ANALYSIS'
  file_mkdir, 'ANALYSIS/Cut'
  file_mkdir, 'ANALYSIS/Aperture'
  spawn, '\rm -f ANALYSIS/Aperture/*.{dat,png}'
  spawn, '\rm -f ANALYSIS/Cut/*.{dat,png}'
  
  loadct, 39, /silent
  
  spawn, '\find . -name "' + object + '*_Flux.fits"', files
  ind = where(strpos(files, 'iter4') gt -1, complement = noind)
  files = files(noind)
  if files(0) eq '' then begin
    print, 'No results!'
    print, 'the Flux map files are required for the analysis procedure'
    retall
  endif
  nlines = n_elements(files)
  
  lines = strarr(nlines)
  ;for i = 0, nlines - 1 do lines(i) = (strsplit(files(i), '/', /extract))(1)
  for i = 0, nlines - 1 do begin
    tmp = ['CII157', 'NII122', 'NII205', 'NIII57', 'OI145', 'OI63', 'OIII88']
    for j = 0, n_elements(tmp) - 1 do if strpos(files(i), tmp(j)) gt -1 then lines(i) = tmp(j)
    if lines(i) eq '' then print, 'LINE NOT FOUND!!!'
  endfor
  
  ind = where( strpos(files, 'UniqueResolution') eq -1 )
  files = files(ind)
  lines = lines(ind)
  nlines = n_elements(files)
  
  print, 'Lines observed:'
  print, lines
  
  ;------------------------------
  
  
  ;z = gauss2d(18, 20, [ 5., 10., nx/10.,  ny/6., nx/2., .6*ny])
  ;z = z + 0.1*randomn(seed, 18, 20)   ;Add random noise, SD = 1
  ;z = rot(z, 40)
  
  tmp = READFITS(files(0), h0,/silent)
  sxdelpar, h0, 'NAXIS3'
  sxaddpar, h0, 'NAXIS', 2
  im0 = reform(tmp(*, *, 0))
  err0 = reform(tmp(*, *, 1))
  ;err0 = READFITS(strtrans(files[0], '_FluxMap', '_NoiseMap'),/silent)
  tmp = reverse(im0, 1)
  S = size(tmp)
  xx = S[1]
  yy = S[2]
  xc = (xx-1)/2.
  yc = (yy-1)/2.
  result = FLTARR(xx, yy, nlines)
  noise = FLTARR(xx, yy, nlines)
  result(*, *, 0) = reverse(im0, 1)
  noise(*, *, 0) = reverse(err0, 1)
  if nlines gt 1 then begin
    for i = 1, nlines - 1 do begin
      tmp = READFITS(files(i), h1,/silent)
      im1 = reform(tmp(*, *, 0))
      ;print, sxpar(h1, 'NAXIS1')
      sxdelpar, h1, 'NAXIS3'
      sxaddpar, h1, 'NAXIS', 2
      hastrom, im1, h1, im1, htmp, h0, missing = !values.f_nan
      err1 = reform(tmp(*, *, 1))
      hastrom, err1, h1, err1, htmp, h0, missing = !values.f_nan
      result(*,*,i) = reverse(im1, 1)
      ;if i eq 0 then result(*, *, i) = z
      noise(*,*,i) = reverse(err1, 1)
    endfor
  endif
  
  extast, h0, astr
  xy2ad, 0, 0, astr, ra0, dec0
  xy2ad, xx, yy, astr, ra1, dec1
  ra = [ra1, ra0]
  dec = [dec0, dec1]
  
  x0 = 0
  y0 = (yy-1)/2.
  x1 = xx-1
  y1 = (yy-1)/2.
  
  sc_size = get_screen_size(resolution=sc_res)
  x_sc = min([0.95*sc_size(0), 1024])
  x_res = sc_res(0)
  y_sc = min([0.95*sc_size(1), 768])
  y_res = sc_res(1)
  
  main = widget_base (title='PACSman_analysis', /column, xsize=x_sc, ysize=y_sc, /base_align_top) ;main window
  wTab = widget_tab(main, uvalue='tab_event') ;tabs
  
  ;------------------------------
  
  wT1 = WIDGET_BASE(wTab, TITLE=' Map analysis ', /COLUMN)
  
  main_panel = widget_base (wT1, /row, xsize=0.99*x_sc, ysize=0.91*y_sc, xpad=0, ypad=0)
  
  left_panel = widget_base (main_panel, /column, xpad=10, ypad=10)
  top_panel = widget_base (left_panel, /row, xsize = 100)
  lines_map = widget_combobox (top_panel, uvalue='line_map', value=files, xsize=100, scr_xsize=430)
  tmp = widget_button(top_panel, value='Explode', uvalue='explode', xsize = 60)
  
  pos_panel = widget_base(left_panel, /row)
  pos_xy_text = widget_label(pos_panel, value='(+......, +......)')
  pos_radec1_text = widget_label(pos_panel, value=' (+......, +......)')
  pos_radec2_text = widget_label(pos_panel, value=' (..........................)')
  flux_text_otf = widget_label(pos_panel, value='-....e-..  ')
  image = widget_draw (left_panel, uvalue='image', xsize=0.5*x_sc, ysize=0.5*x_sc, /button, /motion_events)
  lbl = widget_label(left_panel, value='(plus=aperture center, cross=maximum within aperture, diamond=center of 2D gaussian)')
  lbl = widget_label(left_panel, value='(blue circle=PSF FWHM, black circle=FWHM of 2D gaussian)')
  ;lbl = widget_label(left_panel, value='(left click=change center, right click & drag=change profile cut)')
  radius_panel = widget_base(left_panel, /row)
  lbl = widget_label(radius_panel, value='Radius = ')
  radius_btn = widget_text(radius_panel, value='4.', xsize=6, /editable, uvalue='radius_change')
  radius_left = widget_button(radius_panel, value='<', uvalue='radius_left')
  radius_right = widget_button(radius_panel, value='>', uvalue='radius_right')
  radius_unit_panel = widget_base (radius_panel, column=2, /exclusive)
  radius_unit_pixel_btn = widget_button (radius_unit_panel, value=' (pixel)', uvalue='radius_unit_pixel')
  radius_unit_arcsec_btn = widget_button (radius_unit_panel, value=' (")', uvalue='radius_unit_arcsec')
  lbl = widget_label(radius_panel, value=' ')
  radius_panel_ann = widget_base(radius_panel, /row)
  lbl = widget_label(radius_panel, value=' Annulus = ')
  radius_btn_ann = widget_text(radius_panel, value='5.', xsize=6, /editable, uvalue='radius_ann_change')
  radius_left_ann = widget_button(radius_panel, value='<', uvalue='radius_left_ann')
  radius_right_ann = widget_button(radius_panel, value='>', uvalue='radius_right_ann')
  
  options_panel = widget_base(left_panel, /row)
  ;lbl = widget_label(options_panel, value='Sky subtraction: ')
  ;options_choices = widget_base(options_panel, column=2, /exclusive)
  ;yes_ann_btn = widget_button (options_choices, value='yes', uvalue='yes_annulus')
  ;no_ann_btn = widget_button (options_choices, value='no', uvalue='no_annulus')
  get_objects_btn = widget_button (options_panel, value='Get objects in aperture', uvalue='get_objects')
  
  
  ;gauss2dfit_btn = widget_button(left_panel, value="Gauss2Dfit within aperture", uvalue = "gauss2dfit")
  
  right_panel = widget_base (main_panel, /column, xpad=0, ypad=0)
  right_top_panel = widget_base (right_panel, /row, xpad=0, ypad=0)
  flux_growth = widget_draw (right_top_panel, uvalue='flux_growth', xsize=0.45*x_sc, ysize=0.45*y_sc)
  ;al_legend = widget_draw (right_top_panel, uvalue='none', xsize=0.15*x_sc, ysize=0.15*y_sc)
  sn_growth = widget_draw (right_panel, uvalue='sn_growth', xsize=0.45*x_sc, ysize=0.45*y_sc)
  
  
  ;------------------------------
  
  wT2 = WIDGET_BASE(wTab, TITLE=' Fluxes & surface fits ', /COLUMN)
  
  main_panel2 = widget_base (wT2, /column, xsize=0.99*x_sc, ysize=0.91*y_sc, xpad=0, ypad=0)
  
  image_panel = widget_base(main_panel2, /row)
  image_panel1 = widget_base(image_panel, /column)
  lbl = widget_label(image_panel1, value = 'Data')
  show3_data = widget_draw (image_panel1, uvalue='isurf', xsize=0.49*x_sc, ysize=0.4*x_sc, /button)
  image_panel2 = widget_base(image_panel, /column)
  lbl = widget_label(image_panel2, value = 'Fit')
  show3_fit = widget_draw (image_panel2, uvalue='isurf', xsize=0.49*x_sc, ysize=0.4*x_sc, /button)
  
  result_panel = widget_base(main_panel2, /row, /frame)
  flux_panel = widget_base(result_panel, /column, xsize=0.49*x_sc)
  flux_text = widget_text(flux_panel, value=' ', /scroll, ysize=15, /wrap)
  shape_panel = widget_base(result_panel, /column, xsize=0.49*x_sc)
  shape_text = widget_text(shape_panel, value=' ', /scroll, ysize=15, /wrap)
  
  
  ;------------------------------
  
  wT21 = WIDGET_BASE(wTab, TITLE=' Line spatial cut ', /COLUMN)
  
  main_panel21 = widget_base (wT21, /column, xsize=0.99*x_sc, ysize=0.91*y_sc, xpad=0, ypad=0)
  
  lbl = widget_label(main_panel21, value = 'Only points with flux>error are plotted')
  line_cut = widget_draw (main_panel21, uvalue='line_cut', xsize=0.8*x_sc, ysize=0.8*y_sc)
  
  
  ;------------------------------
  
  
  
  wT3 = WIDGET_BASE(wTab, TITLE=' About ', /COLUMN)
  
  main_panel3 = widget_base (wT3, /column, xsize=0.99*x_sc, ysize=0.91*y_sc, xpad=0, ypad=0)
  
  
  file = filepath(root_dir = ProgramRootDir(), subdir = 'misc/', 'pacsman_logo.bmp')
  if file_test(file) and not keyword_set(makeall) then begin
    test=QUERY_BMP(file, info)
    logo_draw = widget_draw(main_panel3, xsize = info.dimensions(0), ysize = info.dimensions(1))
  endif
    
  ;path = find_with_def('PACSman_analysis.pro', !path)
  ;tmp = strsplit(path, '/', /extract)
  ;path = '/'
  ;for i = 0, n_elements(tmp) - 2 do path = path + tmp[i] + '/'
  ;file = path+'misc/pacsman_logo.bmp'
  ;if file_test(file) then begin
  ;  test=QUERY_BMP(file, info)
  ;  logo_draw = widget_draw(main_panel3, xsize = info.dimensions(0), ysize = info.dimensions(1))
  ;endif
  
  versionfile = filepath(root_dir = ProgramRootDir(), 'pacsman_version.txt')
  readcol, versionfile, version, /silent
  
  ;ind = where(strpos((routine_info(/source))(*).path, 'pacsman_analysis') gt -1)
  ;pathto = (routine_info(/source))(ind(0)).path
  ;readcol, strtrans(pathto, 'pacsman_analysis.pro', 'pacsman_version.txt'), version, /silent
  
  str = ['PACSman version '+a2str(trim(version(0))), $
    '', $
    "IDL suite written by Vianney Lebouteiller & Diane Cormier, Service d'Astrophysique, CEA, Saclay, France", $
    ' ', $
    '- PACSman_fit: line fitting procedure ', $
    '- PACSman_linemap: overplots lines on an image ', $
    '- PACSman_map: projects maps on sky ', $
    '- PACSman_analysis: analyzes maps, flux extraction, ... ', $
    ' ', $
    'URL: http://www.myravian.fr/Homepage/Softwares.html', $
    ' ' $
    ]
  lbl = widget_text(main_panel3, value = str, /wrap, /scroll, xsize = 90, ysize = 20)
  
  ;------------------------------
  ;------------------------------
  
  widget_control, main, /realize
  
  
  widget_control, image, get_value=image_index ;for wset
  widget_control, flux_growth, get_value=flux_growth_index ;for wset
  ;widget_control, legend, get_value=legend_index ;for wset
  widget_control, line_cut, get_value=line_cut_index ;for wset
  widget_control, sn_growth, get_value=sn_growth_index ;for wset
  widget_control, show3_data, get_value=show3_data_index ;for wset
  widget_control, show3_fit, get_value=show3_fit_index ;for wset
  widget_control, logo_draw, get_value=logo_index ;for wset
  
  ;widget_control, no_ann_btn, set_button = 1
  
  if file_test(file) then begin
    device, get_decomposed=current_decomposed
    device, decomposed=1
    wset, logo_index
    ima = read_bmp(file)
    ima = reverse(temporary(ima), 1)
    tv, ima, true=1
    device, decomposed=current_decomposed
  endif
  
  loadct, 39, /silent
  
  tmp = widget_info(image, /geometry)
  drawsize = [tmp.xsize, tmp.ysize]
  
  widget_control, radius_unit_pixel_btn, set_button = 1
  
  widget_control, main, set_uvalue = { $
    bcd_arr: 0, $
    wset_image: image_index, $
    wset_flux_growth: flux_growth_index, $
    ;wset_legend: legend_index, $
    wset_line_cut: line_cut_index, $
    wset_sn_growth: sn_growth_index, $
    wset_show3_data: show3_data_index, $
    wset_show3_fit: show3_fit_index, $
    line_map: lines(0), $
    lines: lines, $
    files: files, $
    flux_text_otf: flux_text_otf, $
    result: result, $
    noise: noise, $
    ;annulus_use: 0, $
    exploded_index: 0, $
    exploded_view: 0, $
    exploded_draw: 0, $
    exploded_text: 0, $
    radius: 4., $
    radius_ann: 5., $
    ra: ra, $
    dec: dec, $
    drawsize: drawsize, $
    image_size: [xx, yy], $
    aper_center: [xc, yc], $
    header: h0, $
    image_index: 0, $
    nx: 0, $
    shape_text: shape_text, $
    ;radius_panel_ann: radius_panel_ann, $
    flux_text: flux_text, $
    radius_btn: radius_btn, $
    radius_btn_ann: radius_btn_ann, $
    radius_unit: 'pixel', $
    pos_xy_text: pos_xy_text, $
    pos_radec1_text: pos_radec1_text, $
    pos_radec2_text: pos_radec2_text, $
    line_cut_start: [x0, y0], $
    wtab: wtab, $
    line_cut_end: [x1, y1] $
    ;contrib: contrib $
    }
    
  widget_control, main, get_uvalue = s
  
  update_plots, s
  
  
  
  ;------------------------------
  
  
  
  ;------------------------------
  plotimage, s
  
  gauss2dfit_image, s
  
  
  
  
  xmanager, 'PACSman_analysis', main, /no_block ; wait for events
  
  
  
end

