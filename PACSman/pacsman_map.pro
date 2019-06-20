;+
; NAME:
;   PACSman_map
;
; AUTHOR:
;   Vianney Lebouteiller, Diane Cormier, CEA/SAp, Saclay, France
;   vianney.lebouteiller@cea.fr
;
;   UPDATED VERSIONs can be found at:
;      http://www.myravian.fr/Homepage/Softwares.html
;
; PURPOSE:
;   Makes projected maps from cube of line fluxes for Herschel/PACS data
;
; MAJOR TOPICS:
;   Map making
;
; CALLING SEQUENCE:
;   PACSman_map
;
; DESCRIPTION:
;
; INPUTS:
;
; RETURNS:
;
;   Nothing. Writes the map in dir.
;
;
; OPTIONAL PARAMETERS:
;
; REFERENCES:
;
; DEPENDENCIES:
;
; IDL ASTROLIB (readfits, writefits, getrot, extast, ad2xy, queryDSS) :
; http://idlastro.gsfc.nasa.gov/contents.html
;
; query2MASS (included) : http://irsa.ipac.caltech.edu
; queryIRAC (included)
; queryMIPS (included)
;
; MODIFICATION HISTORY:
;   07-28-10: calculations are now made in the raster reference frame, not
;   by rotating to RA/DEC axes. Added the calculation of the fraction
;   of subpixels that contribute to the projected grid
;   07-23-10: added error plane for PACSman_fill
;   Written, Jun 2010, VL
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
pro newct, ct, ownct = ownct, ncolors = ncolors;, noct = noct
  ;  help, ownct
  if keyword_set(ownct) then begin
    file_copy, !dir+'/resource/colors/colors1.tbl', '/tmp/pm_colors1.tbl', /allow, /force, /over
    file = '/tmp/pm_colors1.tbl'
    spawn, '\chmod 777 ' + file
  ;  help, file
  endif
  case ct of
    41: begin
      loadct, 33, /sil
      TVLCT, red, green, blue, /GET
      blue(0) = 0
      ;if not keyword_set(noct) then
      modifyct, 41, 'pacsman_map_41', red, green, blue, file = file
      loadct, 41, file = file, /sil, ncolors = ncolors
    end
    46: begin
      loadct, 6, /sil
      TVLCT, red, green, blue, /GET
      mid = 128
      ;red = findgen(256)
      ;blue = reverse(red)
      ;green = [255-2.*findgen(128)]
      ;green = [reverse(green), green]
      
      red = [fltarr(128)+255, 2*reverse(findgen(128))]
      blue = reverse(red)
      green = [255-2.*findgen(128)]
      green = [reverse(green), green]
      
      blue(0) = 0
      red(0) = 0
      green(0) = 0
      
      ;if not keyword_set(noct) then
      modifyct, 41, 'pacsman_map_46', red, green, blue, file = file
      loadct, 41, file = file, /sil, ncolors = ncolors
    end
    21: begin
      loadct, 22, /sil
      TVLCT, red, green, blue, /GET
      blue(0) = 0
      red(0) = 0
      green(0) = 0
      modifyct, 41, 'pacsman_map_47', red, green, blue, file = file
      loadct, 41, file = file, /sil, ncolors = ncolors
    end
  endcase
end

;##############################################################################
;--------------------------------------------------------------------------
function PACSman_simbad_getname, ra, dec

  sss = 'http://simbad.harvard.edu/simbad/sim-script?submit=submit+script&script=output+script=off+console=off%0D%0A'
  ;sss = 'http://simbad.u-strasbg.fr/simbad/sim-script?submit=submit+script&script=output+script=off+console=off%0D%0A'
  sss = sss + 'set+epoch+2000%0D%0A'
  sss = sss + 'format+object+%22SIMBADOBJECT%25-35IDLIST%281%29%22%0D%0A'
  sss = sss + 'query+coo+' + a2str(ra) + '+' + a2str(dec) + '+radius=0.09%0D'
  query = webget(sss)
  if strpos(query.text[0], 'error') gt -1 then begin
    print, 'No object found'
    return, ''
  endif
  object = a2str(strtrans(query.text[0], 'SIMBADOBJECT', ''))
  return, object
  
end
;--------------------------------------------------------------------------
;fill zeros within maps
function PACSman_fill, imag, h

  nx = n_elements(imag[*, 0])
  ny = n_elements(imag[0, *])
  ;print, nx, ny

  ind_zeros = where( finite(imag, /nan), count )
  
  if count eq 0. then return, imag
  
  for boxsize = 35, 7, -2 do begin
  ;boxsize = 25                    ;we check n x n patches only, n needs to be odd
  midboxsize = (boxsize-1) / 2
  for x = midboxsize, nx - midboxsize - 1 do for y = midboxsize, ny - midboxsize - 1 do begin
    box = imag[x-midboxsize:x+midboxsize, y-midboxsize:y+midboxsize]
    ind = where( finite(box, /nan), c )
    if c eq 0 or c eq boxsize^2 then continue ;nothing to fill or all nans
    ;print, x, y, c
    out = [ reform(box[0, 0:boxsize-1]), reform(box[boxsize-1, 0:boxsize-1]), reform(box[0:boxsize-1, 0]), reform(box[0:boxsize-1, boxsize-1])]
    ind = where( finite(out, /nan), c )
    if c gt 0 then continue ;at least one nan in the outer edge of the box
    print, boxsize, x, y
    newbox1 = box
    for x2 = 0, n_elements(box[*, 0]) - 1 do begin
       ind = where( finite(box[x2, *], /nan), c, complement = noind )
       if c gt 0 and n_elements(noind) gt 3 then newbox1[x2, ind] = interpol(box[x2, noind], noind, ind)
    endfor
    newbox2 = box 
    for y2 = 0, n_elements(box[0, *]) - 1 do begin
       ind = where( finite(box[*, y2], /nan), c, complement = noind )
       if c gt 0 and n_elements(noind) gt 3 then newbox2[ind, y2] = interpol(box[noind, y2], noind, ind)
    endfor
    newbox = 0.5 * ( newbox1 + newbox2 )
    imag[x-midboxsize:x+midboxsize, y-midboxsize:y+midboxsize] = newbox
    ;horizontal = imag(*, y)
    ;ind = where( finite(horizontal) )
    ;val1 = interpol(horizontal(ind), ind, x)
    ;vertical = imag(x, *)
    ;ind = where( finite(vertical) )
    ;val2 = interpol(vertical(ind), ind, y)
    ;imag(x, y) = avg([val1, val2])
  endfor
  endfor
  
  
  return, imag
  
end

;______________________________________________________________________________
;fringes + baseline + gaussian model
;function SpectrumModel, x, p, poly_degree = poly_degree
;  result = 0.
;  result += p[6] * sin( p[7]*x + p[8] )
;  for i = 0, poly_degree do result +=  p[i+9]*x^float(i)
;  return, p[0] * exp(-(x-p[1])^2./(2.*p[2]^2.)) + p[3] * exp(-(x-p[4])^2./(2.*p[5]^2.)) + result
;end
;______________________________________________________________________________
;fringes + baseline + gaussian model
function LineModel, x, p
  return, p[0] * exp(-(x-p[1])^2./(2.*p[2]^2.))
end
;-----------------------------------___________________________________________
;polynomial model
;function ContinuumModel, x, p, poly_degree = poly_degree
;  result = 0.
;  result += p[0] * sin( p[1]*x + p[2] )
;  for i = 0, poly_degree do result +=  p[i+3]*x^float(i)
;  return, result
;end
;______________________________________________________________________________
pro PACSman_acknowledge, text, parent_group = parent_group
  if not keyword_set(parent_group) then parent_group = ''
  wait_widget_ws = widget_base(title='PACSman, acknowledge', /column, xsize=400, ysize=50, xoffset=450, yoffset=300, group_leader = parent_group)
  label = widget_label(wait_widget_ws, value=text, /align_left)
  if parent_group ne '' then WIDGET_CONTROL, wait_widget_ws, /MODAL
  widget_control, wait_widget_ws, /realize
end
;##############################################################################
function PACSman_image, s, ra, dec, header = header, status = status

  status = 0
  
  if strpos(s.bgimage, '.fits') gt -1 then begin
  
    image = readfits(s.bgimage, header, /silent)
    extast, header, astr
    if size(astr, /tname) eq 'UNDEFINED' then image = readfits(s.bgimage, header, /silent, exten_no=1)
    
  endif else begin
  
    image = s.bgimage
    
    band_im = image
    tmp = strsplit(band_im, '/', /extract)
    if n_elements(tmp) ne 2 then begin
      print, 'Band could not be found, check help for more details'
      print, 'Input: ' + band_im
      retall
    endif
    tmp = tmp(1)
    
    if not file_test('retrieved_image_'+tmp+'.fits') then begin
    
      print, 'Retrieving image...'
      test = sm_webget('http://www.google.fr/index.html', timeout=1)
      if size(test, /tname) ne 'STRUCT' then if test eq '' then begin
        print, 'Internet connection seems to be down, aborting...'
        return, ''
      endif
      
      ra_center = ra
      dec_center = dec
      object = PACSman_simbad_getname(ra_center, dec_center)
      
      tmp = strsplit(band_im, '/', /extract)
      band_im = tmp(1)
      if band_im eq 'J' or band_im eq 'H' or band_im eq 'K' then query_2mass, [ra_center, dec_center], image, header, imsize=3.*(max(dec)-min(dec))*60., band = band_im else $
        if band_im eq 'ch1' or band_im eq 'ch2' or band_im eq 'ch3' or band_im eq 'ch4' then query_irac, object, image, header, band = band_im else $
        if band_im eq '24um' or band_im eq '70um' or band_im eq '160um' then query_mips, object, image, header, band = band_im else $
        queryDSS, [ra_center, dec_center], image, header, imsize=5.*(max(dec)-min(dec))*60., survey = band_im
        
      if size(image, /n_dim) le 1 then begin
        s.bgimage = ''
      endif else begin
        writefits, 'retrieved_image_'+band_im+'.fits', image, header
      endelse
      
    endif else begin
      image = readfits('retrieved_image_'+tmp+'.fits', header, /silent)
      object = sxpar(header, 'OBJECT')
    endelse
    
  endelse
  
  status = 1
  return, image
  
end
;______________________________________________________________________________
pro plotpng, map, header, xx, yy, filename = filename, s = s, velo = velo, fwhm = fwhm, unit_main = unit_main, unit_px = unit_px, fact = fact

  if keyword_set(velo) then ct = 46 else if keyword_set(fwhm) then ct = 21 else ct = 41
  
  newct, ct, ownct = s.ownct;, ncolors = 200
  extast, header, astr
  xy2ad, 0, 0, astr, ra0, dec0
  xy2ad, xx, yy, astr, ra1, dec1
  pacsman_disp, map, xrange = [ra1, ra0], xtitle = 'RA J2000 (deg)', $
    yrange = [dec0, dec1], position = [0.15, 0.1, 0.95, 0.9], ytitle = 'DEC J2000 (deg)', xsty = 1, ysty = 1, reserve = 1, charsize = 1.5, /squ
  loadct, 0, /sil
  
  pacsman_disp, map, xrange = [ra1, ra0], xtitle = 'RA J2000 (deg)', $
    yrange = [dec0, dec1], position = [0.15, 0.1, 0.95, 0.9], ytitle = 'DEC J2000 (deg)', xsty = 1, ysty = 1, /noplot, reserve = 1, charsize = 1.5, /squ, finalpos = finalpos
  s.finalpos = finalpos
  plots, [s.ra_center, s.ra_center], !y.crange, color = 100
  plots, !x.crange, [s.dec_center, s.dec_center], color = 100
  
  s.ra = [ra0, ra1]
  s.dec = [dec0, dec1]
  
  ncolors = 256
  loc = [0.15, 0.95, 0.95, 0.96]
  bar = BINDGEN(256) # REPLICATE(1B, 10)
  xsize = (loc(2) - loc(0)) * !D.X_VSIZE
  ysize = (loc(3) - loc(1)) * !D.Y_VSIZE
  xstart = loc(0) * !D.X_VSIZE
  ystart = loc(1) * !D.Y_VSIZE
  bar = BYTSCL(bar, TOP=ncolors-1)
  newct, ct, ncolors = ncolors, ownct = s.ownct
  IF !D.NAME EQ 'PS' THEN TV, bar, xstart, ystart, XSIZE=xsize, YSIZE=ysize ELSE TV, CONGRID(bar, xsize, ysize), xstart, ystart
  loadct, 0, /sil
  PLOTS, [loc(0), loc(0), loc(2), loc(2), loc(0)], [loc(1), loc(3), loc(3), loc(1), loc(1)], /NORMA
  plot, /noerase, position = loc, indgen(10), /nodata, xrange = minmax(map, /nan), ysty = 4, xsty = 1, charsize = 1.5
  
  xyouts, 0.01, 0.915, unit_main, /data, /normal, charsize = 1.5
  
  if keyword_set(unit_px) then begin
    axis, /xaxis, xrange = minmax(map*fact, /nan), ysty = 4, xsty = 1, charsize = 1.5
    xyouts, 0.01, 0.97, unit_px, /data, /normal, charsize = 1.5
  endif
  
  if keyword_set(filename) then write_png, filename, tvrd(/true)
end
;______________________________________________________________________________
pro getinp_event, ev ; DBL 8/1/2006
  widget_control,ev.id,get_uvalue=inpid
  widget_control,ev.top,get_uvalue=val_ptr
  if inpid ne 0L then begin
    widget_control, inpid, get_value=val
    *val_ptr=val
  endif else ptr_free,val_ptr
  widget_control, ev.top,/DESTROY
  return
end
function getinp, question, DEFAULT, TITLE=t,MODAL=md,PARENT_GROUP=group, $ ; DBL 8/1/2006
    _EXTRA=e
  if n_elements(t) eq 0 then t='Enter Value'
  if n_elements(md) eq 0 then md=0
  if n_elements(default) eq 0 then default=' '
  if n_elements(group) ne 0 then begin
    g=widget_info(group,/GEOMETRY)
    xoff=g.xoffset+g.scr_xsize/4
    yoff=g.yoffset+g.scr_ysize/4
  endif else begin
    device,get_screen_size=sc
    xoff=sc[0]/3 & yoff=sc[1]/3
  endelse
  val_ptr=ptr_new(DEFAULT)
  if n_elements(group) eq 0 then begin
    base=widget_base(TITLE=t,/COLUMN,xoffset=xoff,yoffset=yoff, uvalue=val_ptr)
  endif else begin
    base=widget_base(TITLE=t,/COLUMN,xoffset=xoff,yoffset=yoff, uvalue=val_ptr,MODAL=md, GROUP=group)
  endelse
  
  rowbase=widget_base(base,/row, /frame)
  input=cw_field(rowbase,value=DEFAULT,TITLE=question,/RETURN_EVENTS,_EXTRA=e)
  widget_control, input, set_uvalue=input ;put the id in
  rbase=widget_base(base,/row, /frame)
  ok=widget_button(rbase,value=' OK ',uvalue=input, xsize=100)
  cancel=widget_button(rbase,value='Cancel', uvalue=0L, xsize=100)
  widget_control, base,/realize,DEFAULT_BUTTON=ok
  widget_control, input,/INPUT_FOCUS
  XManager,'getinp',base ; DBL 8/1/2006
  if NOT ptr_valid(val_ptr) then return,''
  ret=*val_ptr
  ptr_free,val_ptr
  return, ret[0]
end
;______________________________________________________________________________
pro update_plots, s, force = force, all = all, ownct = ownct, ps = ps

  if s.active_file gt -1 or keyword_set(all) then begin
  
    if keyword_set(all) then indices = indgen(n_elements(s.files)) else begin
      widget_control, s.wtab, set_tab_current = 0
      indices = s.active_file
    endelse
    
    for index = 0, n_elements(indices) - 1 do begin
    
      s.active_file = indices(index)
      
      restore, s.files(s.active_file)
      
      ;test VIAN
      ;cube_lres.flux = cube_lres.flux * 0. + cube_lres.ra * 9. ;uniform distribution or ra/dec
      ;cube_lres.error = cube_lres.flux * 0. + cube_lres.ra * 9.
      ;cube_lres.flux = cube_lres.flux * 0. + 9. + (4.*randomu(seed, n_elements(cube_lres.flux))-2.) ;random distribution around 9
      
      ;already in MJy/sr
      ;cube_lres.continuum *= 1e-6
      ;cube_lres.continuum_error *= 1e-6
      
      ;convert in W m-2 sr-1 right away, simpler for later
      fact = 1. / (9.4/3600.*!pi/180.)^2.
      cube_lres.flux *= fact
      cube_lres.error *= fact
      if tag_exist(cube_params(0), 'pacsman_version') then if cube_params(0).pacsman_version ge 3.55 then begin
        print, 'PACSman_version >= 3.55'
        cube_lres.components_flux *= fact
        cube_lres.components_error *= fact
      endif
      cube_lres.integrated_flux *= fact
      cube_lres.integrated_noise *= fact
      
      ;widget_control, s.object_text, get_value = object
      object = s.object_name
      
      tmp = strsplit(s.files(s.active_file), '/', /extract)
      line = tmp(n_elements(tmp)-2)
      if not keyword_set(all) then widget_control, s.line_text, set_value = line
      dir = strtrans(s.files(s.active_file), tmp(n_elements(tmp)-1), '') + '/'
      
      if keyword_set(force) then begin
      
        band = s.band
        position_angle = s.posangle_value
        
      endif else begin ;read band and posangle
      
        band = cube_params(0).band
        ;help, band
        s.band = band
        if not keyword_set(all) then widget_control, s.band_text, set_value = band
        
      ;position_angle = 45. ;cube_params.position_angle
      ;help, position_angle
      ;if finite(position_angle) eq 0 then begin
      ;  if keyword_set(all) then position_angle = 118. else position_angle = float(getinp('Position angle not found, input value in degrees= ', s.posangle_value, parent_group = s.groupid))
      ;endif else s.posangle_value = position_angle
        
      endelse
      
      ;if not keyword_set(all) then widget_control, s.btn_update, sensitive = 1
      
      ;spaxel_size = transpose( [ [8.1, 8.8, 9.4, 8.7, 8.2], $
      ;  [9.4, 9.2, 9.1, 8.4, 7.6], $
      ;  [8.3, 8.4, 8.9, 9.1, 8.8], $
      ;  [8.1, 8.5, 8.9, 8.9, 8.9], $
      ;  [8.4, 8.4, 8.5, 8.9, 9.2] ] )
      ;spaxel_size = fltarr(5, 5) + 9.4
      
      pixel_dx = 9.4/3.             ;pixel size of the final map in arcsec
      ;resolution possible and
      ;because it's a third of a spaxel size
      pixel_dy = pixel_dx
      
      ;flat-field corrections
      file = filepath(root_dir = ProgramRootDir(), subdir = 'calib/', "spaxels_flatfield_"+band+".dat")
      print, 'Reading flat-field:', file
      readcol, file, spec_px, corr, /silent
      corr = transpose(reform(1./corr, 5, 5))
      
      ;PACSman_dir = find_with_def('pacsman_map.pro', !path)
      ;tmp = strsplit(PACSman_dir, '/', /extract)
      ;PACSman_dir = '/'
      ;for i = 0, n_elements(tmp) - 2 do PACSman_dir = PACSman_dir + tmp[i] + '/'
      ;if PACSman_dir eq '/' then PACSman_dir ='./'
      ;print, 'Reading flat-fields in:', PACSman_dir
      ;readcol, PACSman_dir+"/calib/spaxels_flatfield_"+band+".dat", spec_px, corr, /silent
      ;corr = transpose(reform(1./corr, 5, 5))
      
      
      ;if a2str(object(0)) eq 'LMC-N11B' then begin
      ;  help, object
      ;  print, "Manual correction of flat-field for SD data of LMC-N11B"
      ;  case band of
      ;    'red': begin
      ;      ;corr(1, 0) = 0.8
      ;      print, corr(0, 4), corr(1, 0), corr(2, 0), corr(1, 1), corr(2, 2)
      ;      ;corr(0, 4) = 1.
      ;      ;corr(1:2, 0) = 0.5 ;OI 145 and CII 157 in N11B
      ;      ;corr(1, 1) = 0.7 ;OI 145 and CII 157 in N11B
      ;      print, corr(0, 4), corr(1, 0), corr(2, 0), corr(1, 1), corr(2, 2)
      ;    ;corr(0, 3) = 0.7 ;from NII122 in N11B, doesn't agree with CII, due to fringes
      ;    end
      ;    'blue': begin
      ;    ;corr(4, 2) = 0.8 ;from NII 57 in N11B
      ;    ;corr(1, 1) = 0.8 ;from NII 57 in N11B
      ;    end
      ;  endcase
      ;;corr(2, 0) = 10.
      ;endif
      
      if s.niter gt 0 then print, 'number of iterations: ', s.niter
      
      ;spawn, '\ls ' + dir+s.object_name+'_*.{png,fits,sav}', files
      ;if files(0) ne '' then for f = 0, n_elements(files) - 1 do if strpos(files(f), 'Line_map') eq -1 then spawn, '\rm -rf ' + files(f)
      
      wave = cube_lres.wave
      
      ;HERE I CAN CHANGE WHETHER WE WANT THE NOMINAL COMPONENT OR THE 2ND OR THE 3RD, as an option
      help, cube_params(0), /st
      ;help, cube_params(0).pacsman_version
      if tag_exist(cube_params(0), 'pacsman_version') then begin
         if cube_params(0).pacsman_version ge 3.55 then begin
        print, 'PACSman_version >= 3.55'
        help, s.ncomp
        if s.ncomp eq 0 then begin
           if s.integrated eq 1 then begin
              flux = cube_lres.integrated_flux
              noise = cube_lres.integrated_noise
           endif else begin
              flux = cube_lres.flux
              noise = cube_lres.error
           endelse
        endif else begin
           if s.ncomp-1 gt n_elements(cube_lres(0).components_flux) - 1 then begin
              print, "Wrong component number"
              retall
           endif else begin
              print, 'Projecting component: ', s.ncomp
              flux = cube_lres.components_flux(s.ncomp-1)
              noise = cube_lres.components_error(s.ncomp-1)
           endelse
        endelse
     endif else begin 
        print, 'PACSman_version < 3.55'
        if s.integrated eq 1 then begin
           flux = cube_lres.integrated_flux
           noise = cube_lres.integrated_noise
        endif else begin
           flux = cube_lres.flux
           noise = cube_lres.error
        endelse
     endelse
  endif else begin
        print, 'no PACSman version tag'
        if s.integrated eq 1 then begin
           flux = cube_lres.integrated_flux
           noise = cube_lres.integrated_noise
        endif else begin
           flux = cube_lres.flux
           noise = cube_lres.error
        endelse
  endelse
      
      ;if keyword_set(use_calib_error) then begin
      ;case band of
      ;  'red': errcal = 0.2 * flux
      ;  'blue': errcal = 0.1 * flux
      ;endcase
      ;errcal += 0.1 * flux
      errcal = (0.01*s.use_calib_error) * flux
      noise += errcal
      help, 0.01*s.use_calib_error
      
      if s.perturb eq 1 then begin
        print, 'perturb!'
        flux += randomn(seed, n_elements(flux)) * noise
      ;print, flux
      ;stop
      endif
      
      velocity = cube_lres.velocity
      velocity_error = cube_lres.velocity_error
      linewidth = cube_lres.broadening
      linewidth_error = cube_lres.broadening_error
      continuum = cube_lres.continuum
      continuum_error = cube_lres.continuum_error
      
      nrasters = total(cube_params.n_rasters) ;number of rasters
      ;for n = 1, nrasters do for x = 0, 4 do for y = 0, 4 do begin     ;correct by spatial flat-field
      ;  ind = where( cube_lres.spaxel_X eq (x+1) and cube_lres.spaxel_Y eq (y+1) and cube_lres.raster eq n )
      ;flux(ind) *= corr(x, y)
      ;noise(ind) *= corr(x, y)
      ;print, 'nocorr'
      ;endfor
      
      ra = double(cube_lres.ra)
      dec = double(cube_lres.dec)
      s.ra_center = double(0.5*total(minmax(ra)))
      s.dec_center = double(0.5*total(minmax(dec)))
      
      ;##################################################################
      ;print, ra_center, dec_center
      ;help, '2'
      if keyword_set(all) then begin
        set_plot, 'Z'
        ERASE
        device, set_resolution=[800,600], set_pixel_depth = 24, decomposed = 0
      endif else wset, s.wset_image_footprints
      loadct, 0, /silent
      ;help, '3'
      
      edgey = 35./3600.
      edgex = 35./3600. / cos(s.dec_center*!pi/180.)
      
      if s.bgimage ne '' then begin
      
        image = PACSman_image(s, s.ra_center, s.dec_center, header = header, status = status)
        
        if status eq 1 then begin
        
          if get_equinox(header) ne 2000. and get_equinox(header) ne 0. then hprecess, header, 2000.
          getrot, header, rot_angle
          extast, header, astr
          ad2xy, s.ra_center, s.dec_center, astr, x, y
          hrot, image, header, new_image, new_header, rot_angle, x, y, 1
          
          extast, new_header, astr
          ad2xy, max(ra)+edgex, min(dec)-edgey, astr, x1, y1
          ad2xy, min(ra)-edgex, max(dec)+edgey, astr, x2, y2
          
          if x1 gt 0 and x2 lt n_elements(image(*, 0)) and y1 gt 0 and y2 lt n_elements(image(0, *)) then begin
          
            hextract, new_image, new_header, image, header, x1, x2, y1, y2
            
            pacsman_disp, image, xrange = [edgex+max(ra), min(ra)-edgex], xstyle=1, xtitle = 'RA J2000 (deg)', $
              yrange = [min(dec)-edgey, edgey+max(dec)], ystyle=1, ytitle = 'DEC J2000 (deg)', /squ, $
              position = [0.15, 0.1, 0.97, 0.97], /logscale, charsize = 1.5
              
          endif else begin
            print, 'Cut outside image!'
            plot, indgen(10), /nodata, xrange = [edgex+max(ra), min(ra)-edgex], yrange = [min(dec)-edgey, edgey+max(dec)], $
              xtitle = 'RA (J2000)', ytitle = 'DEC (J2000)', position = [0.15, 0.1, 0.97, 0.97]
          endelse
          
        endif
        
      endif else begin ;no image to plot
      
        plot, indgen(10), /nodata, xrange = [edgex+max(ra), min(ra)-edgex], yrange = [min(dec)-edgey, edgey+max(dec)], $
          xtitle = 'RA (J2000)', ytitle = 'DEC (J2000)', position = [0.15, 0.1, 0.97, 0.97]
          
      endelse
      
      plots, [s.ra_center, s.ra_center], !y.crange, color = 100
      plots, !x.crange, [s.dec_center, s.dec_center], color = 100
      
      loadct, 38, /silent
      text = strarr(nrasters)
      cols = intarr(nrasters)
      for n = 1, nrasters do begin ;overplot rasters
        cols(n-1) = (n * 25) mod 255
        ind = where( cube_lres.raster eq n )
        oplot, ra(ind), dec(ind), psym = sym(1), color = cols(n-1), symsize = 0.5
        ind = where( cube_lres.raster eq n and cube_lres.spaxel_X eq 1 and cube_lres.spaxel_Y eq 1 )
        plots, ra(ind), dec(ind), psym = 4, color = cols(n-1), symsize = 2.
        ind = where( cube_lres.raster eq n and cube_lres.spaxel_X eq 1 and cube_lres.spaxel_Y eq 2 )
        plots, ra(ind), dec(ind), psym = 4, color = cols(n-1), symsize = 1.
        text(n-1) = '#' + a2str(n)
        ind = where( cube_lres.raster eq n and cube_lres.spaxel_X eq 3 and cube_lres.spaxel_Y eq 3 )
        xyouts, ra(ind), dec(ind), a2str(n), color = cols(n-1), alignment = 0
      ;wait, 0.07
      endfor
      al_legend, text, textcol = cols, /top, /left, box = 0, charsize = 1
      al_legend, ['spaxel 1, 1', 'spaxel 1, 2'], psym = [4, 4], symsize = [2, 1], box = 0, /top, /right
      write_png, dir+s.object_name+'_LINE'+line+'_Footprints.png', tvrd(/true)
      
      ;##################################################################
      
      ;pa = position_angle * !pi / 180. ;position angle in radians
      
      ;ra and dec are from the spaxels, not the sub-pixels
      dx = double( (ra - s.ra_center) * 3600. * cos(dec*!pi/180.) )
      dy = double( (dec - s.dec_center) * 3600. )
      
      ;turn to orthogonal
      ;rot_matrix = [ [cos(pa), -sin(pa)], [sin(pa), cos(pa)] ]
      ;for i = 0, n_elements(dx) - 1 do begin
      ;  a = [dx(i), dy(i)] # rot_matrix
      ;  a = reform(a)
      ;  dx(i) = -a(0)
      ;  dy(i) = -a(1)
      ;endfor
      
      nx_right = ceil((max(dx, /nan)+sqrt(2.)*pixel_dx)/pixel_dx)
      nx_left = floor((min(dx, /nan)-sqrt(2.)*pixel_dx)/pixel_dx)
      
      ny_top = ceil((max(dy, /nan)+sqrt(2.)*pixel_dy)/pixel_dy)
      ny_bottom = floor((min(dy, /nan)-sqrt(2.)*pixel_dy)/pixel_dy)
      
      mm = long(nx_right - nx_left)
      nn = long(ny_top - ny_bottom)
      
      ;##################################################################
      
      npoints = n_elements(cube_spectra(0).flux)
      ;print, '(debug) array size = ', double(npoints*mm*nn) / 1.e8
      ;print, keyword_set(all)
      if double(npoints*mm*nn) gt 1.e8 or double(npoints*mm*nn) lt 0. then begin
        if not keyword_set(all) then pm_acknowledge, text='no spectral projection will be performed, too many points'
        npoints = 1
      endif
      
      for iter = 0, s.niter do begin
      
        print, 'Creating empty projected grid...'
        grid = replicate({mm: 0, nn: 0, x: 0., y: 0., flux: !values.f_nan, error: !values.f_nan, velocity: !values.f_nan, velocity_error: !values.f_nan, weight: 0., $
          weight_velocity: 0., spatial_resolution: !values.f_nan, ra: double(!values.f_nan), dec: double(!values.f_nan), fwhm: !values.f_nan, fwhm_error: !values.f_nan, $
          weight_fwhm: 0., continuum: !values.f_nan, continuum_error: !values.f_nan, weight_continuum: 0., spec_w: !values.f_nan*dblarr(npoints), spec_f: !values.f_nan*dblarr(npoints)}, long(nn*mm))
        k = long(0)
        for i = 0, mm - 1 do for j = 0, nn - 1 do begin
          grid(k).mm = i
          grid(k).nn = j
          grid(k).x = (i+nx_left+.5) * pixel_dx
          grid(k).y = (j+ny_bottom+.5) * pixel_dy
          k += 1
        endfor
        print, '...done'
        
        ;grid = fltarr(mm, nn, 17)*!values.f_nan      ;x, y, flux, error, velocity, velocity_error, area_ratio, area_ratio_velo, spatial resolution, ra, dec, ;
        ;fwhm, fwhm_error, area_ratio_fwhm, cont, cont_error, area_ratio_cont
        
        ;--------------------------------------------------------------------------
        if iter eq 0 and not keyword_set(all) then begin
          wait, 1.5
          widget_control, s.wtab, set_tab_current = 1
        endif
        if keyword_set(all) then begin
          set_plot, 'Z'
          ERASE
          device, set_resolution=[800,800], set_pixel_depth = 24, decomposed = 0
        endif else wset, s.wset_image_subpixels
        loadct, 0, /silent
        
        ;range = [min([minx-2.*pixel_dx, miny-2.*pixel_dy]), max([maxx+2.*pixel_dx, maxy+2.*pixel_dy])]
        ;needs to be a square
        range = [min([nx_left*pixel_dx, ny_bottom*pixel_dy]), max([nx_right*pixel_dx, ny_top*pixel_dy])]
        
        plot, indgen(10), /nodata, xran = reverse(range), yran = range, xtitle = 'dx (")', ytitle = 'dy (")', xsty = 1, ysty = 1, $
          position = [0.1, 0.1, .97, .97]
        ;pacsman_disp, fltarr(10,10), position = [0.1, 0.1, 1., 1.], xran = reverse(range), yran = range, xtitle = 'dx (")', $
        ;  ytitle = 'dy (")', xsty = 1, ysty = 1, min = 0, max = 10
          
          
        for i = 0, mm-1 do begin
          for j = 0, nn-1 do begin
            loadct, 0, /sil
            xl = (i + nx_left) * pixel_dx
            xr = xl + pixel_dx
            yd = (j + ny_bottom) * pixel_dy
            yu = yd + pixel_dy
            plots, [xl, xr], [yu, yu], linestyle = 0, color = 50
            plots, [xl, xr], [yd, yd], linestyle = 0, color = 50
            plots, [xr, xr], [yd, yu], linestyle = 0, color = 50
            plots, [xl, xl], [yd, yu], linestyle = 0, color = 50
          endfor
        endfor
        
        for i = 0, mm-1 do begin
        
          print, 'i = ' + a2str(i) + '/' + a2str(mm-1)
          
          for j = 0, nn-1 do begin
          
            ;print, 'i = ' + a2str(i) + '/' + a2str(mm-1);, ', j = ' + a2str(j) + '/' + a2str(nn-1)
          
            ind_grid = where( grid.mm eq i and grid.nn eq j )
            
            dxp = grid(ind_grid).x
            dyp = grid(ind_grid).y
            
            coords = [double(dxp/(3600.*cos((s.dec_center+dyp/3600.)*!pi/180.))), double(dyp/3600.)]
            coords += [double(s.ra_center), double(s.dec_center)]
            
            grid(ind_grid).ra = a2str(coords[0])
            grid(ind_grid).dec = a2str(coords[1])
            if i eq mm/2 and j eq nn/2 then begin
              coords_ref = coords
              xy_ref = [i+1, j+1] ;;WCS is 1 when IDL is 0
            endif
            
            for n = 1, nrasters do begin
            
              for x = 0, 4 do for y = 0, 4 do begin
              
                ;index of spaxel for raster n
                ind_spaxel = where( cube_lres.spaxel_X eq (x+1) and cube_lres.spaxel_Y eq (y+1) and cube_lres.raster eq n )
                
                ;gets dx and dy from the spaxel center
                xp = dx(ind_spaxel)
                yp = dy(ind_spaxel)
                
                ;index of cube, to get position angle
                ind_cube = 0
                while n gt total(cube_params(0:ind_cube).n_rasters) do ind_cube += 1
                
                if finite(cube_params(ind_cube).position_angle, /nan) then pa = s.posangle_value * !pi/180. else pa = cube_params(ind_cube).position_angle * !pi/180.
                
                rot_matrix = [ [cos(-pa), -sin(-pa)], [sin(-pa), cos(-pa)] ]
                
                spaxel_size = 3. ;in pixel_dx units
                bl = [-.5*spaxel_size*pixel_dx, -.5*spaxel_size*pixel_dy] # rot_matrix
                ul = [-.5*spaxel_size*pixel_dx, .5*spaxel_size*pixel_dy] # rot_matrix
                br = [.5*spaxel_size*pixel_dx, -.5*spaxel_size*pixel_dy] # rot_matrix
                ur = [.5*spaxel_size*pixel_dx, .5*spaxel_size*pixel_dy] # rot_matrix
                
                bl += [xp, yp]
                ul += [xp, yp]
                br += [xp, yp]
                ur += [xp, yp]
                
                ;loadct, 0, /sil
                ;polyfill, [bl(0), br(0), ur(0), ul(0), ul(0)], [bl(1), br(1), ur(1), ul(1), ul(1)], color = 50
                
                xx = [bl(0), br(0), ur(0), ul(0)]
                yy = [bl(1), br(1), ur(1), ul(1)]
                
                xx /= pixel_dx ;arcsec -> px
                yy /= pixel_dy ;arcsec -> px
                
                distance = sqrt( (i+nx_left+.5 - xp/pixel_dx)^2. + (j+ny_bottom+.5 - yp/pixel_dy)^2. )
                
                polyclip_pm, i+nx_left, j+ny_bottom, xx, yy, status = status
                
                if status eq 1 then begin
                
                  if iter gt 0 then sig = ((1./iter)*flux(ind_spaxel) + iter*grid_ref(ind_grid).flux)/(iter+1./iter) else sig = flux(ind_spaxel)
                  err = noise(ind_spaxel)
                  velo = velocity(ind_spaxel)
                  veloerr = velocity_error(ind_spaxel)
                  fwhm = linewidth(ind_spaxel)
                  fwhmerr = linewidth_error(ind_spaxel)
                  cont = continuum(ind_spaxel)
                  conterr = continuum_error(ind_spaxel)
                  
                  loadct, 38, /sil
                  
                  plots, [xx, xx(0)]*pixel_dx, [yy, yy(0)]*pixel_dy, thick = 2, color = (n * 25) mod 255
                  
                  
                  ;plots, xx*pixel_dx, yy*pixel_dy, ps = 4, color = 150
                  area_ratio = poly_area(xx, yy)
                  weight = area_ratio ;;/ (err / median(noise)) ;/ distance
                  weight = weight(0)
                  
                  if not finite(grid(ind_grid).continuum) then begin ;initialize from NaN to 0
                    ;print, '0! ', ind_grid, n, x, y
                    grid(ind_grid).spec_w = 0.
                    grid(ind_grid).spec_f = 0.
                    grid(ind_grid).flux = 0.
                    grid(ind_grid).error = 0.
                    grid(ind_grid).velocity = 0.
                    grid(ind_grid).velocity_error = 0.
                    grid(ind_grid).fwhm = 0.
                    grid(ind_grid).fwhm_error = 0.
                    grid(ind_grid).continuum = 0.
                    grid(ind_grid).continuum_error = 0.
                  endif
                  
                  indok = where( finite(cube_spectra(ind_spaxel).flux) and cube_spectra(ind_spaxel).flux ne 0. and cube_spectra(ind_spaxel).lambda_obs ne 0., c )
                  ;spectrum in MJy sr-1
                  ;we need to interpolate because of the automatic sampling
                  if c gt 0 then begin
                    grid(ind_grid).spec_f(indok) += weight * interpol(cube_spectra(ind_spaxel).flux(indok), cube_spectra(ind_spaxel).lambda_obs(indok), cube_spectra(0).lambda_obs(indok))
                    grid(ind_grid).spec_w(indok) = cube_spectra(0).lambda_obs(indok)
                  endif
                  
                  grid(ind_grid).flux += sig * weight ;flux, in W m-2 sr-1, so size of pixel doesn't matter
                  grid(ind_grid).error += err^2. * weight ;err * weight ;error
                  
                  ;do we take the quadratic sum or just the sum?
                  ;several spaxels are actually not measuring the same quantity for a given sub-pixel in the projected grid, the fluxes
                  ; are somewhat different because they correspond to a somewhat different coverage, the errors should probably just be the sum then...
                  grid(ind_grid).weight += weight                    ;for the weights normalization
                  
                  if finite(velo) then begin;and table(ind(i)).velocity_error gt 0. then begin
                    grid(ind_grid).velocity += velo * weight ;velocity
                    grid(ind_grid).velocity_error += veloerr^2. * weight ;veloerr * weight ;velocity error
                    grid(ind_grid).weight_velocity += weight ;velocity
                  ;print, grid(ind_grid).velocity, grid(ind_grid).area_ratio_velocity
                  endif
                  
                  if finite(fwhm) then begin; and table(ind(i)).fwhm_error gt 0. then begin
                    grid(ind_grid).fwhm += fwhm * weight ;fwhm
                    grid(ind_grid).fwhm_error += fwhmerr^2. * weight ;fwhmerr * weight ;fwhm error
                    grid(ind_grid).weight_fwhm += weight ;fwhm
                  endif
                  
                  if finite(cont) then begin; and table(ind(i)).continuum_error gt 0. then begin
                    grid(ind_grid).continuum += cont * weight ;cont
                    grid(ind_grid).continuum_error += conterr^2. * weight ;conterr * weight ;cont error
                    grid(ind_grid).weight_continuum += weight ;cont
                  endif
                  
                  loadct, 0, /sil
                  xl = (i + nx_left) * pixel_dx
                  xr = xl + pixel_dx
                  yd = (j + ny_bottom) * pixel_dy
                  yu = yd + pixel_dy
                  plots, [xl, xr], [yu, yu], linestyle = 0, color = 100
                  plots, [xl, xr], [yd, yd], linestyle = 0, color = 100
                  plots, [xr, xr], [yd, yu], linestyle = 0, color = 100
                  plots, [xl, xl], [yd, yu], linestyle = 0, color = 100
                  
                  loadct, 38, /silent
                  if x eq 0 and y eq 0 then plots, xp, yp, psym = 4, color = (n * 25) mod 255, symsize = 2. else $
                    if x eq 0 and y eq 1 then plots, xp, yp, psym = 4, color = (n * 25) mod 255, symsize = 1. else $
                    plots, xp, yp, psym = 6, color = (n * 25) mod 255, symsize = 0.5
                    
                endif  ;contribution
                
              endfor ;spaxel
              
            endfor ;raster
            
            grid(ind_grid).spec_f /= grid(ind_grid).weight
            
            grid(ind_grid).flux /= grid(ind_grid).weight                     ;normalization
            grid(ind_grid).error = sqrt( grid(ind_grid).error / grid(ind_grid).weight )  ;quadratic sum for errors
            ;grid(ind_grid).error /= grid(ind_grid).weight  ;non-quadratic sum for errors
            
            grid(ind_grid).velocity /= grid(ind_grid).weight_velocity                    ;normalization
            grid(ind_grid).velocity_error = sqrt( grid(ind_grid).velocity_error / grid(ind_grid).weight_velocity ) ;quadratic sum for errors
            ;grid(ind_grid).velocity_error /= grid(ind_grid).weight_velocity ;non-quadratic sum for errors
            
            grid(ind_grid).fwhm /= grid(ind_grid).weight_fwhm                     ;normalization
            grid(ind_grid).fwhm_error = sqrt( grid(ind_grid).fwhm_error / grid(ind_grid).weight_fwhm ) ;quadratic sum for errors
            ;grid(ind_grid).fwhm_error /= grid(ind_grid).weight_fwhm ;non-quadratic sum for errors
            
            grid(ind_grid).continuum /= grid(ind_grid).weight_continuum                     ;normalization
            grid(ind_grid).continuum_error = sqrt( grid(ind_grid).continuum_error / grid(ind_grid).weight_continuum ) ;quadratic sum for errors
            ;grid(ind_grid).continuum_error /= grid(ind_grid).weight_continuum ;non-quadratic sum for errors
            
          endfor
          
        endfor ;grid
        
        print, 'Making plots, saving files...'
        
        ;center position
        loadct, 0, /silent
        plots, !x.crange, [0., 0.], linestyle = 0, color = 100
        plots, [0., 0.], !y.crange, linestyle = 0, color = 100
        
        loadct, 0, /silent
        
        outputname = line
        if tag_exist(cube_params(0), 'pacsman_version') then if cube_params(0).pacsman_version ge 3.55 then if s.ncomp gt 0 then outputname += '_comp' + a2str(s.ncomp+1)
        
        ;plot results
        plot, indgen(10), /nodata, xran = reverse(range), yran = range, xtitle = 'dx (")', ytitle = 'dy (")', xsty = 1, $
          ysty = 1, /noerase, position = [0.1, 0.1, .97, .97]
        write_png, dir+s.object_name+'_LINE'+outputname+'_Projected_grid.png', tvrd(/true)
        ;widget_control, s.image_subpixels_text, set_value=dir+s.object_name+'_LINE'+line+'_subpixel_projected_map.png'
        
        
        if keyword_set(all) then begin
          set_plot, 'Z'
          ERASE
          device, set_resolution=[800,800], set_pixel_depth = 24, decomposed = 0
        endif else wset, s.wset_image_subpixels2
        s.dx = minmax(grid.x)
        s.dy = minmax(grid.y)
        imag = fltarr(mm, nn)
        xin = fltarr(mm, nn)
        yin = fltarr(mm, nn)
        for i = 0, mm - 1 do for j = 0, nn - 1 do begin
          ind_grid = where( grid.mm eq i and grid.nn eq j )
          imag(i, j) = grid(ind_grid).flux
          xin(i, j) = grid(ind_grid).x
          yin(i, j) = grid(ind_grid).y
        endfor
        imag = reverse(imag, 1)
        
        pacsman_disp, imag, xin, yin, xrange = reverse(range), yrange = range, xtitle = 'dx (")', ytitle = 'dy (")', $
          xsty = 1, ysty = 1, position = [0.1, 0.1, .97, .97], charsize = 1., finalpos = finalpos
        ;pacsman_disp, imag, 1, 1, xrange = reverse(minmax(grid.x)), yrange = minmax(grid.y), xtitle = 'dx (")', ytitle = 'dy (")', $
        ;  xsty = 1, ysty = 1, position = [0.15, 0.1, .95, .9], charsize = 1.5, /square, finalpos = finalpos
        s.finalpos_grid = finalpos
        
        k = 0
        loadct, 38, /silent
        for n = 1, nrasters do begin
          for x = 0, 4 do for y = 0, 4 do begin
            ind = where( cube_lres.spaxel_X eq (x+1) and cube_lres.spaxel_Y eq (y+1) and cube_lres.raster eq n )
            xp = dx(ind)
            yp = dy(ind)
            if x eq 0 and y eq 0 then plots, xp, yp, psym = 4, color = (n * 25) mod 255, symsize = 2. else $
              if x eq 0 and y eq 1 then plots, xp, yp, psym = 4, color = (n * 25) mod 255, symsize = 1. else $
              plots, xp, yp, psym = 6, color = (n * 25) mod 255, symsize = 0.5
          endfor
        endfor
        write_png, dir+s.object_name+'_LINE'+outputname+'_Projected_grid_flux.png', tvrd(/true)
        
        ;##################################################################
        
        result = grid
        
        nresult = fltarr(mm,nn)
        nnoise = nresult
        velomap = nresult
        velomap_error = nresult
        fwhmmap = nresult
        fwhmmap_error = nresult
        contmap = nresult
        contmap_error = nresult
        for i = 0, mm - 1 do for j = 0, nn - 1 do begin
          ind_grid = where( grid.mm eq i and grid.nn eq j )
          nresult(i, j) = grid(ind_grid).flux
          nnoise(i, j) = grid(ind_grid).error
          velomap(i, j) = grid(ind_grid).velocity
          velomap_error(i, j) = grid(ind_grid).velocity_error
          fwhmmap(i, j) = grid(ind_grid).fwhm
          fwhmmap_error(i, j) = grid(ind_grid).fwhm_error
          contmap(i, j) = grid(ind_grid).continuum
          contmap_error(i, j) = grid(ind_grid).continuum_error
        endfor
        
        ;in W m-2 sr-1
        pixel_sr = pixel_dx^2. * 2.3504e-11 ; 2.3504e-11 sr = 1 arcsec^2.
        
        ;nresult /= pixel_sr
        ;nnoise /= pixel_sr
        
        
        ;MAKE HEADER
        header = strarr(1)
        
        sxaddpar, header, 'SIMPLE', 2
        sxaddpar, header, 'NAXIS', 3
        sxaddpar, header, 'NAXIS1', mm
        sxaddpar, header, 'NAXIS2', nn
        sxaddpar, header, 'NAXIS3', 4
        
        sxaddpar, header, 'OBJECT', object(0)
        sxaddpar, header, 'RADESYS', 'ICRS    ', ' International Celestial Reference System'
        sxaddpar, header, 'EQUINOX', 2000.
        
        sxaddpar, header, 'CRVAL1', coords_ref(0)
        sxaddpar, header, 'CRPIX1', xy_ref(0)
        sxaddpar, header, 'CDELT1', pixel_dx/3600. ;*cos(coords_ref(0)*!pi/180.)
        sxaddpar, header, 'CTYPE1', 'RA---TAN'
        
        sxaddpar, header, 'CRVAL2', coords_ref(1)
        sxaddpar, header, 'CRPIX2', xy_ref(1)
        sxaddpar, header, 'CDELT2', pixel_dy/3600.
        sxaddpar, header, 'CTYPE2', 'DEC--TAN'
        
        sxaddpar, header, 'BUNIT', 'W/m2/sr'
        
        object = strtrans(object, ' ', '')
        sxaddpar, header, 'CROTA1', 0;position_angle
        sxaddpar, header, 'CROTA2', 0;-position_angle
        ;print, header
        
        ind = where(strpos((routine_info(/source))(*).path, 'pacsman_map') gt -1)
        pathto = (routine_info(/source))(ind(0)).path
        sxaddpar, header, 'PACSman_p', pathto
        readcol, strtrans(pathto, 'pacsman_map.pro', 'pacsman_version.txt'), version, /silent
        version = version(0)
        sxaddpar, header, 'PACSman_v', a2str(trim(version))
        sxaddpar, header, 'Date', systime()
        if tag_exist(cube_params(0), 'pacsman_version') then if cube_params(0).pacsman_version ge 3.51 then begin
          sxaddpar, header, 'HIPE', trim(cube_params(0).hipe_version)
          str = ''
          for si = 0, n_elements(cube_params) - 1 do str += trim(cube_params(si).obsid) + ' '
          sxaddpar, header, 'OBS_ID', trim(str)
        endif
        
        ;WRITE FITS FILES
        ;SPECTRA
        spec_w = fltarr(mm, nn, npoints)
        spec_f = spec_w
        np = 0
        for i = 0, mm - 1 do for j = 0, nn - 1 do begin
          ind_grid = where( grid.mm eq i and grid.nn eq j )
          ind_tmp = where( grid(ind_grid).spec_w gt 0. and finite(grid(ind_grid).spec_w), c )
          if c gt np then np = c
          spec_w(i, j, *) = grid(ind_grid).spec_w
          spec_f(i, j, *) = grid(ind_grid).spec_f
        endfor
        if np eq 0 then begin
          acknowledge, text = ['Spectrum seems to be dominated by NaNs']
          retall
        endif
        spec_w = spec_w(*, *, 0:np-1)
        spec_f = spec_f(*, *, 0:np-1)
        
        header_ref = header
        sxaddpar, header, 'EXTEND', 'T'
        sxaddpar, header, 'EXTNAME', 'flux (MJy sr-1)'
        writefits, dir+s.object_name+'_LINE'+line+'_spectra.fits', spec_f, header
        sxdelpar, header, 'EXTEND'
        header_exten = header
        sxaddpar, header_exten, 'XTENSION', 'IMAGE', before = 'SIMPLE'
        sxaddpar, header_exten, 'EXTNAME', 'wavelength (observed)'
        writefits, dir+s.object_name+'_LINE'+line+'_spectra.fits', spec_w, header_exten, /append
        header_exten = header
        sxaddpar, header_exten, 'XTENSION', 'IMAGE', before = 'SIMPLE'
        sxaddpar, header_exten, 'EXTNAME', 'wavelength (rest)'
        writefits, dir+s.object_name+'_LINE'+line+'_spectra.fits', spec_w/(1.+cube_params.redshift), header_exten, /append
        
        header = header_ref
        bigresult = [[[nresult]], [[nnoise]], [[velomap]], [[velomap_error]], [[fwhmmap]], [[fwhmmap_error]], [[contmap]], [[contmap_error]]]
        writefits, dir+s.object_name+'_LINE'+outputname+'_Flux_raw.fits', bigresult, header
        ;writefits, dir+s.object_name+'_LINE'+outputname+'_VeloMap_raw.fits', velomap, header
        ;writefits, dir+s.object_name+'_LINE'+outputname+'_FWHMMap_raw.fits', fwhmmap, header
        ;writefits, dir+s.object_name+'_LINE'+outputname+'_SpatialResolution_raw.fits', reform(grid(*, *, 5)), header
        
        nresult = PACSman_fill(nresult, header)
        nnoise = PACSman_fill(nnoise, header)
        
        bigresult = [[[nresult]], [[nnoise]], [[velomap]], [[velomap_error]], [[fwhmmap]], [[fwhmmap_error]], [[contmap]], [[contmap_error]]]
        writefits, dir+s.object_name+'_LINE'+outputname+'_Flux.fits', bigresult, header
        
        ;make DS9 region file
        pacsman_makereg2, dir+s.object_name+'_LINE'+outputname+'_Flux.fits'
        
        sxaddpar, header, 'BUNIT', 'W/m2/px'
        writefits, dir+s.object_name+'_LINE'+outputname+'_Flux_wm-2px-1.fits', [[[nresult * pixel_sr]], [[nnoise * pixel_sr]], [[velomap]], [[velomap_error]], [[fwhmmap]], [[fwhmmap_error]], [[contmap]], [[contmap_error]]], header
        
        sxaddpar, header, 'NAXIS', 2
        sxdelpar, header, 'NAXIS3'
        sxaddpar, header, 'BUNIT', 'W/m2/sr'
        angle = 5.
        hrot, nresult, header, nresult2, htmp, angle, -1, -1, 1, missing = !values.f_nan;, int = 2.
        hrot, nnoise, header, nnoise2, htmp, angle, -1, -1, 1, missing = !values.f_nan;, int = 2.
        hrot, velomap, header, velomap2, htmp, angle, -1, -1, 1, missing = !values.f_nan;, int = 2.
        hrot, velomap_error, header, velomap_error2, htmp, angle, -1, -1, 1, missing = !values.f_nan;, int = 2.
        hrot, fwhmmap, header, fwhmmap2, htmp, angle, -1, -1, 1, missing = !values.f_nan;, int = 2.
        hrot, fwhmmap_error, header, fwhmmap_error2, htmp, angle, -1, -1, 1, missing = !values.f_nan;, int = 2.
        hrot, contmap, header, contmap2, htmp, angle, -1, -1, 1, missing = !values.f_nan;, int = 2.
        hrot, contmap_error, header, contmap_error2, htmp, angle, -1, -1, 1, missing = !values.f_nan;, int = 2.
        angle = -5.
        hrot, nresult2, header, nresult2, htmp, angle, -1, -1, 1, missing = !values.f_nan;, int = 2.
        hrot, nnoise2, header, nnoise2, htmp, angle, -1, -1, 1, missing = !values.f_nan;, int = 2.
        hrot, velomap2, header, velomap2, htmp, angle, -1, -1, 1, missing = !values.f_nan;, int = 2.
        hrot, velomap_error2, header, velomap_error2, htmp, angle, -1, -1, 1, missing = !values.f_nan;, int = 2.
        hrot, fwhmmap2, header, fwhmmap2, htmp, angle, -1, -1, 1, missing = !values.f_nan;, int = 2.
        hrot, fwhmmap_error2, header, fwhmmap_error2, htmp, angle, -1, -1, 1, missing = !values.f_nan;, int = 2.
        hrot, contmap2, header, contmap2, htmp, angle, -1, -1, 1, missing = !values.f_nan;, int = 2.
        hrot, contmap_error2, header, contmap_error2, htmp, angle, -1, -1, 1, missing = !values.f_nan;, int = 2.
        ind = where( finite(nresult) and not finite(nresult2), c )
        if c gt 0 then begin
          nresult2(ind) = nresult(ind)
          nnoise2(ind) = nnoise(ind)
          velomap2(ind) = velomap(ind)
          velomap_error2(ind) = velomap_error(ind)
          fwhmmap2(ind) = fwhmmap(ind)
          fwhmmap_error2(ind) = fwhmmap_error(ind)
          contmap2(ind) = contmap(ind)
          contmap_error2(ind) = contmap_error(ind)
        endif
        nresult = nresult2
        nnoise = nnoise2
        velomap = velomap2
        velomap_error = velomap_error2
        fwhmmap = fwhmmap2
        fwhmmap_error = fwhmmap_error2
        contmap = contmap2
        contmap_error = contmap_error2
        bigresult = [[[nresult]], [[nnoise]], [[velomap]], [[velomap_error]], [[fwhmmap]], [[fwhmmap_error]], [[contmap]], [[contmap_error]]]
        sxaddpar, header, 'NAXIS', 3
        sxaddpar, header, 'NAXIS3', n_elements(bigresult(0, 0, *))
        writefits, dir+s.object_name+'_LINE'+outputname+'_Flux_smooth.fits', bigresult, header
        
        ;-------------------
        if keyword_set(all) then begin
          set_plot, 'Z'
          ERASE
          device, set_resolution=[800,800], set_pixel_depth = 24, decomposed = 0
        endif else begin
          wset, s.wset_image_map
          widget_control, s.wtab, set_tab_current = 3
        endelse
        if s.niter eq 0 then wait, 1.5
        
        map = reverse(nresult, 1)
        plotpng, map, header, mm - 1, nn - 1, s = s, filename = dir+s.object_name+'_LINE'+outputname+'_Flux.png', $
          unit_main = 'W m-2 sr-1', unit_px = 'W m-2 px-1', fact = pixel_sr
          
          
        ;-------------------
        if keyword_set(all) then begin
          set_plot, 'Z'
          ERASE
          device, set_resolution=[800,800], set_pixel_depth = 24, decomposed = 0
        endif else wset, s.wset_image_errmap
        newct, 41, ownct = s.ownct
        ;tvscl, congrid(reverse(nresult, 1), sizex, sizey)
        
        map = reverse(nnoise, 1)
        plotpng, map, header, mm - 1, nn - 1, s = s, filename = dir+s.object_name+'_LINE'+outputname+'_Error.png', $
          unit_main = 'W m-2 sr-1', unit_px = 'W m-2 px-1', fact = pixel_sr
          
          
        ;-------------------
        if keyword_set(all) then begin
          set_plot, 'Z'
          ERASE
          device, set_resolution=[800,800], set_pixel_depth = 24, decomposed = 0
        endif else wset, s.wset_image_detmap
        newct, 41, ownct = s.ownct
        ;tvscl, congrid(reverse(nresult, 1), sizex, sizey)
        
        map = reverse(nresult/nnoise, 1)
        plotpng, map, header, mm - 1, nn - 1, s = s, filename = dir+s.object_name+'_LINE'+outputname+'_Detection.png', unit_main = ''
        
        ;-------------------
        if keyword_set(all) then begin
          set_plot, 'Z'
          ERASE
          device, set_resolution=[800,800], set_pixel_depth = 24, decomposed = 0
        endif else wset, s.wset_velo_map
        ;tvscl, congrid(reverse(nresult, 1), sizex, sizey)
        
        map = reverse(velomap, 1)
        ind = where( reverse(nresult, 1) lt 3.*reverse(nnoise, 1), c)
        if c gt 0 then map(ind) = !values.f_nan
        plotpng, map, header, mm - 1, nn - 1, s = s, filename = dir+s.object_name+'_LINE'+outputname+'_Velocity.png', /velo, unit_main = '[km s'+textoidl('^{-1}')+']'
        
        ;-------------------
        if keyword_set(all) then begin
          set_plot, 'Z'
          ERASE
          device, set_resolution=[800,800], set_pixel_depth = 24, decomposed = 0
        endif else wset, s.wset_fwhm_map
        ;tvscl, congrid(reverse(nresult, 1), sizex, sizey)
        
        map = reverse(fwhmmap, 1)
        ind = where( reverse(nresult, 1) lt 3.*reverse(nnoise, 1), c)
        if c gt 0 then map(ind) = !values.f_nan
        plotpng, map, header, mm - 1, nn - 1, s = s, filename = dir+s.object_name+'_LINE'+outputname+'_FWHM.png', /fwhm, unit_main = '[km s'+textoidl('^{-1}')+']'
        
        
        ;-------------------
        if keyword_set(all) then begin
          set_plot, 'Z'
          ERASE
          device, set_resolution=[800,800], set_pixel_depth = 24, decomposed = 0
        endif else wset, s.wset_cont_map
        ;tvscl, congrid(reverse(nresult, 1), sizex, sizey)
        
        map = reverse(contmap, 1)
        plotpng, map, header, mm - 1, nn - 1, s = s, filename = dir+s.object_name+'_LINE'+outputname+'_Continuum.png', unit_main = '[Jy]'
        
        ;-------------------
        if keyword_set(all) then begin
          set_plot, 'Z'
          ERASE
          device, set_resolution=[800,800], set_pixel_depth = 24, decomposed = 0
        endif else wset, s.wset_linecont_map
        ;tvscl, congrid(reverse(nresult, 1), sizex, sizey)
        
        wave = avg(cube_lres.wave, /nan)
        contmap_um = (1.e6*contmap*1E-26*(2.9979E14/wave)/wave);contmap in MJy/sr
        contmap_error_um = (1.e6*contmap_error*1E-26*(2.9979E14/wave)/wave)
        
        map = nresult/contmap_um
        ind = where( contmap lt 3.*contmap_error or nresult lt 3.*nnoise or map lt 0., c )
        if c gt 0 then map(ind) = !values.f_nan
        map = reverse(map, 1)
        
        plotpng, map, header, mm - 1, nn - 1, s = s, filename = dir+s.object_name+'_LINE'+outputname+'_LineToContinuum.png', unit_main = ''
        
        if s.niter gt 0 then grid_ref = grid; save for next iteration
        
      endfor ;iter
      
      plot_linemap, s, all = all, ps = ps
      
    endfor ;makeall
    
    set_plot, 'X'
    
  endif else begin
    print, 'Some inputs are missing...'
  endelse
  
end
;------------------------------------------------------------------------------
pro plot_linemap, s, all = all, ps = ps
  ;###########################################

  ps = 1
  restore, s.files(s.active_file)
  
  tmp = strsplit(s.files(s.active_file), '/', /extract)
  line = tmp(n_elements(tmp)-2)

  ;dir = line + '/'
  dir = strtrans(s.files(s.active_file), tmp(n_elements(tmp)-1), '')
  cube_name = strsplit(s.files(s.active_file), '/', /extract)
  cube_name = cube_name[n_elements(cube_name)-1]
  cube_name = strtrans(cube_name, '.sav', '')

  if total(cube_params.n_rasters) eq 1 and not file_test(dir + 'footprint_results.sav') then pacsman_footprint, lines = line, nsrc = s.nsrc, cube = cube_name
    
  ;if keyword_set(all) then begin
  ;  if keyword_set(ps) then begin
  set_plot, 'PS'
  outputname = line
  if tag_exist(cube_params(0), 'pacsman_version') then if cube_params(0).pacsman_version ge 3.55 then if s.ncomp gt 0 then outputname += '_comp' + a2str(s.ncomp+1)
  device, filename = dir+s.object_name+'_LINE'+outputname+'_Line_map.eps', /color, xsize = 15, ysize = 15, /encaps
  ;  endif else begin
  ;    set_plot, 'Z'
  ;    ERASE
  ;    device, set_resolution=[800,800], set_pixel_depth = 24, decomposed = 0
  ;  endelse
  ;endif else wset, s.wset_line_map
  
  ra = cube_lres.ra
  dec = cube_lres.dec
  n = size(ra, /n_elements)
  ra = reform(ra, n)
  dec = reform(dec, n)
  
  
  spec_w = cube_spectra.lambda
  spec_f = cube_spectra.flux / ( 1.e-6 / (9.4*!pi/180./3600.)^2. );=> Jy
  spec_e = cube_spectra.error / ( 1.e-6 / (9.4*!pi/180./3600.)^2. );=> Jy
  params = cube_fitparams.params
  
  ra_center = 0.5*total(minmax(ra))
  dec_center = 0.5*total(minmax(dec))
  
  ;dra = ( max(ra) - min(ra) ) / cos( dec_center*!pi/180. )
  ;ddec = max(dec) - min(dec)
  
  if not keyword_set(steps) then steps = [20. / 3600., 20. / 3600.] ;steps in arcsec for plotting lines
  
  edgey = 5./3600.
  edgex = 5./3600. / cos(dec_center*!pi/180.)
  
  !p.position = [0.2, 0.15, 1., 1.]
  !p.region = [0.2, 0.15, 1., 1.]
  case s.use_bgimage of
    0: begin
      loadct, 0, /sil
      plot, indgen(10), /nodata, xran = [max(ra)+edgex, min(ra)-edgex], yran = [min(dec)-edgey, max(dec)+edgey], xsty = 1, ysty = 1, xtitle = 'RA J2000 (deg)', ytitle = 'DEC J2000 (deg)'
      loadct, s.ct_lines, /sil
    end
    1: begin
      outputname = line
      if tag_exist(cube_params(0), 'pacsman_version') then if cube_params(0).pacsman_version ge 3.55 then if s.ncomp gt 0 then outputname += '_comp' + a2str(s.ncomp+1)
      imag = readfits(dir+s.object_name+'_LINE'+outputname+'_Flux.fits', header, /sil)
      imag = reform(imag(*, *, 0))
      imag = reverse(imag, 1)
      sxaddpar, header, 'NAXIS', 2
      sxdelpar, header, 'NAXIS3'
      if get_equinox(header) ne 2000. then hprecess, header, 2000.
      getrot, header, rot_angle
      extast, header, astr
      ad2xy, ra_center, dec_center, astr, x, y
      hrot, imag, header, new_im, new_header, rot_angle, x, y, 1, missing = !values.f_nan
      extast, new_header, astr
      ad2xy, min(ra), min(dec), astr, x1, y1
      ad2xy, max(ra), max(dec), astr, x2, y2
      if x1 gt 0 and x2 lt n_elements(imag(*, 0)) and y1 gt 0 and y2 lt n_elements(imag(0, *)) then begin
        hextract, new_im, new_header, imag, header, x1, x2, y1, y2
        loadct, s.ct_image, /sil
        pacsman_disp, imag, xran = [max(ra)+edgex, min(ra)-edgex], yran = [min(dec)-edgey, max(dec)+edgey], xstyle=1, xtitle = 'RA J2000 (deg)', ystyle=1, ytitle = 'DEC J2000 (deg)', aspect = 1, $
          position = [0.2, 0.15, 1., 1.], charsize = 1.2, /squ;, /logscale
      endif else print, 'Cut outside image!'
    end
    2: begin
      imag = PACSman_image(s, ra_center, dec_center, header = header, status = status)
      if status eq 1 then begin
        if get_equinox(header) ne 2000. and get_equinox(header) ne 0. then hprecess, header, 2000.
        getrot, header, rot_angle
        extast, header, astr
        ad2xy, ra_center, dec_center, astr, x, y
        hrot, imag, header, new_im, new_header, rot_angle, x, y, 1, missing = !values.f_nan
        extast, new_header, astr
        ad2xy, max(ra)+edgex, min(dec)-edgey, astr, x1, y1
        ad2xy, min(ra)-edgex, max(dec)+edgey, astr, x2, y2
        if x1 gt 0 and x2 lt n_elements(imag(*, 0)) and y1 gt 0 and y2 lt n_elements(imag(0, *)) then begin
          hextract, new_im, new_header, imag, header, x1, x2, y1, y2
          if keyword_set(ps) then begin
            loadct, 0, /sil
            tmp = max(imag, /nan)-imag
          endif else begin
            loadct, s.ct_image, /sil
            tmp = imag
          endelse
          tmp_sort = tmp[reverse(sort(tmp))]
          pacsman_disp, tmp, xran = [max(ra)+edgex, min(ra)-edgex], yran = [min(dec)-edgey, max(dec)+edgey], xstyle=1, xtitle = 'RA J2000 (deg)', ystyle=1, ytitle = 'DEC J2000 (deg)', aspect = 1, $
            position = [0.2, 0.15, 1., 1.], charsize = 1.2, /squ, tickcolor, min = min(tmp_sort[0:0.999*n_elements(tmp)], /nan), max = max(tmp, /nan);, /logscale
        endif else print, 'Cut outside image!'
      endif
    end
  endcase
  
  ;plots, ra_center, dec_center, psym = 4
  loadct, s.ct_lines, /sil
  
  pa = cube_params(0).position_angle * !pi/180. ;VIAN
  rot_matrix = [ [cos(-pa), -sin(-pa)], [sin(-pa), cos(-pa)] ]
  pixel_dx = 9./3. / 3600.
  pixel_dy = 9./3. / 3600.
  bl = [-1.5*pixel_dx, -1.5*pixel_dy] # rot_matrix
  ul = [-1.5*pixel_dx, 1.5*pixel_dy] # rot_matrix
  br = [1.5*pixel_dx, -1.5*pixel_dy] # rot_matrix
  ur = [1.5*pixel_dx, 1.5*pixel_dy] # rot_matrix
  
  bl(0) /= cos(dec_center*!pi/180.)
  ul(0) /= cos(dec_center*!pi/180.)
  br(0) /= cos(dec_center*!pi/180.)
  ur(0) /= cos(dec_center*!pi/180.)
  
  for i = 0, n_elements(ra) - 1 do begin
    ;plots, ra(i), dec(i), ps = 4
  
    bl2 = bl + [ra(i), dec(i)]
    ul2 = ul + [ra(i), dec(i)]
    br2 = br + [ra(i), dec(i)]
    ur2 = ur + [ra(i), dec(i)]
    plots, [bl2(0), br2(0), ur2(0), ul2(0), bl2(0)], [bl2(1), br2(1), ur2(1), ul2(1), bl2(1)], color = 0, line = 0, thick=4
    plots, [bl2(0), br2(0), ur2(0), ul2(0), bl2(0)], [bl2(1), br2(1), ur2(1), ul2(1), bl2(1)], color = 255, line = 0
    
  endfor
  
  ;plot lines
  
  case s.line_step of
    0: arr = indgen(n) ;spaxel
    1 : begin ;raster
      arr = indgen(n)
      arr = arr(where(arr mod 25 eq 0)) + 12 ;central pixel of a raster
    end
    2: begin ;mosaic
      k = 0
      arr = fltarr(1e5)
      for rai = min(ra), max(ra), 0.1*(max(ra)-min(ra)) do begin
        for deci = min(dec), max(dec), 0.1*(max(dec)-min(dec)) do begin
          dist = sqrt( (rai-ra)^2. + (deci-dec)^2. )
          ind = sort(dist)
          arr(k) = ind(0)
          k += 1
        endfor
      endfor
      arr = arr(1:k-1)
    end
  endcase
  
  
  ;s.line_plot = 0
  case s.line_plot of
    1: begin ; plot fits
    
      scale = steps(1)/max(params(0, arr))
      for j = 0, n_elements(arr) - 1 do begin
        i = arr(j)
        w = reform(spec_w(*, i))
        ind = where(finite(w), c)
        w = w(ind)
        w = w - cube_params.lambda_rest ;median(w)
        w = -w ;for RA plot
        
        ind = where( abs(w) lt 0.5*max(w) )
        w = w(ind)
        
        y = LineModel(w, reform(params(0:2, i)))
        w /= max(w)
        w *= 4./3600. / cos(dec(i)*!pi/180.)
        y /= max(params(0, arr))
        y *= 3./3600.
        col = 255. ;* params(0, i) / max(params(0, arr))
        
        det = cube_lres(i).flux / cube_lres(i).error
        
        if det gt 2. then begin
          baseline = dec(i)-2./3600.
          oplot, ra(i)+w, baseline+y, color = 0., thick = 10
          oplot, ra(i)+w, baseline+y, color = 255., thick = 1
          xyouts, ra(i), dec(i)+max(y, /nan), a2str(trim(0.1*round(det*10.))), color = 0.
        endif
        
      endfor
    end
    0: begin ;plot data
    
      ;prepare flux array
      scale = 0.
      ;for i = 0, n_elements(spec_f(0, *))-1 do scale = max([scale, (spec_f(*, i)-median(spec_f(*, i)))])
      for i = 0, n_elements(spec_f(0, *))-1 do begin
        w = reform(spec_w(*, i))
        f = reform(spec_f(*, i))
        ind = where(finite(w) and finite(f), c)
        if c gt 0 then begin
          w = w(ind)
          f = f(ind)
          w = w - cube_params(0).lambda_rest; median(w) ;w is already rest
          w = -w ;for RA plot
          deltalambda = cube_params(0).lambda_rest * sqrt(cube_params(0).reference_fwhm^2.+cube_params(0).broadening(1)^2.) / 3.e5
          ind = where( abs(w) lt 3.*deltalambda, c )
          if c gt 0 then begin
            w = w(ind)
            f = f(ind)
            scale = max([scale, f])
          endif
        endif
      endfor
      ;help, scale
      scale *= 2.5
      
      for j = 0, n_elements(arr) - 1 do begin
        i = arr(j)
        w = reform(spec_w(*, i))
        f = reform(spec_f(*, i))
        ind = where(finite(w) and finite(f), c)
        if c gt 0 then begin
          w = w(ind)
          f = f(ind)
          w = w - cube_params(0).lambda_rest; median(w) ;w is already rest
          w = -w ;for RA plot
          
          ;I18
          deltalambda = cube_params(0).lambda_rest * sqrt(cube_params(0).reference_fwhm^2.+cube_params(0).broadening(1)^2.) / 3.e5
          ind = where( abs(w) lt 3.*deltalambda )
          w = w(ind)
          f = f(ind)
          ;print, cube_params.lambda_rest, w
          
          w /= max(w)
          w *= 4./3600. / cos(dec(i)*!pi/180.) ;.1 * (max(ra)-min(ra)) ;spaxel size
          ;f = f - median(f)
          f /= scale ;max(params(0, *))
          f *= 12./3600. ;spaxel size
          col = 255. ;* max(f) / max(params(0, *))
          
          baseline = dec(i)-2./3600.
          oplot, ra(i)+w, baseline+0.*f, color = 0., thick = 4
          oplot, ra(i)+w, baseline+0.*f, color = 255., thick = 1
          ;help, w, f
          ;f = median(f, 3)
          
          oplot, ra(i)+w, baseline+f, color = 0., thick = 10
          oplot, ra(i)+w, baseline+f, color = 255., thick = 1
        ;oplot, ra(i)+w+0.1/3600., dec(i)+f-0.1/3600., color = 255., thick = 1
        endif
      endfor
    end
  endcase
  
  loadct, 39, /sil
  l = cube_params(0).lambda_rest*(1.+cube_params(0).redshift)
  fwhm = interpol([9.2, 9.2, 9.2, 9.5, 10., 11.5, 13., 14.5], [50, 70, 80, 90, 120, 150, 180, 200], l)
  
  offd = 0.5*fwhm/3600.
  offr = offd / cos( median(cube_lres.dec) * !pi/180. )
  
  if file_test(dir + 'footprint_results.sav') then begin
    restore, dir + 'footprint_results.sav'
    ;if min(result.brightest eq result.psfbrightest) ne 0 then begin
      plots, result.psfposition(0)*[1., 1.], result.psfposition(1)+[-offd, offd], thick = 5, color = 230
      ;plots, res.position(0)*[1., 1.], res.position(1)+[-off, off], thick = 6
      plots, result.psfposition(0)+[-offr, offr], result.psfposition(1)*[1., 1.], thick = 5, color = 230
      ;plots, res.position(0)+[-off, off], res.position(1)*[1., 1.], thick = 6
      tvellipse, offr, offd, result.psfposition(0), result.psfposition(1), color = 230, /data, thick = 2

      if n_elements(result.psfposition) eq 4 then begin
      plots, result.psfposition(2)*[1., 1.], result.psfposition(3)+[-offd, offd], thick = 5, color = 230
      ;plots, res.position(0)*[1., 1.], res.position(1)+[-off, off], thick = 6
      plots, result.psfposition(2)+[-offr, offr], result.psfposition(3)*[1., 1.], thick = 5, color = 230
      ;plots, res.position(0)+[-off, off], res.position(1)*[1., 1.], thick = 6
      tvellipse, offr, offd, result.psfposition(2), result.psfposition(3), color = 230, /data, thick = 2
      endif

    ;endif else tvellipse, offr, offd, ra_center, dec_center, color = 230, /data, thick = 2
  endif else begin
    tvellipse, offr, offd, ra_center, dec_center, color = 230, /data, thick = 2
  endelse
  
  
  ;if keyword_set(all) then begin
  if keyword_set(ps) then device, /close else write_png, dir+s.object_name+'_LINE'+outputname+'_Line_map.png', tvrd(/true)
  ;endif
  ;if not keyword_set(all) then widget_control, s.line_map_text, set_value=dir+s.object_name+'_LINE'+line+'_Line_map.png'
  
  set_plot, 'X'
end
;______________________________________________________________________________
pro image_click, ev, s, grid = grid

  if keyword_set(grid) then finalpos = s.finalpos_grid else finalpos = s.finalpos
  
  if ev.press ne 1 and ev.press ne 4 then return
  ;update_plots, s, /force
  x = ev.x - finalpos(0)*s.drawsize(0)
  y = ev.y - finalpos(1)*s.drawsize(1)
  
  ra = s.ra(1) - x * (s.ra(1)-s.ra(0)) / ( (finalpos(2)-finalpos(0)) * s.drawsize(0) )
  dec = s.dec(0) + y * (s.dec(1)-s.dec(0)) / ( (finalpos(3)-finalpos(1)) * s.drawsize(1) )
  print, 'RA/DEC: ', ra, dec
  
  restore, s.files(s.active_file)
  ;cube = readfits(s.files(s.active_file), /silent)
  distance = sqrt( ((ra-cube_lres.ra)*cos(dec*!pi/180.))^2. + (dec-cube_lres.dec)^2. ) * 3600.
  ind_dist = where( distance lt 0.45*sqrt(2)*9.4, count )
  
  tmp = strsplit(s.files(s.active_file), '/', /extract)
  line = tmp(n_elements(tmp)-2)
  dir = line + '/'
  
  if count eq 0 then return ;no spaxel contributes to the sub-pixel
  
  ;ind = indsort(0)
  ;  c = 1
  ;endif
  ;if distance(indsort(0)) gt 15. then break
  
  restore, s.files(s.active_file)
  ;tmp = readfits(strtrans(s.files(s.active_file), '_res.fits', '_spectra.fits'), /silent)
  w = cube_spectra.lambda
  f = cube_spectra.flux
  e = cube_spectra.error
  
  ;finding min and max
  minw = 1.e5
  maxw = -1.e5
  minf = 1.e5
  maxf = -1.e5
  for i = 0, count - 1 do begin
    ww = cube_spectra(ind_dist(i)).lambda
    ff = cube_spectra(ind_dist(i)).flux
    minw = min([minw, min(ww, /nan)], /nan)
    maxw = max([maxw, max(ww, /nan)], /nan)
    minf = min([minf, min(ff, /nan)], /nan)
    maxf = max([maxf, max(ff, /nan)], /nan)
  endfor
  spawn, '\ls '+line+'/*_spectra.fits', file ;this is for the range
  if file(0) ne '' then begin
    spec_f = readfits(file(0), h, /sil)
  endif else spec_f = 0.
  minf = min([minf, min(spec_f, /nan)], /nan)
  maxf = max([maxf, max(spec_f, /nan)], /nan)
  
  xtit = textoidl('\lambda-\lambda_{rest} [\mum]')
  if ev.press eq 1 then begin
    window, 0
    loadct, 0, /sil
    plot, indgen(10), xran = [minw, maxw], yran = [minf, maxf], /nodata, title = 'Spectra corresponding to position: ' + a2str(ra) + ' / ' + a2str(dec), $
      xtit = xtit, ytit = 'Flux density [MJy sr-1]', xsty = 1
    plots, [0, 0], !y.crange, linestyle = 2, color = 100
    plots, !x.crange, [0, 0], linestyle = 0, color = 100
    loadct, 38, /sil
  endif
  legend_str = strarr(count)
  
  for i = 0, count - 1 do begin
    if ev.press eq 4 then spawn, 'open ./'+line+'/PLOTS_*/raster'+a2str(cube_spectra(ind_dist(i)).raster)+'_'+a2str(cube_spectra(ind_dist(i)).spaxel_X)+a2str(cube_spectra(ind_dist(i)).spaxel_Y)+'.eps'
    
    legend_str(i) = 'Raster ' + a2str(cube_spectra(ind_dist(i)).raster) + ', spaxel ' + a2str(cube_spectra(ind_dist(i)).spaxel_X) + ' ' + a2str(cube_spectra(ind_dist(i)).spaxel_Y)
    
    ww = cube_spectra(ind_dist(i)).lambda
    ff = cube_spectra(ind_dist(i)).flux
    ind2 = where( finite(ww) and finite(ff) )
    if ev.press eq 1 then oplot, ww(ind2), ff(ind2), ps = 10, color = 50 + 50 * i
  endfor
  if ev.press eq 1 then begin
    spawn, '\ls '+line+'/*_spectra.fits', file
    if file(0) ne '' then begin
      spec_f = readfits(file(0), h, /sil)
      spec_w = readfits(file(0), h, /sil, exten_no=2)
      extast, h, astr
      ad2xy, ra, dec, astr, xx, yy
      xx = nint(xx)
      yy = nint(yy)
      loadct, 0, /sil
      oplot, spec_w(xx, yy, *), spec_f(xx, yy, *), ps = 10, color = 250
      al_legend, 'Rebinned spectrum', color = 250, textcolor = 250, box = 0, /top, /left
    endif
    loadct, 38, /sil
    al_legend, legend_str, color = 50 + 50*indgen(c), textcolor = 50 + 50*indgen(c), box = 0, /top, /right
  endif
  
  
end
;______________________________________________________________________________
;______________________________________________________________________________
pro PACSman_map_event, ev
  ;event handler
  widget_control, ev.id, get_uvalue=uvalue
  widget_control, ev.top, get_uvalue = s
  CASE uvalue OF
    'none':
    'tab_event':
    'image_subpixels2': begin
      image_click, ev, s, /grid
    end
    'image_map': begin
      image_click, ev, s
    end
    'velo_map': begin
      widget_control, s.image_map, send_event=ev
    end
    'fwhm_map': begin
      widget_control, s.image_map, send_event=ev
    end
    'posangle_change': begin
      widget_control, s.posangle_input, get_value=tmp
      s.position_angle = tmp
      update_plots, s
    end
    'file_list': begin
      s.active_file = ev.index - 1
      print, 'S.ACTIVE_FILE'
      help, s.active_file
    ;s.band = ''
    ;widget_control, s.band_list, set_combobox_select = 0
    ;s.bgimage = ''
    ;widget_control, s.bgimage_list, set_combobox_select = 0
    ;widget_control, s.posangle_text, set_value = ''
    end
    'niter_list': begin
      case ev.index of
        0: s.niter = 0
        1: s.niter = 1
        2: s.niter = 4
        3: s.niter = 6
        4: s.niter = 10
        5: s.niter = 20
      endcase
    end
    ;'band_list': begin
    ;  case ev.index of
    ;    0: s.band = ''
    ;    1: s.band = 'blue'
    ;    2: s.band = 'red'
    ;  endcase
    ;end
    ;'update': begin
    ;  update_plots, s, /force
    ;end
    'make': begin
      update_plots, s
    end
    'makeall': begin
      update_plots, s, /all
    end
    'object_name': begin
      widget_control, s.object_text, get_value=tmp
      tmp = strtrans(tmp, ' ', '_')
      tmp = strtrans(tmp, '\\', '_')
      tmp = strtrans(tmp, '\/', '_')
      tmp = strtrans(tmp, '%', '_')
      widget_control, s.object_text, set_value=tmp
      s.object_name = tmp
      print, 'Object name changed to ' + tmp
    end
    'autofind_object_ned': begin
      restore, s.files(s.active_file)
      ra = cube_lres.ra
      dec = cube_lres.dec
      ra_center = 0.5*total(minmax(ra))
      dec_center = 0.5*total(minmax(dec))
      sss = "http://ned.ipac.caltech.edu/cgi-bin/nph-objsearch?search_type=Near+Position+Search&of=xml_main&"
      sss = sss + "ra=" + a2str(ra_center) + "&dec=" + a2str(dec_center) + "&sr=" + a2str(0.015)
      query = webget(sss)
      ind_data = where(strmid(query.text, 0, 11) eq '<TABLEDATA>')
      if ind_data[0] eq -1 then begin
        print, 'no sources!'
      endif else begin
        objects_str = strarr(1)
        r = ind_data + 1
        while strmid(query.text[r], 0, 4) eq '<TR>' do begin
          st = query.text[r+1]
          num = nint( strmid(st, 4, strpos(st, '</TD>') - 4) )
          st = query.text[r+2]
          name = strmid(st, 4, strpos(st, '</TD>') - 4)
          st = query.text[r+3]
          ra = double( strmid(st, 4, strpos(st, '</TD>') - 4) )
          st = query.text[r+4]
          dec = double( strmid(st, 4, strpos(st, '</TD>') - 4) )
          st = query.text[r+5]
          type = strmid(st, 4, strpos(st, '</TD>') - 4)
          ;objects_str = [objects_str, 'NED ' + a2str(num) + ' | ' + a2str(ra) + ' | ' + a2str(dec) + ' | ' + name + ' | ' + type]
          objects_str = [objects_str, name]
          r = r + 19
        endwhile
        objects_str = objects_str[1:n_elements(objects_str)-1]
        print, objects_str
        tmp = objects_str(0)
        tmp = strtrans(tmp, ' ', '_')
        tmp = strtrans(tmp, '\\', '_')
        tmp = strtrans(tmp, '\/', '_')
        tmp = strtrans(tmp, '%', '_')
        widget_control, s.object_text, set_value=tmp
        s.object_name = tmp
      endelse
    end
    'autofind_object_simbad': begin
      restore, s.files(s.active_file)
      ra = cube_lres.ra
      dec = cube_lres.dec
      ra_center = 0.5*total(minmax(ra))
      dec_center = 0.5*total(minmax(dec))
      tmp = sm_simbad_getname(ra_center, dec_center)
      if tmp eq '' then begin
        PACSman_acknowledge, 'Object not found...'
      endif else begin
        tmp = strtrans(tmp, ' ', '_')
        tmp = strtrans(tmp, '\\', '_')
        tmp = strtrans(tmp, '\/', '_')
        tmp = strtrans(tmp, '%', '_')
        widget_control, s.object_text, set_value=tmp
        s.object_name = tmp
      endelse
    end
    'display_res_change': begin
      s.display_resolution = abs(s.display_resolution - 1)
      update_plots, s
    end
    'bgimage_list': begin
      case ev.index of
        0: s.bgimage = ''
        1: begin
          file = dialog_pickfile(/must_exist, filter=['*.fits', '*'], title='Choose an image to read')
          if file eq '' then break
          s.bgimage = file
        end
        2: s.bgimage = 'DSS/2b'
        3: s.bgimage = 'DSS/2r'
        4: s.bgimage = 'DSS/1'
        5: s.bgimage = 'DSS/2i'
        6: s.bgimage = '2MASS/J'
        7: s.bgimage = '2MASS/H'
        8: s.bgimage = '2MASS/K'
        9: s.bgimage = 'IRAC/ch1'
        10: s.bgimage = 'IRAC/ch2'
        11: s.bgimage = 'IRAC/ch3'
        12: s.bgimage = 'IRAC/ch4'
        13: s.bgimage = 'MIPS/24um'
        14: s.bgimage = 'MIPS/70um'
        15: s.bgimage = 'MIPS/160um'
      endcase
      widget_control, s.bgimage_text, set_value = 'Background image selected: ' + s.bgimage
    end
    'calib_list': begin
      case ev.index of
        0: s.use_calib_error = 0
        1: s.use_calib_error = 10
        2: s.use_calib_error = 20
      endcase
      print, 'Using calibration error [%]: ' + a2str(s.use_calib_error)
    end
    'comp_change': begin
      widget_control, s.comp_text, get_value = tmp
      s.ncomp = tmp
      print, 'changed component number to:', s.ncomp
    end
    'line_map_proceed': begin
      plot_linemap, s
    end
    'line_plot_choice': begin
      s.line_plot = ev.index
    end
    'line_step_choice': begin
      s.line_step = ev.index
    end
    'use_bgimage_nothing': begin
      s.use_bgimage = 0
    end
    'use_bgimage_map': begin
      s.use_bgimage = 1
    end
    'use_bgimage_image': begin
      s.use_bgimage = 2
    end
    'ct_image': begin
      widget_control, s.ct_image_box, get_value = tmp
      if tmp lt 0 then tmp = 0
      if tmp gt 39 then tmp = 39
      s.ct_image = tmp
      widget_control, s.ct_image_box, set_value = a2str(tmp)
      print, 'Color table changed to ' + a2str(tmp)
    end
    'ct_lines': begin
      widget_control, s.ct_lines_box, get_value = tmp
      if tmp lt 0 then tmp = 0
      if tmp gt 39 then tmp = 39
      s.ct_lines = tmp
      widget_control, s.ct_lines_box, set_value = a2str(tmp)
      print, 'Color table changed to ' + a2str(tmp)
    end
    ;'use_correctionfactor': begin
    ;  s.use_correctionfactor = abs( s.use_correctionfactor - 1 )
    ;  print, s.use_correctionfactor
    ;end
    'quit': widget_control, ev.top, /destroy
  ENDCASE
  if uvalue ne 'quit' then widget_control, ev.top, set_uvalue = s
end
;______________________________________________________________________________
;______________________________________________________________________________
;MAIN PROGRAM
pro PACSman_map, load_image = load_image, noedge = noedge, object = object, makeall = makeall, integrated = integrated, ownct = ownct, noct = noct, $
    help = help, ps = ps, use_calib_error = use_calib_error, perturb = perturb, ncomp = ncomp, nsrc = nsrc;, script = script,
    
  ;set_plot, 'X'
  ;cleanplot, /sil
    
  if keyword_set(help) then begin
    print, 'Calling sequence.'
    print, ''
    print, "Minimum calling sequence example: "
    print, "  pacsman_map              {brings up the GUI}"
    print, ""
    print, "Return FITS files and PNG snapshots. FITS planes:"
    print, " flux, error on flux, radial velocity, error on radial velocity, line FWHM, error on FWHM, continuum, error on continuum"
    print, ""
    print, "Optional :"
    print, "  load_image             = 'DSS/2b'            {Load background image. Choices: DSS/2b,2r,1,2i - 2MASS/J,H,K - IRAC/ch1,ch2,ch3,ch4 - MIPS/24um,70um,160um}"
    print, "  object                 = 'xxxxxxxx'          {Object name for output files}"
    print, "  ncomp                  = 0                   {what component to project [default = 0, i.e., nominal line]}"
    print, ""
    print, "Optional switches:"
    print, "  /makeall               {builds all the maps that can be made from the .sav files from pacsman_fit}"
    print, "  /ownct                 {color table will be saved and used in the /tmp folder}"
    print, "  /noct                  {no color table}"
    print, "  /use_calib_error       {adding a calibration error in % [default = 10]}"
    print, "  /perturb               {perturbing flux map for MC error calculation}"
    print, ""
    retall
  endif
  
  !p.multi = 0
  !p.position = 0
  
  if not keyword_set(nsrc) then nsrc = 1

  if not keyword_set(makeall) then begin
    device,true_color=24
    device,decomposed=0
    ;device,retain=2
  endif
  
  if keyword_set(use_calib_error) then begin
    if use_calib_error eq 1 then use_calib_error = 10.
    print, 'Using calibration error [%]: ' + a2str(use_calib_error)
  endif else use_calib_error = 0.
  
  if not keyword_set(ncomp) then ncomp = 0
  
  print, 'Searching for cubes...'
  files = file_search('./', 'cube*.sav', count = count_files)
  if count_files eq 0 then begin
    print, 'No valid files found!'
    retall
  endif
  ind = where( strpos(files, 'datacloud') eq -1 and strpos(files, 'Grid') eq -1 and strpos(files, '_3x3') eq -1 and strpos(files, '_5x5') eq -1)
  files = files(ind)
  
  if not keyword_set(object) then object = 'OBJECT'
  
  ;##################################################################
  
  if keyword_set(makeall) then begin
    sc_size = [900, 900]
    sc_res = [2.6e-2, 2.6e-2]
  endif else sc_size = get_screen_size(resolution=sc_res)
  x_sc = min([0.95*sc_size(0), 1024])
  x_res = sc_res(0)
  y_sc = min([0.95*sc_size(1), 768])
  y_res = sc_res(1)
  
  s = { $
    wset_image_footprints: 0, $
    wset_image_subpixels: 0, $
    wset_image_subpixels2: 0, $
    wset_image_map: 0, $
    wset_image_detmap: 0, $
    wset_image_errmap: 0, $
    wset_velo_map: 0, $
    wset_fwhm_map: 0, $
    wset_cont_map: 0, $
      wset_linecont_map: 0, $
      integrated: keyword_set(integrated), $
    ;image_subpixels_text: image_subpixels_text, $
    ;image_subpixels2_text: image_subpixels2_text, $
    ;velo_map_text: velo_map_text, $
    ;fwhm_map_text: fwhm_map_text, $
    wset_line_map: 0, $
    image_map: 0, $
    image_detmap: 0, $
    velo_map: 0, $
    fwhm_map: 0, $
    cont_map: 0, $
    linecont_map: 0, $
    comp_text: 0, $
    line_plot: 0, $
    wtab: 0, $
    ;use_correctionfactor_box: use_correctionfactor_box, $
    ;use_correctionfactor: 1, $
    line_step: 0, $
    line_map_text: 0, $
    ;velo_map_text2: velo_map_text2, $
    ;fwhm_map_text2: fwhm_map_text2, $
    ;image_map_text: image_map_text, $
    ;image_map_text2: image_map_text2, $
    ;image_detmap_text: image_detmap_text, $
    ;image_detmap_text2: image_detmap_text2, $
    ;image_footprints_text: image_footprints_text, $
    files: files, $
    noedge: keyword_set(noedge), $
    ;object: object, $
    display_resolution: 0, $
    cube_name: '', $
    band: '', $
    active_file: -1, $
    groupid: 0, $
    ct_lines: 0, $
    ct_image: 39, $
    ct_lines_box: 0, $
    ct_image_box: 0, $
    drawsize: [800, 800], $
    bgimage: '', $
    btn_update: 0, $
    unit: 'wm-2sr-1', $
    bgimage_list: 0, $
    bgimage_text: 0, $
    calib_list: 0, $
    ;band_list: band_list, $
    ra: [0., 0.], $
    dec: [0., 0.], $
    ra_center: 0.d, $
    dec_center: 0.d, $
    dx: [0., 0.], $
    dy: [0., 0.], $
    object_name: object, $
    object_text: 0, $
    band_text: 0, $
    line_text: 0, $
    posangle_value: 0., $
    use_bgimage: 0, $
    ownct: keyword_set(ownct), $
    niter: 0, $
    finalpos: fltarr(4.), $
    finalpos_grid: fltarr(4.), $
    use_calib_error: use_calib_error, $
    ncomp: ncomp, $
    nsrc: nsrc, $
    perturb: keyword_set(perturb) $
    }
    
  if keyword_set(load_image) then begin
    s.use_bgimage = 2
    s.bgimage = load_image ;in case this is a fits file
  endif
  
  if not keyword_set(makeall) then begin
  
    main = widget_base(title='PACSman_map', /column, xsize=x_sc, ysize=y_sc, /base_align_top) ;main window
    wTab = widget_tab(main, uvalue='tab_event') ;tabs
    
    
    ;------------------------------
    
    wT1 = WIDGET_BASE(wTab, TITLE=' Parameters & Footprints ', /COLUMN)
    
    main_panel1 = widget_base (wT1, /row, scr_xsize=0.99*x_sc, scr_ysize=0.99*y_sc)
    
    left_panel1 = widget_base(main_panel1, /column, scr_xsize=0.3*x_sc, scr_ysize=0.98*y_sc)
    lbl = widget_label(left_panel1, value='Input parameters', /align_left)
    lbl = widget_label(left_panel1, value=' ', /align_left)
    file_list = widget_combobox (left_panel1, uvalue='file_list', value=['Choose map...', files], xsize=100)
    ;band_list = widget_combobox (left_panel1, uvalue='band_list', value=['Choose band...', 'blue', 'red'], xsize=100)
    lbl = widget_label(left_panel1, value=' ', /align_left)
    
    tmp12 = widget_base(left_panel1, /column, /frame)
    lbl = widget_label(tmp12, value='Line:', /align_left)
    line_text = widget_label(tmp12, value='...                                       ', /align_left)
    
    tmp11 = widget_base(left_panel1, /column, /frame)
    lbl = widget_label(tmp11, value='Band:', /align_left)
    band_text = widget_label(tmp11, value='...                                       ', /align_left)
    
    ;tmp1 = widget_base(left_panel1, /column, /frame)
    ;lbl = widget_label(tmp1, value='Position angle:', /align_left)
    ;posangle_text = widget_label(tmp1, value='...                                       ', /align_left)
    
    ;tmp13 = widget_base(left_panel1, /column, /frame)
    ;lbl = widget_label(tmp13, value='Use flux calibration correction factor:', /align_left)
    ;lbl = widget_label(tmp13, value='> calTree = getCalTree(obs=obs)', /align_left)
    ;lbl = widget_label(tmp13, value='> print calTree.spectrometer', /align_left)
    ;lbl = widget_label(tmp13, value='If time dependency < 13 the factor is needed', /align_left)
    ;use_factor_base = widget_base(tmp13, /nonexclusive, column=1)
    ;use_correctionfactor_box = widget_button(use_factor_base, value='Use correction factor', uvalue='use_correctionfactor')
    
    lbl = widget_label(left_panel1, value=' ', /align_left)
    
    tmp2 = widget_base(left_panel1, /column, /frame)
    
    lbl = widget_label(tmp2, value='Optional: choose background image',  /align_left)
    bgimage_list = widget_combobox(tmp2, uvalue='bgimage_list', value=['Choose image...', 'Browse for a FITS file', 'DSS/2b', 'DSS/2r', 'DSS/1', 'DSS/2i', '2MASS/J', '2MASS/H', '2MASS/K', $
      'IRAC/ch1', 'IRAC/ch2', 'IRAC/ch3', 'IRAC/ch4', 'MIPS/24um', 'MIPS/70um', 'MIPS/160um'], xsize=100)
    bgimage_text = widget_text(tmp2, value='...', /align_left, xsize=40)
    lbl = widget_label(tmp2, value=' ', /align_left)
    
    lbl = widget_label(tmp2, value='Optional: object name (hit enter afterwards):', /align_left)
    tmp0 = widget_base(tmp2, /row)
    object_text = widget_text(tmp0, uvalue='object_name', value=object, /editable, xsize=23)
    lbl = widget_label(tmp0, value=', or ', /align_left)
    btn = widget_button(tmp0, uvalue='autofind_object_ned', value='NED', ysize=10)
    btn = widget_button(tmp0, uvalue='autofind_object_simbad', value='SIMBAD', ysize=10)
    lbl = widget_label(tmp2, value=' ', /align_left)
    
    lbl = widget_label(tmp2, value='Add calibration error?', /align_left)
    calib_list = widget_combobox(tmp2, uvalue='calib_list', value=['no', '10%', '20%'], xsize=40)
    
    lbl = widget_label(tmp2, value='Component number (press enter to validate)', /align_left)
    comp_text = widget_text(tmp2, uvalue='comp_change', value='0', /editable, xsize=5)
    
    ;lbl = widget_label(tmp2, value='Beta: niter',  /align_left)
    ;niter_list = widget_combobox(tmp2, uvalue='niter_list', value=['0', '1', '2'], xsize=100)
    
    lbl = widget_label(left_panel1, value=' ', /align_left)
    btn = widget_button(left_panel1, uvalue='make', value=' Make map ', xsize=0.99*x_sc, ysize=30)
    btn = widget_button(left_panel1, uvalue='makeall', value=' Make all maps ', xsize=0.99*x_sc, ysize=30)
    ;btn_update = widget_button(left_panel1, uvalue='update', value=' Update maps ', xsize=0.99*x_sc, ysize=30, sensitive=0)
    lbl = widget_label(left_panel1, value=' ', /align_left)
    btn_quit = widget_button(left_panel1, uvalue='quit', value=' Quit ', xsize=0.99*x_sc, ysize=30)
    
    right_panel1 = widget_base(main_panel1, /column, scr_xsize=0.7*x_sc, scr_ysize=0.98*y_sc)
    image_footprints = widget_draw (right_panel1, uvalue='image_footprints', xsize=0.65*x_sc, ysize=0.85*y_sc)
    ;tmp = widget_label(right_panel1, value='Image is saved at:', /align_left)
    ;image_footprints_text = widget_text(right_panel1, value=' ', /align_left, xsize=100)
    
    
    ;------------------------------
    
    wT2 = WIDGET_BASE(wTab, TITLE=' Projected grid ', /COLUMN)
    
    main_panel2 = widget_base (wT2, /row, xsize=0.99*x_sc, ysize=0.99*y_sc)
    
    left_panel2 = widget_base(main_panel2, /column, xsize=0.3*x_sc)
    ;lbl = widget_label(left_panel2, value='Display spatial sampling (beta)', /align_left)
    ;display_res_base = widget_base(left_panel2, /nonexclusive, column=1)
    ;display_res = widget_button(display_res_base, value='', uvalue='display_res_change')
    
    right_panel2 = widget_base(main_panel2, /column, xsize=0.7*x_sc, ysize=0.99*y_sc)
    tmp = widget_label(right_panel2, value='Note: subpixels are 9.4/3 arcsec size. Squares within indicate no spaxel contributes to the supixel.', /align_left)
    image_subpixels = widget_draw (right_panel2, uvalue='image_subpixels', xsize=0.65*x_sc, ysize=0.85*y_sc)
    ;tmp = widget_label(right_panel2, value='Image is saved at:', /align_left)
    ;image_subpixels_text = widget_text(right_panel2, value=' ', /align_left, xsize=100)
    
    
    ;------------------------------
    
    wT22 = WIDGET_BASE(wTab, TITLE=' Projected grid + Flux ', /COLUMN)
    
    main_panel22 = widget_base (wT22, /row, xsize=0.99*x_sc, ysize=0.99*y_sc)
    
    left_panel22 = widget_base(main_panel22, /column, xsize=0.3*x_sc)
    lbl = widget_label(left_panel22, value='left click: show spectra', /align_left)
    lbl = widget_label(left_panel22, value='right click: show spectra and open fits', /align_left)
    ;display_res_base = widget_base(left_panel22, /nonexclusive, column=1)
    ;display_res = widget_button(display_res_base, value='', uvalue='display_res_change')
    
    right_panel22 = widget_base(main_panel22, /column, xsize=0.7*x_sc, ysize=0.99*y_sc)
    tmp = widget_label(right_panel22, value='Note: subpixels are 9.4/3 arcsec size.', /align_left)
    image_subpixels2 = widget_draw (right_panel22, uvalue='image_subpixels2', xsize=0.65*x_sc, ysize=0.85*y_sc, /button)
    ;tmp = widget_label(right_panel22, value='Image is saved at:', /align_left)
    ;image_subpixels2_text = widget_text(right_panel22, value=' ', /align_left, xsize=100)
    
    
    
    ;------------------------------
    
    wT3 = WIDGET_BASE(wTab, TITLE=' Flux ', /COLUMN)
    
    main_panel3 = widget_base (wT3, /row, xsize=0.99*x_sc, ysize=0.99*y_sc)
    
    
    left_panel3 = widget_base(main_panel3, /column, xsize=0.3*x_sc)
    lbl = widget_label(left_panel3, value='Flux in [W m-2 px-1] (top), [W m-2 sr-1] (bottom)', /align_left)
    lbl = widget_label(left_panel3, value='with 9.4"/3 pixel size', /align_left)
    lbl = widget_label(left_panel3, value='', /align_left)
    lbl = widget_label(left_panel3, value='left click: show spectra', /align_left)
    lbl = widget_label(left_panel3, value='right click: show spectra and open fits', /align_left)
    ;display_res_base = widget_base(left_panel22, /nonexclusive, column=1)
    ;display_res = widget_button(display_res_base, value='', uvalue='display_res_change')
    
    right_panel3 = widget_base(main_panel3, /column, xsize=0.7*x_sc, ysize=0.99*y_sc)
    tmp = widget_label(right_panel3, value='', /align_left)
    image_map = widget_draw (right_panel3, uvalue='image_map', xsize=0.65*x_sc, ysize=0.85*y_sc, /button)
    ;tmp = widget_label(right_panel3, value='Image is saved at:', /align_left)
    ;image_map_text = widget_text(right_panel3, value=' ', /align_left, xsize=100)
    ;image_map_text2 = widget_text(right_panel3, value=' ', /align_left, xsize=100)
    
    
    ;------------------------------
    
    wT31 = WIDGET_BASE(wTab, TITLE=' Error ', /COLUMN)
    
    main_panel31 = widget_base (wT31, /row, xsize=0.99*x_sc, ysize=0.99*y_sc)
    
    left_panel31 = widget_base(main_panel31, /column, xsize=0.3*x_sc)
    lbl = widget_label(left_panel31, value='Flux in [W m-2 px-1] (top), [W m-2 sr-1] (bottom)', /align_left)
    lbl = widget_label(left_panel31, value='with 9.4"/3 pixel size', /align_left)
    lbl = widget_label(left_panel31, value='', /align_left)
    lbl = widget_label(left_panel31, value='left click: show spectra', /align_left)
    lbl = widget_label(left_panel31, value='right click: show spectra and open fits', /align_left)
    ;display_res_base = widget_base(left_panel22, /nonexclusive, column=1)
    ;display_res = widget_button(display_res_base, value='', uvalue='display_res_change')
    
    right_panel31 = widget_base(main_panel31, /column, xsize=0.7*x_sc, ysize=0.99*y_sc)
    tmp = widget_label(right_panel31, value='', /align_left)
    image_errmap = widget_draw (right_panel31, uvalue='image_map', xsize=0.65*x_sc, ysize=0.85*y_sc, /button)
    ;tmp = widget_label(right_panel31, value='Image is saved at:', /align_left)
    ;image_detmap_text = widget_text(right_panel31, value=' ', /align_left, xsize=100)
    ;image_detmap_text2 = widget_text(right_panel31, value=' ', /align_left, xsize=100)
    
    
    ;------------------------------
    
    wT32 = WIDGET_BASE(wTab, TITLE=' Detection ', /COLUMN)
    
    main_panel32 = widget_base (wT32, /row, xsize=0.99*x_sc, ysize=0.99*y_sc)
    
    left_panel32 = widget_base(main_panel32, /column, xsize=0.3*x_sc)
    lbl = widget_label(left_panel32, value='Detecion level in [sigma]', /align_left)
    lbl = widget_label(left_panel32, value='', /align_left)
    lbl = widget_label(left_panel32, value='left click: show spectra', /align_left)
    lbl = widget_label(left_panel32, value='right click: show spectra and open fits', /align_left)
    ;display_res_base = widget_base(left_panel22, /nonexclusive, column=1)
    ;display_res = widget_button(display_res_base, value='', uvalue='display_res_change')
    
    right_panel32 = widget_base(main_panel32, /column, xsize=0.7*x_sc, ysize=0.99*y_sc)
    tmp = widget_label(right_panel32, value='', /align_left)
    image_detmap = widget_draw (right_panel32, uvalue='image_map', xsize=0.65*x_sc, ysize=0.85*y_sc, /button)
    ;tmp = widget_label(right_panel31, value='Image is saved at:', /align_left)
    ;image_detmap_text = widget_text(right_panel31, value=' ', /align_left, xsize=100)
    ;image_detmap_text2 = widget_text(right_panel31, value=' ', /align_left, xsize=100)
    
    
    
    ;------------------------------
    
    wT35 = WIDGET_BASE(wTab, TITLE=' Radial velocity ', /COLUMN)
    
    
    main_panel35 = widget_base (wT35, /row, xsize=0.99*x_sc, ysize=0.99*y_sc)
    
    left_panel35 = widget_base(main_panel35, /column, xsize=0.3*x_sc)
    lbl = widget_label(left_panel35, value='Radial velocity in [km s-1]', /align_left)
    lbl = widget_label(left_panel35, value='', /align_left)
    lbl = widget_label(left_panel35, value='left click: show spectra', /align_left)
    lbl = widget_label(left_panel35, value='right click: show spectra and open fits', /align_left)
    ;display_res_base = widget_base(left_panel22, /nonexclusive, column=1)
    ;display_res = widget_button(display_res_base, value='', uvalue='display_res_change')
    
    right_panel35 = widget_base(main_panel35, /column, xsize=0.7*x_sc, ysize=0.99*y_sc)
    tmp = widget_label(right_panel35, value='', /align_left)
    velo_map = widget_draw (right_panel35, uvalue='velo_map', xsize=0.65*x_sc, ysize=0.85*y_sc, /button)
    ;tmp = widget_label(right_panel35, value='Image is saved at:', /align_left)
    ;velo_map_text = widget_text(right_panel35, value=' ', /align_left, xsize=100)
    
    
    ;------------------------------
    ;------------------------------
    
    wT36 = WIDGET_BASE(wTab, TITLE=' FWHM_int ', /COLUMN)
    
    
    main_panel36 = widget_base (wT36, /row, xsize=0.99*x_sc, ysize=0.99*y_sc)
    
    
    left_panel36 = widget_base(main_panel36, /column, xsize=0.3*x_sc)
    lbl = widget_label(left_panel36, value='FWHM in [km s-1]', /align_left)
    lbl = widget_label(left_panel36, value='', /align_left)
    lbl = widget_label(left_panel36, value='left click: show spectra', /align_left)
    lbl = widget_label(left_panel36, value='right click: show spectra and open fits', /align_left)
    ;display_res_base = widget_base(left_panel22, /nonexclusive, column=1)
    ;display_res = widget_button(display_res_base, value='', uvalue='display_res_change')
    
    right_panel36 = widget_base(main_panel36, /column, xsize=0.7*x_sc, ysize=0.99*y_sc)
    tmp = widget_label(right_panel36, value='', /align_left)
    fwhm_map = widget_draw (right_panel36, uvalue='image_map', xsize=0.65*x_sc, ysize=0.85*y_sc, /button)
    ;tmp = widget_label(right_panel36, value='Image is saved at:', /align_left)
    ;fwhm_map_text = widget_text(right_panel36, value=' ', /align_left, xsize=100)
    ;fwhm_map_text2 = widget_text(right_panel36, value=' ', /align_left, xsize=100)
    
    ;------------------------------
    ;------------------------------
    
    wT37 = WIDGET_BASE(wTab, TITLE=' Cont ', /COLUMN)
    
    
    main_panel37 = widget_base (wT37, /row, xsize=0.99*x_sc, ysize=0.99*y_sc)
    
    
    left_panel37 = widget_base(main_panel37, /column, xsize=0.3*x_sc)
    ;lbl = widget_label(left_panel37, value='FWHM in [km s-1]', /align_left)
    lbl = widget_label(left_panel37, value='', /align_left)
    lbl = widget_label(left_panel37, value='left click: show spectra', /align_left)
    lbl = widget_label(left_panel37, value='right click: show spectra and open fits', /align_left)
    ;display_res_base = widget_base(left_panel22, /nonexclusive, column=1)
    ;display_res = widget_button(display_res_base, value='', uvalue='display_res_change')
    
    right_panel37 = widget_base(main_panel37, /column, xsize=0.7*x_sc, ysize=0.99*y_sc)
    tmp = widget_label(right_panel37, value='', /align_left)
    cont_map = widget_draw (right_panel37, uvalue='image_map', xsize=0.65*x_sc, ysize=0.85*y_sc, /button)
    ;tmp = widget_label(right_panel36, value='Image is saved at:', /align_left)
    ;fwhm_map_text = widget_text(right_panel36, value=' ', /align_left, xsize=100)
    ;fwhm_map_text2 = widget_text(right_panel36, value=' ', /align_left, xsize=100)
    
    
    ;------------------------------
    ;------------------------------
    
    wT38 = WIDGET_BASE(wTab, TITLE=' Line/Cont ', /COLUMN)
    
    
    main_panel38 = widget_base (wT38, /row, xsize=0.99*x_sc, ysize=0.99*y_sc)
    
    
    left_panel38 = widget_base(main_panel38, /column, xsize=0.3*x_sc)
    ;lbl = widget_label(left_panel37, value='FWHM in [km s-1]', /align_left)
    lbl = widget_label(left_panel38, value='', /align_left)
    lbl = widget_label(left_panel38, value='left click: show spectra', /align_left)
    lbl = widget_label(left_panel38, value='right click: show spectra and open fits', /align_left)
    ;display_res_base = widget_base(left_panel22, /nonexclusive, column=1)
    ;display_res = widget_button(display_res_base, value='', uvalue='display_res_change')
    
    right_panel38 = widget_base(main_panel38, /column, xsize=0.7*x_sc, ysize=0.99*y_sc)
    tmp = widget_label(right_panel38, value='', /align_left)
    linecont_map = widget_draw (right_panel38, uvalue='image_map', xsize=0.65*x_sc, ysize=0.85*y_sc, /button)
    ;tmp = widget_label(right_panel36, value='Image is saved at:', /align_left)
    ;fwhm_map_text = widget_text(right_panel36, value=' ', /align_left, xsize=100)
    ;fwhm_map_text2 = widget_text(right_panel36, value=' ', /align_left, xsize=100)
    
    
    ;------------------------------
    
    
    wT5 = WIDGET_BASE(wTab, TITLE=' Line map ', /COLUMN)
    
    
    main_panel5 = widget_base (wT5, /row, xsize=0.99*x_sc, ysize=0.99*y_sc)
    
    
    left_panel5 = widget_base(main_panel5, /column, xsize=0.3*x_sc)
    
    tmp = widget_label(left_panel5, value='What to plot:', /align_left)
    line_plot_cb = widget_combobox (left_panel5, uvalue='line_plot_choice', value=[' data ', ' fits '], xsize=100)
    tmp = widget_label(left_panel5, value='', /align_left)
    tmp = widget_label(left_panel5, value='Spatial step:', /align_left)
    line_step_cb = widget_combobox (left_panel5, uvalue='line_step_choice', value=[' each spaxel ', ' central spaxels ', ' mosaic '], xsize=100)
    btn_linemap = widget_button(left_panel5, uvalue='line_map_proceed', value=' Create EPS file ', xsize=0.99*x_sc, ysize=30)
    tmp = widget_label(left_panel5, value='', /align_left)
    tmp = widget_label(left_panel5, value='', /align_left)
    tmp = widget_label(left_panel5, value='Background: ', /align_left)
    base_btn = widget_base(left_panel5, /exclusive)
    btn1 = widget_button(base_btn, uvalue='use_bgimage_nothing', value=' Nothing ', xsize=0.99*x_sc, ysize=30)
    btn2 = widget_button(base_btn, uvalue='use_bgimage_map', value=' Map ', xsize=0.99*x_sc, ysize=30)
    btn3 = widget_button(base_btn, uvalue='use_bgimage_image', value=' Image ', xsize=0.99*x_sc, ysize=30)
    tmp = widget_label(left_panel5, value='', /align_left)
    tmp = widget_label(left_panel5, value='', /align_left)
    ;tmp = widget_label(left_panel5, value='Color table for the lines', /align_left)
    ;ct_lines_box = widget_text(left_panel5, uvalue='ct_lines', value='0', scr_xsize=0.1*x_sc, /editable)
    ;tmp = widget_label(left_panel5, value='', /align_left)
    ;tmp = widget_label(left_panel5, value='', /align_left)
    tmp = widget_label(left_panel5, value='Color table for the background image', /align_left)
    ct_image_box = widget_text(left_panel5, uvalue='ct_image', value='39', scr_xsize=0.1*x_sc, /editable)
    
    right_panel5 = widget_base(main_panel5, /column, xsize=0.7*x_sc, ysize=0.99*y_sc)
    tmp = widget_label(right_panel5, value='Note: subpixels are 9.4/3 arcsec size.', /align_left)
    line_map = widget_draw (right_panel5, uvalue='line_map', xsize=0.65*x_sc, ysize=0.85*y_sc)
    tmp = widget_label(right_panel5, value='Image is saved at:', /align_left)
    line_map_text = widget_text(right_panel5, value=' ', /align_left, xsize=100)
    
    
    
    
    
    ;------------------------------
    
    
    wT4 = WIDGET_BASE(wTab, TITLE=' About ', /COLUMN)
    main_panel4 = widget_base (wT4, /column, xsize=0.99*x_sc, ysize=0.91*y_sc)
    
    file = filepath(root_dir = ProgramRootDir(), subdir = 'misc/', 'pacsman_logo.bmp')
    if file_test(file) and not keyword_set(makeall) then begin
      test=QUERY_BMP(file, info)
      logo_draw = widget_draw(main_panel4, xsize = info.dimensions(0), ysize = info.dimensions(1))
    endif
    
    ;path = find_with_def('PACSman_map.pro', !path)
    ;tmp = strsplit(path, '/', /extract)
    ;path = '/'
    ;for i = 0, n_elements(tmp) - 2 do path = path + tmp[i] + '/'
    ;file = path+'misc/pacsman_logo.bmp'
    ;if file_test(file) and not keyword_set(makeall) then begin
    ;  test=QUERY_BMP(file, info)
    ;  logo_draw = widget_draw(main_panel4, xsize = info.dimensions(0), ysize = info.dimensions(1))
    ;endif
    
    ;ind = where(strpos((routine_info(/source))(*).path, 'pacsman_map') gt -1)
    ;pathto = (routine_info(/source))(ind(0)).path
    ;readcol, strtrans(pathto, 'pacsman_map.pro', 'pacsman_version.txt'), version, /silent
    
    versionfile = filepath(root_dir = ProgramRootDir(), 'pacsman_version.txt')
    readcol, versionfile, version, /silent
    
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
    lbl = widget_text(main_panel4, value = str, /wrap, /scroll, xsize = 90, ysize = 20)
    
    ;------------------------------
    
    
    widget_control, main, /realize
    
    
    widget_control, image_footprints, get_value=image_footprints_index ;for wset
    widget_control, image_subpixels, get_value=image_subpixels_index ;for wset
    widget_control, image_subpixels2, get_value=image_subpixels_index2 ;for wset
    widget_control, image_map, get_value=image_map_index ;for wset
    widget_control, image_errmap, get_value=image_errmap_index ;for wset
    widget_control, image_detmap, get_value=image_detmap_index ;for wset
    widget_control, velo_map, get_value=velo_map_index ;for wset
    widget_control, fwhm_map, get_value=fwhm_map_index ;for wset
    widget_control, cont_map, get_value=cont_map_index ;for wset
    widget_control, linecont_map, get_value=linecont_map_index ;for wset
    if file_test(file) and not keyword_set(makeall) then widget_control, logo_draw, get_value=logo_index ;for wset
    widget_control, line_map, get_value=line_map_index ;for wset
    
    if keyword_set(load_image) then widget_control, btn3, set_button = 1 else widget_control, btn1, set_button = 1
    
    if file_test(file) and not keyword_set(makeall) then begin
      device, get_decomposed=current_decomposed
      wset, logo_index
      ima = read_bmp(file)
      ima = reverse(temporary(ima), 1)
      tv, ima, true=1
      device, decomposed=current_decomposed
    endif
    
    drawsize = widget_info(image_map, /geometry)
    drawsize = [drawsize.xsize, drawsize.ysize]
    
    s.wset_image_footprints = image_footprints_index
    s.wset_image_subpixels = image_subpixels_index
    s.wset_image_subpixels2 = image_subpixels_index2
    s.wset_image_map = image_map_index
    s.wset_image_detmap = image_detmap_index
    s.wset_image_errmap = image_errmap_index
    s.wset_velo_map = velo_map_index
    s.wset_fwhm_map = fwhm_map_index
    s.wset_cont_map = cont_map_index
    s.wset_linecont_map = linecont_map_index
    s.wset_line_map = line_map_index
    s.image_map = image_map
    s.image_detmap = image_detmap
    s.velo_map = velo_map
    s.fwhm_map = fwhm_map
    s.cont_map = cont_map
    s.linecont_map = linecont_map
    s.wtab = wtab
    s.line_map_text = line_map_text
    s.files = files
    s.groupid = main
    ;s.ct_lines_box = ct_lines_box
    s.ct_image_box = ct_image_box
    s.drawsize = drawsize
    ;s.btn_update = btn_update
    s.bgimage_list = bgimage_list
    s.bgimage_text = bgimage_text
    s.calib_list = calib_list
    s.object_name = object
    s.object_text = object_text
    s.comp_text = comp_text
    s.band_text = band_text
    s.line_text = line_text
    
    widget_control, main, set_uvalue = s
  ;##################################################################
    
  endif
  ;widget_control, s.use_correctionfactor_box, set_button = 1
  
  if not keyword_set(load_image) then load_image = '' else begin
    case load_image of
      'DSS/2b': n = 2
      'DSS/2r': n = 3
      'DSS/1': n = 4
      'DSS/2i': n = 5
      '2MASS/J': n = 6
      '2MASS/H': n = 7
      '2MASS/K': n = 8
      'IRAC/ch1': n = 9
      'IRAC/ch2': n = 10
      'IRAC/ch3': n = 11
      'IRAC/ch4': n = 12
      'MIPS/24um': n = 13
      'MIPS/70um': n = 14
      'MIPS/160um': n = 15
      else: n = -1
    endcase
    if n eq - 1 then load_image = '' else begin
      s.bgimage = load_image
      if not keyword_set(makeall) then begin
        widget_control, s.bgimage_text, set_value = 'Background image selected: ' + s.bgimage
        widget_control, s.bgimage_list, set_combobox_select = n
      endif
    endelse
  endelse
  
  if not keyword_set(makeall) then begin
    case use_calib_error of
      '0': n = 0
      '10': n = 1
      '20': n = 2
      else:
    endcase
    widget_control, s.calib_list, set_combobox_select = n
  endif
  
  if keyword_set(makeall) then begin
    update_plots, s, /all, ps = ps
  endif else begin
    widget_control, main, set_uvalue = s
    xmanager, 'PACSman_map', main, /no_block ; wait for events
  endelse
  
;##################################################################
  
  
  
;##################################################################
  
  
end

