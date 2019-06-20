;+
; NAME:
;   PACSman_fit
;
; AUTHOR:
;   Vianney Lebouteiller, Diane Cormier, CEA/SAp, Saclay, France
;   vianney.lebouteiller@cea.fr
;
;   UPDATED VERSIONs can be found at:
;      http://www.myravian.fr/Homepage/Softwares.html
;
; PURPOSE:
;   Performs line fitting for Herschel/PACS spectroscopic maps
;
; MAJOR TOPICS:
;   Line Fitting
;
; CALLING SEQUENCE:
;   PACSman_fit, line, n_rasters, cube_name, cube_nod = cube_nod,
;   redshift = redshift, display = display,
;   skipfit = skipfit, use_errors = use_errors,
;   poly_degree = poly_degree, lineshift = lineshift, integrate = integrate,
;   frac = frac, fringes = fringes, strong_constraints = strong_constraints,
;   mic = mic, addlocal = addlocal, band = band, order = order
;
; DESCRIPTION:
;
;   Performs a line fitting on the cloud of data points output by the
;   HIPE PACS pipeline for PACS spectroscopy. Some input parameters
;   need to given to the program (which line to fit, how many
;   rasters, etc...).
;   Plots of fits are saved in a directory for later inspection.
;   The fit results are also saved in a text file.
;   The program saves the cube (5x5xNRASTERS) of line fluxes, errors,
;   RA, DEC, etc...
;
;   The input has to be a PACS cube. If the data is in frames, then
;   specFrames2PacsCube can be used in HIPE to make the cubes which
;   can be then exported to work with PACSman_fit.
;
; MANDATORY INPUTS:
;   line (string) - line among the available lines in the program. More lines
;          can be added. Available lines by default are 'CII157',
;          'OIII88', 'OI63', 'NII122', 'OI145', 'NIII57', and
;          'NII205'. This is also the name of the directory where the
;          cubes (that were exported from HIPE in a FITS
;         format) are expected to be.
;
;   n_cubes (integer) - number of cubes
;
;   cube_name (string) - names of cubes (before the number of the
;                        raster). Example: cube_name is 'cubeA' if
;                        files are named cubeA1.fits, cubeA2.fits, etc...
;
; RETURNS:
;
;   Nothing. Output the plots in a directory where the cubes are and create cubes of line
;   fluxes, errors, RA, DEC, etc... for later use, for instance with
;   PACSman_map to create a map projected on the sky.
;
;
; CALLING SEQUENCE:
;
; see main routine below or simply type pacsman_fit with no input parameters
;
; REFERENCES:
;
; DEPENDENCIES:
;
; MPFIT : http://cow.physics.wisc.edu/~craigm/idl/idl.html
;
; IDL ASTROLIB (frebin, readfits, writefits, sxpar) :
; http://idlastro.gsfc.nasa.gov/contents.html
;
; binary.pro (included) written by David Fanning : http://www.dfanning.com/
; The procedure is included with the PACSman package under the name
; PACSman_binary.pro, with a few modifications for handling of
; Herschel/PACS data
;
; MODIFICATION HISTORY:
;   07-23-10: added cubes of rebinned spectra
;   Written, Jun 2010, VL
;   v3.3: added combine keyword for combining spectra of the 3x3
;   v3.4: better error propagation, correcting fluxes for 1/(1+z)
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
;---------------------------------------------------------------------------------
function getextension, file, reqextension
  n = 0
  !ERROR_STATE.CODE=0
  extension = ''
  while !ERROR_STATE.CODE eq 0 and extension ne reqextension do begin
    test = readfits(file, h, exten_no=n, /silent)
    extension = strlowcase(strtrim(sxpar(h, 'EXTNAME')))
    if extension eq 'signal' then begin
      print, 'Could it be that the inputs are frames, not cubes?'
      close, /all
      retall
    endif
    n += 1
  endwhile
  ;help, file, n, extension, reqextension
  if !ERROR_STATE.CODE ne 0 then return, -1 else return, n-1
end
;---------------------------------------------------------------------------------
FUNCTION str_round, input
  x = input
  ;stupid function to round the number in the string...
  if strpos(x, 'NaN') gt -1 or strpos(x, 'Inf') gt - 1 or finite(x) eq 0 then return, 'NaN' else begin
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

;---------------------------------------------------------------------------------
;+
; NAME:
;  PACman_binary
;
; PURPOSE:
;
;   This function is used to display a binary representation of byte,
;   integer, and long integer values.
;
; AUTHOR:
;
;   FANNING SOFTWARE CONSULTING
;   David Fanning, Ph.D.
;   1645 Sheely Drive
;   Fort Collins, CO 80526 USA
;   Phone: 970-221-0438
;   E-mail: davidf@dfanning.com
;   Coyote's Guide to IDL Programming: http://www.dfanning.com/
;
; CATEGORY:
;
;   Utilities
;
; CALLING SEQUENCE:
;
;   output = Binary(theNumber)
;
; RETURN VALUE:
;
;   output:        A string array of 0s and 1s to be printed (normally), in a
;                  binary representation of the number. The number is represented with
;                  the highest bits on the left and the lowest bits on the right,
;                  when printed with the PRINT command.
;
; ARGUMENTS:
;
;  theNumber:      The number for which the user wants a binary representation.
;                  It must be BYTE, INT, or LONG.
;
; KEYWORDRS:
;
;  COLOR:          If this keyword is set, the binary representation always
;                  contains 24 bits of output.
;
;  SEPARATE:       If this keyword is set, the output is separated with space
;                  between each group of eight bits.
;
; EXAMPLE:
;
;  IDL> Print, Binary(24B)
;          0 0 0 1 1 0 0 0
;  IDL> Print, Binary(24L)
;          0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 0 0 0
;  IDL> Print, Binary(24L, /COLOR)
;          0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 0 0 0
;  IDL> Print, Binary(24L, /COLOR, /SEPARATE)
;          0 0 0 0 0 0 0 0    0 0 0 0 0 0 0 0    0 0 0 1 1 0 0 0
;
; MODIFICATION HISTORY:
;
;  Written by: David W. Fanning, November 10, 2007.
;  Fixed a problem with error handling. 13 March 2008. DWF.
;  Added noreverse for handling of Herschel/PACS data. Jun 2010. VL
;-
;******************************************************************************************;
;  Copyright (c) 2008, by Fanning Software Consulting, Inc.                                ;
;  All rights reserved.                                                                    ;
;                                                                                          ;
;  Redistribution and use in source and binary forms, with or without                      ;
;  modification, are permitted provided that the following conditions are met:             ;
;                                                                                          ;
;      * Redistributions of source code must retain the above copyright                    ;
;        notice, this list of conditions and the following disclaimer.                     ;
;      * Redistributions in binary form must reproduce the above copyright                 ;
;        notice, this list of conditions and the following disclaimer in the               ;
;        documentation and/or other materials provided with the distribution.              ;
;      * Neither the name of Fanning Software Consulting, Inc. nor the names of its        ;
;        contributors may be used to endorse or promote products derived from this         ;
;        software without specific prior written permission.                               ;
;                                                                                          ;
;  THIS SOFTWARE IS PROVIDED BY FANNING SOFTWARE CONSULTING, INC. ''AS IS'' AND ANY        ;
;  EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES    ;
;  OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT     ;
;  SHALL FANNING SOFTWARE CONSULTING, INC. BE LIABLE FOR ANY DIRECT, INDIRECT,             ;
;  INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED    ;
;  TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;         ;
;  LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND             ;
;  ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT              ;
;  (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS           ;
;  SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.                            ;
;******************************************************************************************;
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
;______________________________________________________________________________
;converts binary mask into a boolean mask
function convert2bool, maskin

  nn = n_elements(maskin(0, 0, *)) * 32 ;31 or 32 ? 32 for chop/nod at least
  maskout = intarr(5, 5, nn)
  for x = 0, 4 do for y = 0, 4 do begin
    arr = intarr(1)
    for i = long(0), n_elements(maskin(x, y, *)) - 1 do arr = [arr, nint(PACSman_binary(maskin(x, y, i), /noreverse))] ;somehow the mask extensions are the transposed matrix (5x5)
    maskout(x, y, *) = arr[1:nn]
  endfor
  return, maskout
end
;______________________________________________________________________________
;print on screen and on file
pro write_message, textarr, lun
  for nn = 0, n_elements(textarr) - 1 do begin
    printf, lun, a2str(textarr[nn])
    print, textarr[nn]
  endfor
end
;______________________________________________________________________________
;baseline + asymetric gaussian model
;function SpectrumModel, x, p, poly_degree = poly_degree;
;
;  result = 0.
;  for i = 0, poly_degree do result +=  p[i+4]*x^float(i)
;
;  return, result + $
;    (p[0] * exp(-(x-p[1])^2./(2.*p[2]^2.)))
;end
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
pro spec_smooth, lambda, flux, lambda_smooth, flux_smooth, error_smooth, $
    lambda_ref = lambda_ref, bin = bin, error_indiv = error_indiv, line_range = line_range, getbin = getbin, goodinds = goodinds, lun = lun
  ;error_indiv is the error on individual points in the datacloud (~cloud dispersion)
  ;error_smooth is the error in each wavelength bin
    
  if not keyword_set(lambda_ref) then lambda_ref = lambda
  bin = nint(bin*10000.)*0.0001 ;round bin value to 3 digits
  
  if keyword_set(getbin) then begin
  
    if keyword_set(lun) then write_message, 'Initial bin size [um]: ' + a2str(bin), lun else print, 'Initial bin size [um]: ', bin
    loop = 0
    repeat begin
    
      test = 1
      
      lambda_smooth = dblarr( 1+floor((max(lambda_ref, /nan)-min(lambda_ref, /nan)) / bin) )
      
      ;first check bin size, in order to have much sampling at least around line
      ninbins = intarr(n_elements(lambda_smooth))
      ;print, line_range
      k = long(0)
      for l = min(lambda_ref, /nan), max(lambda_ref, /nan), bin do begin
        if abs(l) gt line_range then continue
        ind = where( (lambda ge l) and (lambda lt l+bin) and finite(flux), c )
        ninbins[k] = c
        k += 1
      endfor
      if k eq 0 then begin
        ;print, 'Some wavelength bins have small number of measurements, consider decrease the sampling. '
        ;print, 'Suggested sampling: '
        ;print, ninbins
        test = 0
        bin *= 1.5
      endif else begin
        ninbins = ninbins[0:k-1]
        if min(ninbins, /nan) lt 20. then begin
          ;print, 'Some wavelength bins have small number of measurements, consider decrease the sampling. '
          ;print, 'Suggested sampling: '
          ;print, ninbins
          test = 0
          bin *= 1.5
        endif
      endelse
      
      loop += 1
      if loop gt 10. then begin
        print, 'There was a problem...'
        retall
      endif
      
    endrep until test eq 1
    ;more conservative bin
    bin *= 1.5
    if keyword_set(lun) then write_message, 'Final bin size [um]: ' + a2str(bin), lun else print, 'Final bin size [um]: ', bin
  endif
  lambda_smooth = dblarr( 1+floor((max(lambda_ref, /nan)-min(lambda_ref, /nan)) / bin) )
  
  ;then do flux_smooth
  flux_smooth = lambda_smooth
  ;error_smooth = 0. * flux_smooth
  k = long(0)
  for l = min(lambda_ref, /nan), max(lambda_ref, /nan), bin do begin
    lambda_smooth(k) = l + 0.5*bin
    ;ind = where( (lambda ge l) and (lambda lt l+bin) and finite(flux), c )
    ind = where( (lambda ge l-0.5*bin) and (lambda lt l+1.5*bin) and finite(flux), c )
    if c lt 5 then flux_smooth(k) = !values.f_nan else flux_smooth(k) = median( flux(ind) )
    k += 1
  endfor
  
  error_indiv = 0. * flux_smooth
  ninbins = intarr(n_elements(flux_smooth))
  ;finally do the errors, using flux_smooth to calculate dispersion around line shape
  k = long(0)
  goodinds = intarr(1)
  for l = min(lambda_ref, /nan), max(lambda_ref, /nan), bin do begin
  
    if finite(flux_smooth(k), /nan) then begin
      ;error_smooth(k) = !values.f_nan
      error_indiv(k) = !values.f_nan
      k += 1
      continue
    endif
    
    ;ind = where( (lambda ge l) and (lambda lt l+bin) and finite(flux), c )
    ind = where( (lambda ge l-0.5*bin) and (lambda lt l+1.5*bin) and finite(flux), c )
    if c lt 5 then begin
      ;error_smooth(k) = !values.f_nan
      error_indiv(k) = !values.f_nan
    endif else begin
      ninbins(k) = c ;n_elements(ind)
      
      ;error_smooth(k) = 0.5*median(abs(flux(ind) - flux_smooth(k)))/0.6744897501960817D ;median absolute deviation. 0.67 is inverseErf(0.5) * sqrt(2), to approx the sigma
      ;error_smooth(k) /= sqrt(ninbins(k))
      
      ;http://www.stat.berkeley.edu/~bradluen/stat2/lecture28.pdf
      tmp = flux(ind)-interpol(flux_smooth, lambda_smooth, lambda(ind))
      val = robust_sigma(tmp, /zero, goodvec = goodvec); we use goodvec to keep the indices of the non-clipped pixels for the chisq determination
      if val eq -1 then begin
        print, 'replacing by stdev value'
        val = stdev(tmp)
        goodinds = [goodinds, ind]
      endif else goodinds = [goodinds, ind[goodvec]]
      error_indiv(k) = val
    ;val /= sqrt(ninbins(k))
    ;error_smooth(k) = val
      
    endelse
    k += 1
  endfor
  goodinds = goodinds[1:n_elements(goodinds)-1]
  ;we might have some duplicates if we choose +/-1.5 bin
  goodinds = goodinds[uniq(goodinds, sort(goodinds))]
  
  ;take care of NaNs
  ;we need to replace them because errors are used for weights later (error_indiv), and most importantly for
  ;sigma-clipping (error_smooth)
  ind = where( not finite(error_indiv), c, complement = noind)
  if c gt 0 then begin
    if noind(0) ne -1 then error_indiv(ind) = interpol(error_indiv(noind), lambda_smooth(noind), lambda_smooth(ind))
  endif
  error_indiv = abs( smooth(error_indiv, 3, /edge_truncate, /nan) )
  error_smooth = error_indiv / sqrt(ninbins)
  
end
;______________________________________________________________________________
;______________________________________________________________________________
;______________________________________________________________________________
;MAIN PROGRAM
pro PACSman_fit, line, n_cubes, cube_name, cube_nod = cube_nod, redshift = redshift, display = display, $
    use_errors = use_errors, poly_degree = poly_degree, lineshift = lineshift, $
    fringes = fringes, strong_constraints = strong_constraints, mic = mic, addcomponents = addcomponents, components_fwhm = components_fwhm, band = band, order = order, $
    no_correction_factor = no_correction_factor, use_correction_factor = use_correction_factor, broadening = broadening, constraints_wavelength = constraints_wavelength, $
    savedatacloud = savedatacloud, constraints_continuum = constraints_continuum, sampling = sampling, viewrange = viewrange, $
    find = find, central3by3combine = central3by3combine, central5by5combine = central5by5combine, pls_corr = pls_corr, align = align, rebinned = rebinned, obsid = obsid, help = help, cleanmask = cleanmask, $
    use_velomap = use_velomap, montecarlo = montecarlo, blend = blend, combine_rasters = combine_rasters, subtract_rasters = subtract_rasters, dontdo_plots_publi = dontdo_plots_publi
  ;asymetric = asymetric
    
    
    
    
  ;==========================================================================================
  ;==========================================================================================
  ;==========================================================================================
  ;CHANGE HERE THE WAVELENGTHS AND LABELS OF THE LINES
  line_label      = ['CII157', 'OI63',     'OIII88', 'OI145',     'NII122',  'NIII57',  'NII205',     'OIII52',    'OH79', 'OH119', 'OIV26', 'OH65', 'OH84', 'CO144', 'OH1192',  'A2',   'A3',   'A4']
  line_wavelength = [ 157.7409, 63.183705, 88.356,   145.525439, 121.89806,    57.33,    205.1782,    51.8145,    79.,       119., 25.8903,    65.,   84.,    144.784, 119.220, 119.441, 118.663, 119.848] ;rest-frame lambda in um
  ;from NIST line_wavelength = [ 157.6790, 63.1852, 88.3564,   145.535, 121.800,    57.317,    205.300,    51.8145, 79., 0., 0., 0., 0. ] ;rest-frame lambda in um
  line_label      = [line_label, 'CII157indiv', 'OI63indiv', 'OIII88indiv']
  line_wavelength = [line_wavelength, 157.7409, 63.183705, 88.356]
  ;==========================================================================================
  ;==========================================================================================
  ;==========================================================================================
  
  
  
  ;==========================================================================================
  ;==========================================================================================
  ;==========================================================================================
  if keyword_set(help) or not keyword_set(line) or not keyword_set(n_cubes) or not keyword_set(cube_name) then begin
    print, 'Calling sequence.'
    print, ''
    print, "Minimum calling sequence example: "
    print, "  pacsman_fit, 'CII157'    {line label},"
    print, "                3          {number of raster positions for a map (number of cubes, e.g., 10 or cubeA1...A10.fits).},"
    print, "               'cubeA'     {cube name(s)}"
    print, ""
    print, "Available lines (edit pacsman_fit.pro to change the labels and wavelengths): "
    print, line_label
    print, ""
    print, "Optional parameters:"
    print, "  cube_nod               = 'cubeB'             {if chop/nod observation, name of other nod cubes}"
    print, "  redshift               = 0.        {z}       {Redshift. Line flux is modified accordingly}"
    print, "  lineshift              = 0.        {km s-1}  {artificial velocity shift assumed to be due to wavelength calibration uncertainties. Line flux is not modified}"
    print, "  viewrange              = [10, 10]  {FWHM}    {how many FHWMs on either side of the line will be used only for the plots}"
    print, "  poly_degree            = 2...5               {polynomial degree for continuum fit, must be <=5. Default is 2}"
    print, "  broadening             = [0, 500]  {km s-1}  {allowed intrinsic min/max line broadening. If only one value given, sets the maximum broadening}"
    print, "  constraints_continuum  = [2, 6, 6] {FWHM}    {first parameter is +/- range around the line, the 2 other parameters define the continuum range}"
    print, "  constraints_wavelength = 0.5       {FWHM}    {allowed shift around the expected central wavelength}"
    print, "  addcomponents          = [0, ...]  {km s-1}  {add other component. Values are the radial velocities. }"
    print, "  components_fwhm        = [0, ...]  {km s-1}  {*Intrinsic* FWHM of the other components}"
    print, "  sampling               = 10.                 {bin size is calculated as fwhm/sampling for the smoothed spectra (used for errors and plots). Bin size is usually calculated automatically unless sampling is specified}"
    print, "  subtract_rasters      =  [1, ...]            {subtract the spectra from input array of rasters, when a position in the map is free of signal. Uses the saved file from a previous run}"
    print, "  band                   = 'red'               {force band value, for early data}"
    print, "  pls_corr               = 0..3                {Don't use unless you know exactly why you use it. Use pacsman_footprint to calculate total fluxes. NOTE: choice is automatically made if central3by3combine is set. }"
    print, "                                                   0: no point-source correction, 1: correct for fraction in central spaxel, 2: correct for fraction in 3x3, 3: correct for fraction in 5x5}"
    print, "  order                  = '1'                 {force order value, for early data}"
    print, "  /find     or      find = 2.        {FWHM}    {DEACTIVATED. finds the source around +/- value by a first guess call to gaussfit rather than assuming a centroid from the redshift}"
    print, ""
    print, "Optional switches:"
    print, "  /combine_rasters       {combine the rasters that are spatially coinciding, i.e., combine spectra from cycles before fitting}"
    print, "  /nonegativefluxes      {negative line fluxes will be overridden to 0}"
    print, "  /obsid                 {obsid wil be appended to output directories}"
    print, "  /cleanmask             {force recalculation of mask files}"
    print, "  /use_correction_factor {depends on the calibration version number, i.e., time dependence <13 in getCalTree(obs=obs).spectrometer}"
    print, "  /strong_constraints    {tighter constraints on FWHM and radial velocity parameters, to use for more reliable upper limits}"
    print, "  /display               {display fits, won't save them}"
    print, "  /fringes               {add a sinusoidal component to model the continuum}"
    print, "  /savedatacloud         {saves the full data cloud. The smoothed spectrum is always saved by default}"
    print, "  /rebinned              {works on the rebinned cubes instead of the data clouds}"
    print, "  /central3by3combine    {combine the central 3x3 spaxel spectra before fitting, see also /align switch}"
    print, "  /central5by5combine    {combine all (5x5) the spaxel spectra before fitting, see also /align switch}"
    print, "  /align                 {align spaxel spectra of the 3x3 central spaxels before combining, only used in conjunction with /central3by3combine switch}"
    print, "  /dontdo_plots_publi    {By default, prepare line fit plots for publications, use this option to skip}"
    print, "  /montecarlo            {estimate errors using Monte-Carlo algorithm [SLOW!], can also give montecarlo=n_iterations (defaut=100)}"
    print, "  /blend                 {option for blended lines (when 2nd and/or 3rd components are set). Will force the line fluxes to be positive}"
    print, "  /use_velomap           {ONLY FOR TESTS: use the velocity map from a previous run as a way to better estimate the 1st guess on the line velocity}"
    print, ''
    help, line, cube_name
    if keyword_set(line) then begin
      print, ''
      print, 'Previous call:'
      print, ''
      if file_test(line+'/calling_parameters.txt') then spawn, '\cat '+line+'/calling_parameters.txt'
    endif
    retall
  endif
  ;==========================================================================================
  ;==========================================================================================
  ;==========================================================================================
  
  
  
  
  ;==========================================================================================
  ;==========================================================================================
  ;==========================================================================================
  ;some checks
  if keyword_set(components_fwhm) and keyword_set(addcomponents) then if n_elements(components_fwhm) ne n_elements(addcomponents) then begin
    print, 'components_fwhm and addcomponents arrays dont have the same number of elements. '
    retall
  endif
  if keyword_set(addcomponents) and not keyword_set(components_fwhm) then components_fwhm = addcomponents * 0.
  
  if (n_cubes eq 1) and keyword_set(combine_rasters) then begin
    print, '/!\ The /combine_rasters keyword was set, but only one cube was given'
    retall
  endif
  
  ;check is line exists
  ind = where( line_label eq strtrim(line), c )
  line_index = ind(0)
  if c eq 0 then begin
    print, 'Line not found. If more lines are needed, modify the program located at: '
    which, 'pacsman_fit'
    print, ''
    print, 'Input:           ', line
    print, 'Available lines: ', line_label
    retall
  endif
  
  if not keyword_set(mic) then mic = '' else mic = 'mic'
  if not keyword_set(obsid) then obsid = '' else obsid = '_obsid' + a2str(obsid)
  line = line_label(line_index) + mic + obsid
  
  dir = line
  output_name = cube_name       ;output for writing results
  if keyword_set(cube_nod) then output_name += cube_nod
  if keyword_set(central3by3combine) then output_name += '_3x3'
  if keyword_set(central5by5combine) then output_name += '_5x5'
  
  if file_test(line+'/calling_parameters.txt') then spawn, '\mv ' + line+'/calling_parameters.txt ' + line+'/calling_parameters_previous.txt'
  
  ;save results from previous run if some rasters have to be subtracted in a 2nd run
  if keyword_set(subtract_rasters) then begin
    if file_test(dir+'/'+output_name+'.sav') eq 0 then begin
      print, '/!\ The previous result saved file was not found: ' +  dir+'/'+output_name+'.sav'
      retall
    endif
    restore, filename = dir+'/'+output_name+'.sav' ;restore cube_spectra saved from previous iteration
    cube_spectra_before_subtracting = cube_spectra
    spawn, '\mv ' + dir+'/'+output_name+'.sav ' + dir+'/'+output_name+'_before_subtracting.sav'
    spawn, '\mv ' + dir+'/'+line_label(line_index)+'_fitresults_'+output_name+'.txt ' + dir+'/'+line_label(line_index)+'_fitresults_'+output_name+'_before_subtracting.txt'
    spawn, '\mv ' + dir + '/PLOTS_' + output_name + ' ' + dir + '/PLOTS_' + output_name + '_before_subtracting'
    if file_test(dir + '/PLOTSpubli_' + output_name, /dir) then spawn, '\mv ' + dir + '/PLOTSpubli_' + output_name + ' ' + dir + '/PLOTSpubli_' + output_name + '_before_subtracting'
    spawn, '\mv ' + line+'/calling_parameters.txt ' + line+'/calling_parameters_before_subtracting.txt'
  endif
  
  ;save calling parameters
  str = "pacsman_fit"
  if keyword_set(line) then str += ", '" + line_label(line_index) + "'"
  if keyword_set(n_cubes) then str += ", " + a2str(n_cubes)
  if keyword_set(cube_name) then str += ", '" + cube_name + "'"
  if keyword_set(cube_nod) then str += ", cube_nod='" + cube_nod + "'"
  if keyword_set(redshift) then str += ", redshift=" + a2str(trim(redshift))
  if keyword_set(broadening) then if n_elements(broadening) eq 1 then str += ", broadening="+a2str(broadening) else str += ", broadening=["+a2str(trim(broadening(0)))+","+a2str(trim(broadening(1)))+"]"
  if keyword_set(constraints_wavelength) then str += ", constraints_wavelength=" + a2str(constraints_wavelength)
  if keyword_set(viewrange) then begin
    if n_elements(viewrange) eq 1 then str += ", viewrange="+a2str(trim(viewrange(0))) $
    else str += ", viewrange=["+a2str(trim(viewrange(0)))+","+a2str(trim(viewrange(1)))+"]"
  endif
  if keyword_set(fringes) then str += ", /fringes"
  if keyword_set(strong_constraints) then str += ", /strong_constraints"
  if keyword_set(addcomponents) then begin
    str += ", addcomponents=["
    for ncomp = 0, n_elements(addcomponents) - 1 do begin
      str += a2str(trim(addcomponents(ncomp)))
      if ncomp lt n_elements(addcomponents) - 1 then str += ","
    endfor
    str += "]"
  endif
  if keyword_set(components_fwhm) then begin
    str += ", components_fwhm=["
    for ncomp = 0, n_elements(components_fwhm) - 1 do begin
      str += a2str(trim(components_fwhm(ncomp)))
      if ncomp lt n_elements(components_fwhm) - 1 then str += ","
    endfor
    str += "]"
  endif
  if keyword_set(use_correction_factor) then str += ", /use_correction_factor"
  if keyword_set(sampling) then str += ", sampling=" + a2str(trim(sampling))
  if keyword_set(find) then str += ", /find"
  if keyword_set(central3by3combine) then str += ", /central3by3combine"
  if keyword_set(central5by5combine) then str += ", /central5by5combine"
  if keyword_set(combine_rasters) then str += ", /combine_rasters"
  if keyword_set(mic) then str += ", /mic"
  if keyword_set(use_velomap) then str += ", /use_velomap"
  if keyword_set(obsid) then str += ", obsid = '"+a2str(strtrans(obsid, '_obsid', ''))+"'"
  if keyword_set(monte) then str += ", /montecarlo"
  ;if keyword_set(plots_publi) then str += ", /plots_publi"
  
  if keyword_set(band) then str += ", band=" + band
  if keyword_set(order) then str += ", order=" + a2str(order)
  
  print, 'string to save in pacsman_call.txt: ', str
  openw, lun2, line+'/calling_parameters.txt', /get_lun
  printf, lun2, str
  openw, lun, dir+'/'+line_label(line_index)+'_fitresults_'+output_name+'.txt', width = 400, /get_lun
  write_message, str, lun
  
  if keyword_set(constraints_continuum) then begin
    printf, lun2, 'constraints_continuum: ', constraints_continuum
    write_message, 'constraints_continuum:', lun
    for i = 0, n_elements(constraints_continuum) - 1 do write_message, a2str(constraints_continuum(i)), lun
  endif
  close, /all
  
  openw, lun, dir+'/'+line_label(line_index)+'_fitresults_'+output_name+'.txt', width = 400, /get_lun
  ;==========================================================================================
  ;==========================================================================================
  ;==========================================================================================
  
  
  
  
  ;==========================================================================================
  ;==========================================================================================
  ;==========================================================================================
  ;INPUT PARAMETERS
  
  ;************************************
  ;PARAMETERS THAT CAN BE CHANGED
  if not keyword_set(constraints_continuum) then constraints_continuum = [2., 6., 6.]
  if n_elements(constraints_continuum) eq 2 then constraints_continuum = [constraints_continuum, constraints_continuum(1)]
  
  if not keyword_set(sampling) then begin
    sampling = 30.
    getbin = 1 ;for spec_smooth, binsize will be automatically calculated
  endif
  sampling = float(sampling)
  
  if keyword_set(strong_constraints) then begin
    if not keyword_set(broadening) then broadening = [0., 10.] ; allowed kms-1 to add to the fwhm at the observed wavelength
    if not keyword_set(constraints_wavelengths) then constraints_wavelength = 0.1   ; allowed shift in central_wavelength in units of fwhm
  endif else begin
    if not keyword_set(broadening) then broadening = [0., 50.] ; allowed kms-1 to add to the fwhm at the observed wavelength
    if n_elements(broadening) eq 1 then broadening = [0., broadening]
    if not keyword_set(constraints_wavelength) then constraints_wavelength = 0.5   ; allowed shift in central_wavelength in units of fwhm
  endelse
  
  ;if keyword_set(addsecond) then if n_elements(addsecond) eq 1 then addsecond = [addsecond, 0.]
  ;if keyword_set(addthird) then if n_elements(addthird) eq 1 then addthird = [addthird, 0.]
  
  if constraints_wavelength gt constraints_continuum(0) then begin
    str = "/!\ Warning: the range allowed for the line center wavelength is overlaps with the range for the continuum 1st guess."
    write_message, str, lun
  endif
  
  if keyword_set(redshift) then begin
    redshift = redshift(0)
    str = ["/!\ Redshift will be used and line fluxes will be modified accordingly.", $
      "    Use the 'lineshift' keyword to shift the line with no change for the flux."]
    write_message, str, lun
  endif else redshift = 0.
  if not keyword_set(lineshift) then lineshift = 0.
  if not keyword_set(ignore_errors) then ignore_errors = 1 ;force non-treatment of error plane until errors are well estimated by the pipeline
  
  
  ;************************************
  ;PARAMETERS THAT CANNOT BE CHANGED
  n_iter_noise = 100.
  if keyword_set(montecarlo) then if montecarlo gt 1 then n_iter_noise = montecarlo else n_iter_noise = 100 ;number of iterations to calculate error in integrated flux
  n_sigma_clipping = 15 ;200; datacloud is *very* dispersed, we don't wanna throw away everything just yet
  
  
  ;************************************
  ;INITIALIZE SOME OTHER STUFF
  
  ;cleanplot, /sil
  set_plot, 'X'
  !p.font = 0
  !p.thick = 2
  ;loadct, 0, /sil
  ;cleanplot, /sil
  !p.position = 0;[0.15, 0.15, 0.95, 0.9]
  !p.multi = 0
  ;loadct, 0, /sil
  ;DefSysV, '!RNG', Obj_New('RandomNumberGenerator')
  cspeed = 299792.458 ;kms
  
  ;resolution of PNG plots
  reso = [640, 480]
  
  ;if keyword_set(dsm) then dir = './PACS_spec/fits/' + line_label(line_index) + 'mic' else dir = './' + line_label(line_index) + mic
  
  spawn, '\rm ' + './' + line_label(line_index) + '/*_{fitparams,fiterrors,spectra,res,datacloud}.fits'
  spawn, '\rm ' + './' + line_label(line_index) + '/*_{fitparams,fiterrors,spectra,res,datacloud}.sav'
  
  ;test header
  test = readfits(dir+'/'+cube_name+'1.fits', h, /silent)
  ;position angle
  position_angle = !values.f_nan
  tmp = sxpar(h, 'POSANGLE', count = c) ;position angle in degrees
  if c gt 0 then position_angle = float(tmp)
  
  
  ;************************************
  ;find band and order
  if not keyword_set(band) then begin
    ;check band
    band = a2str(sxpar(h, 'BAND', count = c))
    help, band
    if c eq 0 or band eq 'N/A' then begin
      band = ''
      read, band, format = '(a)', prompt = 'Band (red/blue) : '
      ;check order
      order = 0
      if band eq 'red' then order = 1 else read, order, format = '(i)', prompt = 'Order (1 for red band 2/3 for blue band): '
    endif else begin
      if band eq 'B3A' or strmid(band, 0, 10) eq 'Blue short' then begin
        band = 'blue'
        order = 3
      endif
      if band eq 'B2B' or strmid(band, 0, 9) eq 'Blue long' then begin
        band = 'blue'
        order = 2
      endif
      if band eq 'B2A' then begin
        band = 'blue'
        order = 2
      endif
      if band eq 'R1' or strmid(band, 0, 11) eq 'Red channel' then begin
        band = 'red'
        order = 1
      endif
      if strpos(band, 'Red channel') gt -1 then order = 1 ; else if strpos(band, 'order=') gt -1 then order = nint(strmid(band, strpos(band, 'order=')+6, 1)) else order = 1
      if strpos(band, 'Blue channel') gt -1 then if strpos(band, 'order=') gt -1 then order = nint(strmid(band, strpos(band, 'order=')+6, 1)) else read, order, format = '(i)', prompt = 'Order (1 for red band 2/3 for blue band): '
    endelse
  endif else begin
    if band eq 'red' then order = 1 ;there's only one order (1) in band R1
    if ((band ne 'red') and (band ne 'blue')) or ((order ne 1) and (order ne 2) and (order ne 3)) or (band eq 'red' and order ne 1) or (band eq 'blue' and order eq 1) then begin
      str = 'Wrong values for band (red,blue) and/or order (1,2,3)'
      write_message, str, lun
      help, band, order
    endif
  endelse
  
  
  ;************************************
  ;check resolution
  case order of
    1: begin
      wave = [105, 158, 175, 210] ;um
      R = [318, 239, 212, 140] ;km s-1
    end
    2: begin
      wave = [75, 90] ;um
      R = [156, 121] ;km s-1
    end
    3: begin
      wave = [55, 60, 72] ;um
      R = [114, 98, 55]
    end
    else: begin
      print, 'order not found'
      close, /all
      retall
    end
  endcase
  
  lambda_rest = line_wavelength(line_index)
  lambda_obs = lambda_rest * (1.+redshift)
  str = ['Lambda_rest, Lambda_obs [um]: ' + a2str(lambda_rest) + ' ' + a2str(lambda_obs)]
  write_message, str, lun
  
  lambda_obs += lambda_obs*lineshift/cspeed
  fwhm_kms = interpol(R, wave, lambda_obs)
  str = ['Expected FWHM [km s-1]: ' + a2str(fwhm_kms)]
  write_message, str, lun
  
  ;scale the expected continuum window
  constraints_continuum *= sqrt( broadening(1)^2. + fwhm_kms^2. )/fwhm_kms
  str = ['Continuum [xFWHMs]: ', string(constraints_continuum)]
  write_message, str, lun
  
  if not keyword_set(viewrange) then viewrange = constraints_continuum(1:2) else begin
    if n_elements(viewrange) eq 1 then viewrange = [viewrange[0], viewrange[0]]
    if max(viewrange le 0.) eq 1 then begin
      str = 'viewrange parameter cannot be negative. See help. '
      write_message, str, lun
      close, /all
      retall
    endif
  endelse
  
  broadening_sigma = broadening / (2.*sqrt(2.*alog(2.))) ;reason in sigma from now on
  
  fwhm_tot = lambda_obs * sqrt(fwhm_kms^2.+broadening(1)^2.)/cspeed ;i km s-1; for rebinning purposes
  
  file_delete, dir + '/PLOTS_' + output_name + '/', /recursive, /allow_nonexistent, /quiet
  file_mkdir, dir + '/PLOTS_' + output_name + '/'
  file_mkdir, dir + '/PLOTS_' + output_name + '/raw'
  file_mkdir, dir + '/PLOTS_' + output_name + '/continuum_normalized'
  file_mkdir, dir + '/PLOTS_' + output_name + '/continuum_subtracted'
  file_mkdir, dir+'/PLOTS_' + output_name + '/background'
  
  if not keyword_set(dontdo_plots_publi) then begin
    file_delete, dir + '/PLOTSpubli_' + output_name + '/', /recursive, /allow_nonexistent, /quiet
    file_mkdir, dir + '/PLOTSpubli_' + output_name + '/'
    file_mkdir, dir + '/PLOTSpubli_' + output_name + '/raw'
    file_mkdir, dir + '/PLOTSpubli_' + output_name + '/continuum_normalized'
    file_mkdir, dir + '/PLOTSpubli_' + output_name + '/continuum_subtracted'
  endif
  
  if not keyword_set(poly_degree) then poly_degree = 2
  if poly_degree gt 5 then begin
    str = " /!\ poly_degree must be <= 5"
    write_message, str, lun
    close, /all
    retall
  endif
  ;extra_degree = 3 ;(for fringes)
  
  
  ;************************************
  ;CREATE STRUCTURES
  cube_params = {line: line_label(line_index), reference_fwhm: float(fwhm_kms), band: band, order: nint(order), position_angle: float(position_angle), lambda_rest: float(lambda_rest), redshift: double(redshift), fringes: keyword_set(fringes), $
    poly_degree: nint(poly_degree), broadening: float(broadening), constraints_continuum: float(constraints_continuum), range: keyword_set(range), addcomponents: keyword_set(addcomponents), addthird: keyword_set(addthird), $
    lineshift: float(lineshift), strong_constraints: keyword_set(strong_constraints), n_rasters: nint(0), n_cubes: nint(n_cubes), pacsman_version: '', date:'', hipe_version:'', obsid: ''}
    
  if n_elements(addcomponents) gt 0 then tmp = fltarr(n_elements(addcomponents)) else tmp = 0.
  
  cube_lres = replicate({spaxel_X: nint(0), spaxel_Y: nint(0), raster: nint(0), flux: double(0.), error: double(0.), components_flux: double(tmp), components_error: double(tmp), velocity: float(0.), velocity_error: float(0.), fwhm: float(0.), fwhm_error: float(0.), $
    chisq: float(0.), ra: double(0.), dec: double(0.), wave: float(0.), integrated_flux: double(0.), integrated_noise: double(0.), continuum: float(0.), continuum_error: float(0.), eqw: float(0.), broadening: float(0.), broadening_error: float(0.)}, 25*n_cubes)
    
  cube_groups = replicate({filename: '', group: nint(-1)}, 25*n_cubes)
  
  ;if keyword_set(fullrange) then npoints = 5000. else npoints = 500. ;for smoothed spectrum
  npoints = 8000
  cube_spectra = replicate({spaxel_X: nint(0), spaxel_Y: nint(0), raster: nint(0), lambda: fltarr(npoints)*!values.f_nan, flux: fltarr(npoints)*!values.f_nan, $
    normalized_flux: fltarr(npoints)*!values.f_nan, cont_subtracted_flux: fltarr(npoints)*!values.f_nan, lambda_obs: fltarr(npoints)*!values.f_nan, error: fltarr(npoints)*!values.f_nan, ra: double(0.), dec: double(0.)}, 25*n_cubes)
    
  ;cube_fitparams = replicate({spaxel_X: 0, spaxel_Y: 0, raster: 0, params: fltarr(9+poly_degree+1), params_errors: fltarr(9+poly_degree+1)}, 5*5*n_rasters)
    
  indparams_fringes     = 3+3*n_elements(addcomponents) + indgen(3)
  indparams_cont        = 3+3*n_elements(addcomponents) + 3 + indgen(poly_degree+1)
  indparams_contfringes = 3+3*n_elements(addcomponents) + indgen(3+poly_degree+1)
  nparams               = 3+3*n_elements(addcomponents) + 3 + poly_degree+1
  cube_fitparams = replicate({spaxel_X: nint(0), spaxel_Y: nint(0), raster: nint(0), params: fltarr(nparams), params_errors: fltarr(nparams)}, 25*n_cubes)
  
  
  
  ;************************************
  k = 0
  for n = 1, n_cubes do for x = 1, 5 do for y = 1, 5 do begin
    cube_lres(k).raster = n
    cube_lres(k).spaxel_X = x
    cube_lres(k).spaxel_Y = y
    cube_spectra(k).spaxel_X = x
    cube_spectra(k).spaxel_Y = y
    cube_spectra(k).raster = n
    cube_fitparams(k).spaxel_X = x
    cube_fitparams(k).spaxel_Y = y
    cube_fitparams(k).raster = n
    k += 1
  endfor
  
  versionfile = filepath(root_dir = ProgramRootDir(), 'pacsman_version.txt')
  readcol, versionfile, version, /silent
  
  ;ind = where(strpos((routine_info(/source))(*).path, 'pacsman_fit') gt -1)
  ;pathto = (routine_info(/source))(ind(0)).path
  ;str = ['PACSman file: ', pathto]
  ;write_message, str, lun
  ;readcol, strtrans(pathto, 'pacsman_fit.pro', 'pacsman_version.txt'), version, /silent
  str = ['PACSman version: ', a2str(trim(version(0)))]
  write_message, str, lun
  cube_params.pacsman_version = a2str(trim(version(0)))
  cube_params.date = systime()
  
  test = readfits(dir+'/'+cube_name+a2str(1)+'.fits', h, /sil)
  cube_params.hipe_version = trim(sxpar(h, 'CREATOR')) + '/' + trim(sxpar(h, 'HCSS____'))
  cube_params.obsid = trim(sxpar(h, 'OBS_ID'))
  ;==========================================================================================
  ;==========================================================================================
  ;==========================================================================================
  
  
  
  
  
  
  
  
  ;==========================================================================================
  ;==========================================================================================
  ;==========================================================================================
  ;PLS correction?
  if not keyword_set(pls_corr) then pls_corr = 0
  if keyword_set(central3by3combine) then pls_corr = 2
  if keyword_set(central5by5combine) then pls_corr = 3
  PointSourceLossCalibFile = filepath(root_dir = ProgramRootDir(), subdir = 'calib/', 'PACS_PointSourceCorrection.dat')
  readcol, PointSourceLossCalibFile, pls_lambda, pls_frac_in1, tmp, pls_frac_in3x3, pls_frac_in5x5, /sil
  ;case pls_corr of
  ;  1: begin
  ;    pls_frac_lambda = [50, 60, 70, 80, 90, 100, 180, 220]
  ;    pls_frac = [0.75, 0.71, 0.70, 0.695, 0.69, 0.675, 0.45, 0.335]
  ;  end
  ;  2: begin
  ;    calibfile3by3 = filepath(root_dir = ProgramRootDir(), subdir = 'calib/', 'Pointsource_3x3_over_total.txt')
  ;    readcol, calibfile3by3, pls_frac3by3_lambda, pls_frac3by3, /silent
  ;    pls_frac_lambda = pls_frac3by3_lambda
  ;    pls_frac = pls_frac3by3
  ;  end
  ;  else:
  ;endcase
  ;==========================================================================================
  ;==========================================================================================
  ;==========================================================================================
  
  
  
  
  
  ;###########################################
  ;READ CUBE'S EXTENSIONS
  ;!ERROR_STATE.CODE=0
  ;extensions = strarr(25)
  ;n = 0
  ;while !ERROR_STATE.CODE eq 0 do begin
  ;  test = readfits(dir+'/'+cube_name+'1.fits', h, exten_no=n, /silent)
  ;  extensions(n) = strlowcase(strtrim(sxpar(h, 'EXTNAME')))
  ;  ;print, extensions(n), size(test, /dim)
  ;  n += 1
  ;endwhile
  ;extensions = extensions(0:n-1)
  ;;print, extensions
  ;ind = where( extensions eq 'signal', c )
  ;if c gt 0 then begin
  ;  print, 'Could it be that the inputs are frames, not cubes?'
  ;  retall
  ;endif
  
  ;if keyword_set(cube_nod) then begin
  ;  !ERROR_STATE.CODE=0
  ;  extensions_nod = strarr(25)
  ;  n = 0
  ;  while !ERROR_STATE.CODE eq 0 do begin
  ;    test = readfits(dir+'/'+cube_nod+'1.fits', h, exten_no=n, /silent)
  ;    extensions_nod(n) = strlowcase(strtrim(sxpar(h, 'EXTNAME')))
  ;    n += 1
  ;  endwhile
  ;  extensions_nod = extensions_nod(0:n-1)
  ;endif
  
  ;###########################################
  ;POSITION ANGLE
  tmp = readfits(dir+'/'+cube_name+a2str(1)+'.fits', cube_header, /silent)
  ;posangle = sxpar(h, 'POSANGLE', count = c)
  
  ;###########################################
  ;HOW MANY CYCLES
  spawn, '\find ' + dir + ' -maxdepth 1 -name "' + cube_name+'?.fits" -o -name "' + cube_name+'??.fits"', res
  n_files = n_elements(res)
  if (n_files gt n_cubes) and not keyword_set(combine_rasters) then begin
    str = ['/!\ The number of files (' + cube_name+'?.fits) is larger than the input number of cubes but the /combine_rasters keyword was not set, is that normal?']
    write_message, str, lun
  endif
  ;n_cycles = n_elements(res) / n_rasters
  ;if n_elements(res) gt n_rasters then n_cycles = 1 ;in case someone calls it with only a subset of available raster positions ;this causes a bug
  
  ;str = ['Number of cycles: ' + a2str(n_cycles)]
  ;write_message, str, lun
  ;cube_params.n_cycles = n_cycles
  
  
  ;==========================================================================================
  ;==========================================================================================
  ;==========================================================================================
  ;POPULATES RA AND DEC CUBES
  str = ['Building RA and DEC cubes']
  write_message, str, lun
  str = '   - Fraction of readouts >1" from the median pointing (if any): '
  write_message, str, lun
  for n = 1, n_cubes do begin ;loop on cubes
    file = dir+'/'+cube_name+a2str(n)+'.fits'
    ;file = dir+'/'+cube_name+a2str(1+n_cycles*(n-1))+'.fits'
    ra_cube = double(readfits(file, h, exten_no=getextension(file, 'ra'), /silent))
    dec_cube = double(readfits(file, h, exten_no=getextension(file, 'dec'), /silent))
    if keyword_set(cube_nod) then begin
      filenod = dir+'/'+cube_nod+a2str(n)+'.fits'
      ;file = dir+'/'+cube_nod+a2str(1+n_cycles*(n-1))+'.fits'
      ra_nod = double(readfits(filenod, h, exten_no=getextension(filenod, 'ra'), /silent))
      dec_nod = double(readfits(filenod, h, exten_no=getextension(filenod, 'dec'), /silent))
    endif
    
    for x = 0, 4 do for y = 0, 4 do begin ;loop on spatial pixels
      ;mid = n_elements(ra_cube(x, y, *)) / 2
      ind = where( cube_lres.spaxel_X eq (x+1) and cube_lres.spaxel_Y eq (y+1) and cube_lres.raster eq n )
      cube_groups[ind].filename = file
      if keyword_set(cube_nod) then begin ;if 2 nods, take the average RA, DEC
        ;cube_lres(ind).ra = 0.5 * ( double(0.5 * total( minmax(ra_cube(x, y, *), /nan) )) + double(0.5 * total( minmax(ra_nod(x, y, *), /nan) )) )
        ;cube_lres(ind).dec = 0.5 * ( double(0.5 * total( minmax(dec_cube(x, y, *), /nan) )) + double(0.5 * total( minmax(dec_nod(x, y, *), /nan) )) )
        ra1 = double( median(ra_cube(x, y, *)) )
        dec1 = double( median(dec_cube(x, y, *)) )
        ;now we check for outliers. We do it for each spaxel because the offsets might be due to a rotation of the footprint, not necessarilly a shift
        dtest = sqrt( abs( ra_cube(x, y, *) - ra1 )^2. + abs( dec_cube(x, y, *) - dec1 )^2. ) * 3600.
        indoutliers = where( dtest gt 1., coutliers )
        readouts = n_elements(reform(ra_cube(x, y, *)))
        fracoutliers = nint(100. * coutliers / readouts)
        str = '      - ['+(dir+'/'+cube_name+a2str(n)+'.fits') + ' / ' + a2str(x+1) + ', ' + a2str(y+1) + '] : ' + a2str(fracoutliers) + '%'
        if fracoutliers gt 1. then write_message, str, lun
        if fracoutliers gt 25. then begin
          str = ['      /!\ WARNING: '+a2str(fracoutliers)+'% of the readouts are >1" away from the median pointing (normal: 25% for cn, 0% for un): ']
          write_message, str, lun
        endif
        ra2 = double( median(ra_nod(x, y, *)) )
        dec2 = double( median(dec_nod(x, y, *)) )
        ;now we check for outliers. We do it for each spaxel because the offsets might be due to a rotation of the footprint, not necessarilly a shift
        dtest = sqrt( abs( ra_nod(x, y, *) - ra2 )^2. + abs( dec_nod(x, y, *) - dec2 )^2. ) * 3600.
        indoutliers = where( dtest gt 1., coutliers )
        readouts = n_elements(reform(ra_nod(x, y, *)))
        fracoutliers = nint(100. * coutliers / readouts)
        str = '      - ['+(dir+'/'+cube_name+a2str(n)+'.fits') + ' / ' + a2str(x+1) + ', ' + a2str(y+1) + '] : ' + a2str(fracoutliers) + '%'
        if fracoutliers gt 1. then write_message, str, lun
        if fracoutliers gt 25. then begin
          str = ['      /!\ WARNING: '+a2str(fracoutliers)+'% of the readouts are >1" away from the median pointing (normal: 25% for cn, 0% for un): ']
          write_message, str, lun
        endif
        ;check median distance between the 2 nods
        dtest = sqrt( abs( ra1-ra2 )^2. + abs( dec1-dec2 )^2. ) * 3600.
        if dtest gt 2. then begin
          str = ['      /!\ WARNING: the median coordinates of the 2 nods are '+str_round(a2str(dtest))+'" away from each other: ['+dir+'/'+cube_name+a2str(n)+'.fits / ' + a2str(x+1) + ', ' + a2str(y+1) + ']' ]
          write_message, str, lun
        endif
        cube_lres(ind).ra = avg( [ra1, ra2] )
        cube_lres(ind).dec = avg( [dec1, dec2] )
      ;do we create a mask to consider only the fraction of readouts that is well pointed?
      endif else begin
        ;cube_lres(ind).ra = double(0.5 * total( minmax(ra_cube(x, y, *), /nan) )) ;ra_cube(x, y, mid)
        ;cube_lres(ind).dec = double(0.5 * total( minmax(dec_cube(x, y, *), /nan) )) ;dec_cube(x, y, mid)
        cube_lres(ind).ra = double( median(ra_cube(x, y, *)) )
        cube_lres(ind).dec = double( median(dec_cube(x, y, *)) )
        ;now we check for outliers. We do it for each spaxel because the offsets might be due to a rotation of the footprint, not necessarilly a shift
        dtest = sqrt( abs( ra_cube(x, y, *) - cube_lres(ind).ra )^2. + abs( dec_cube(x, y, *) - cube_lres(ind).dec )^2. ) * 3600.
        indoutliers = where( dtest gt 1., coutliers )
        readouts = n_elements(reform(ra_cube(x, y, *)))
        fracoutliers = nint( 100. * coutliers / readouts )
        str = '      - ['+(dir+'/'+cube_name+a2str(n)+'.fits') + ' / ' + a2str(x+1) + ', ' + a2str(y+1) + '] : ' + a2str(fracoutliers) + '%'
        if fracoutliers gt 1. then write_message, str, lun
        if fracoutliers gt 25. then begin
          str = ['      /!\ WARNING: '+a2str(fracoutliers)+'% of the readouts are >1" away from the median pointing (normal: 25% for cn, 0% for un): ']
          write_message, str, lun
        endif
      endelse
      cube_spectra(ind).ra = cube_lres(ind).ra
      cube_spectra(ind).dec = cube_lres(ind).dec
    ;do we create a mask to consider only the fraction of readouts that is well pointed?
    endfor
  endfor
  if keyword_set(cube_nod) then help, temporary(ra_nod), temporary(dec_nod)
  ;writefits, dir+'/'+output_name+'_ra.fits', double(cube_lra)
  ;writefits, dir+'/'+output_name+'_dec.fits', double(cube_ldec)
  ;==========================================================================================
  ;==========================================================================================
  ;==========================================================================================
  
  
  
  
  
  
  ;==========================================================================================
  ;==========================================================================================
  ;==========================================================================================
  if keyword_set(combine_rasters) then begin
    str = 'Concatenating raster positions that are spatially coinciding...'
    write_message, str, lun
    ;---------------------------------------------------------------------------------
    ;CONCATENE SPATIALLY-COINCIDING RASTERS?
    ;get central spaxels for each raster position and calculate which cubes are co-spatial
    ind_centralspaxels = where( cube_lres.spaxel_X eq 3 and cube_lres.spaxel_Y eq 3, count_centralspaxels )
    group, cube_lres[ind_centralspaxels].ra, cube_lres[ind_centralspaxels].dec, 1./3600., ngroup
    n_rasters = max(ngroup) + 1
    ;cube_fitparams is taken care of after the for loop
    cube_lres_old = cube_lres
    cube_lres = cube_lres_old[0:25*n_rasters-1]
    cube_spectra_old = cube_spectra
    cube_spectra = cube_spectra_old[0:25*n_rasters-1]
    for i = 0, n_rasters - 1 do begin
      ind_ingroup = where( ngroup eq i, c_ingroup )
      ;update cube_groups
      str = ['Group ' + a2str(i+1) + ': ' + cube_groups[ind_centralspaxels[ind_ingroup]].filename]
      write_message, str, lun
      for j = 0, c_ingroup - 1 do cube_groups[ind_centralspaxels[ind_ingroup[j]]-12:ind_centralspaxels[ind_ingroup[j]]+12].group = i+1
      ind_ingroup = ind_ingroup[0] ;take only first cycle for populating spaxel coordinates
      cube_lres[i*25:i*25+25-1] = cube_lres_old[ind_centralspaxels[ind_ingroup]-12:ind_centralspaxels[ind_ingroup]+12]
      ;update raster (group) number?
      cube_lres[i*25:i*25+25-1].raster = i+1
      cube_spectra[i*25:i*25+25-1] = cube_spectra_old[ind_centralspaxels[ind_ingroup]-12:ind_centralspaxels[ind_ingroup]+12]
    endfor
    help, temporary(cube_lres_old), temporary(cube_spectra_old)
    cube_fitparams = cube_fitparams[0:25*n_rasters-1] ;not yet populated so we can just trim
    str = 'Number of independent groups (i.e., raster positions): ' + a2str(n_rasters)
    write_message, str, lun
    n_cycles = float(n_cubes)/float(n_rasters)
    if n_cycles ne nint(n_cycles) then begin
      str = ['/!\ The number of cycles is not an integer, pointing problems?', '/!\ The spatially-coinciding raster positions will still be combined']
      write_message, str, lun
    endif else begin
      str = 'Number of cycles                                     : ' + a2str(n_cubes/n_rasters)
      write_message, str, lun
    endelse
  endif else n_rasters = n_cubes
  cube_params.n_rasters = n_rasters
  ;==========================================================================================
  ;==========================================================================================
  ;==========================================================================================
  
  
  
  ;###########################################
  ;###########################################
  ;###########################################
  ;###########################################
  ;###########################################
  ;RASTERS LOOP
  ;###########################################
  ;###########################################
  ;###########################################
  ;###########################################
  ;###########################################
  ;###########################################
  for n = 1, n_rasters do begin ;loop on rasters
  
    if n_rasters eq 1 then raster_str = '' else raster_str = '_raster' + a2str(n)
    
    str = ['___________________________________', '___________________________________', 'Analyzing raster '+a2str(n)]
    write_message, str, lun
    
    ;tmp = readfits(dir+'/'+cube_name+a2str(n)+'.fits', h, exten_no=(where(extensions eq 'mask'))[0], /silent)
    
    if keyword_set(rebinned) then begin
       extflux = 'image'
       extwave = 'wavegrid'
       extnoise = 'stddev'
    endif else begin
       extflux = 'flux'
       extwave = 'wave'
       extnoise = 'noise'
    endelse
    ;extflux = 'image' else extflux = 'flux'
    ;if keyword_set(rebinned) then extwave = 'wavegrid' else extwave = 'wave'
    
    
    ;==========================================================================================
    ;==========================================================================================
    ;==========================================================================================
    str = 'Rearranging cubes for raster '+a2str(n)+'...'
    write_message, str, lun
    ;CONCATENE SPATIALLY-COINCIDING RASTERS?
    if keyword_set(combine_rasters) then begin
      ind_ingroup = where( cube_groups.group eq n, c_ingroup )
      filenames = cube_groups[ind_ingroup].filename
      filenames = filenames[uniq(filenames)]
      file = filenames[0]
    endif else file = dir+'/'+cube_name+a2str(n)+'.fits'
    str = '...' + file
    write_message, str, lun
    ;flux
    flux_cube = double(readfits(file, h, exten_no=getextension(file, extflux), /silent))
    ;readouts
    readouts = n_elements(reform(flux_cube(0, 0, *)))
    ;lambda
    lambda_cube = double(readfits(file, h, exten_no=getextension(file, extwave), /silent))
    if keyword_set(rebinned) then begin
      tmp = fltarr(5, 5, readouts)
      for k=0,4 do for p=0,4 do tmp(k,p,*) = lambda_cube
      lambda_cube = tmp
    endif
    ;noise
    nextnoise = getextension(file, extnoise)
    if nextnoise ne -1 then noise_cube = double(readfits(file, h, exten_no=nextnoise, /silent)) else noise_cube = 0. * flux_cube ;no noise extension, replaced by fittedflux?
    
    ;CONCATENE SPATIALLY-COINCIDING RASTERS?
    if keyword_set(combine_rasters) then begin
      if n_elements(filenames) gt 1 then begin
        ;print, 'Concatenating cubes at same position'
        for j = 1, n_elements(filenames) - 1 do begin
          file = filenames[j]
          str = '......' + file
          write_message, str, lun
          ;flux
          tmp = double(readfits(file, h, exten_no=getextension(file, extflux), /silent))
          flux_cube = [[[flux_cube]], [[tmp]]]
          ;readouts
          readouts = n_elements(reform(flux_cube(0, 0, *)))
          ;lambda
          lambda_tmp = double(readfits(file, h, exten_no=getextension(file, extwave), /silent))
          if keyword_set(rebinned) then begin
            tmp = fltarr(5, 5, readouts)
            for k=0,4 do for p=0,4 do tmp(k,p,*) = lambda_tmp
            lambda_tmp = tmp
          endif
          lambda_cube = [[[lambda_cube]], [[lambda_tmp]]]
          ;noise
          nextnoise = getextension(file, extnoise)
          if nextnoise ne -1 then begin
            tmp = double(readfits(file, h, exten_no=nextnoise, /silent))
            noise_cube = [[[noise_cube]], [[tmp]]]
          endif else noise_cube = 0. * flux_cube
        endfor
      endif
    endif
    ;==========================================================================================
    ;==========================================================================================
    ;==========================================================================================
    
    
    
    ;==========================================================================================
    ;==========================================================================================
    ;==========================================================================================
    ;COMBINE NODS
    if keyword_set(cube_nod) then begin ;combine flux and errors between the 2 nods if needed
    
      str = 'Rearranging cubes for raster '+a2str(n)+'... [nod]'
      write_message, str, lun
      
      ;CONCATENE SPATIALLY-COINCIDING RASTERS?
      if keyword_set(combine_rasters) then begin
        ind_ingroup = where( cube_groups.group eq n, c_ingroup )
        filenames = strtrans(cube_groups[ind_ingroup].filename, cube_name, cube_nod)
        filenames = filenames[uniq(filenames)]
        file = filenames[0]
      endif else file = dir+'/'+cube_nod+a2str(n)+'.fits'
      str = '...' + file
      write_message, str, lun
      ;flux
      flux_nod = double(readfits(file, h, exten_no=getextension(file, extflux), /silent))
      ;readouts
      readouts_nod = n_elements(reform(flux_nod(0, 0, *)))
      ;lambda
      lambda_nod = double(readfits(file, h, exten_no=getextension(file, extwave), /silent))
      if keyword_set(rebinned) then begin
        tmp = fltarr(5, 5, readouts_nod)
        for k=0,4 do for p=0,4 do tmp(k,p,*) = lambda_nod
        lambda_nod = tmp
      endif
      ;noise
      nextnoise = getextension(file, extnoise)
      if nextnoise ne -1 then noise_nod = double(readfits(file, h, exten_no=nextnoise, /silent)) else noise_nod = 0. * flux_nod
      
      ;CONCATENE SPATIALLY-COINCIDING RASTERS?
      if keyword_set(combine_rasters) then begin
        if n_elements(filenames) gt 1 then begin
          ;print, 'Concatenating cubes at same position [nod]'
          for j = 1, n_elements(filenames) - 1 do begin
            file = strtrans(filenames[j], cube_name, cube_nod)
            str = '......' + file
            write_message, str, lun
            ;flux
            tmp = double(readfits(file, h, exten_no=getextension(file, extflux), /silent))
            flux_nod = [[[flux_nod]], [[tmp]]]
            ;readouts
            readouts_nod = n_elements(reform(flux_cube(0, 0, *)))
            ;lambda
            lambda_tmp = double(readfits(file, h, exten_no=getextension(file, extwave), /silent))
            if keyword_set(rebinned) then begin
              tmp = fltarr(5, 5, readouts_nod)
              for k=0,4 do for p=0,4 do tmp(k,p,*) = lambda_tmp
              lambda_tmp = tmp
            endif
            lambda_nod = [[[lambda_nod]], [[lambda_tmp]]]
            ;noise
            nextnoise = getextension(file, extnoise)
            if nextnoise ne -1 then begin
              tmp = double(readfits(file, h, exten_no=nextnoise, /silent))
              noise_nod = [[[noise_nod]], [[tmp]]]
            endif else noise_nod = 0. * flux_nod
          endfor
        endif
      endif
      if size(tmp, /tname) ne 'UNDEFINED' then help, temporary(tmp)
      if size(lambda_tmp, /tname) ne 'UNDEFINED' then help, temporary(lambda_tmp)
      
      test = max(abs(lambda_nod-lambda_cube))
      if test gt 1.5e-3 then write_message, " /!\ There are differences in the wavelength elements nodA/nodB. Maximum difference in [um]: " + a2str(test), lun
      if n_elements(flux_cube) ne n_elements(flux_nod) then begin
        str = " nod A and nod B don't have the same number of elements !"
        write_message, str, lun
        stop
      endif
      flux_cube = 0.5 * ( flux_cube + flux_nod(*, *, 0:n_elements(flux_cube(0, 0, *))-1) )
      noise_cube = 0.5 * sqrt( noise_cube^2. + noise_nod^2. )
      help, temporary(flux_nod), temporary(lambda_nod), temporary(noise_nod)
    endif ;combine nods
    ;==========================================================================================
    ;==========================================================================================
    ;==========================================================================================
    
    
    
    ;==========================================================================================
    ;==========================================================================================
    ;==========================================================================================
    ;pre-flight performance correction factor?
    if keyword_set(no_correction_factor) then begin
      str = "no_correction_factor keyword is obsolete, see help"
      write_message, str, lun
      close, /all
      retall
    endif
    if keyword_set(use_correction_factor) then begin
      str = "Using correction factor (1.1 for blue, 1.3 for red)"
      write_message, str, lun
      case band of
        'blue': cfact = 1.3
        'red': cfact = 1.1
      endcase
      flux_cube = flux_cube / cfact
      noise_cube = noise_cube / cfact
    endif
    ;==========================================================================================
    ;==========================================================================================
    ;==========================================================================================
    
    
    
    
    ;deredshift wavelength scale
    lambda_cube = (lambda_cube - lambda_obs*lineshift/cspeed) / ( 1. + redshift )
    fwhm_um = lambda_rest * fwhm_kms/cspeed ;um
    
    ;readouts = n_elements(reform(flux_cube(0, 0, *)))
    str = ['Number of data points: '+a2str(25*readouts)]
    write_message, str, lun
    
    
    
    
    ;==========================================================================================
    ;==========================================================================================
    ;==========================================================================================
    ;MASK CREATION
    if n_rasters gt 1 then spawn, '\rm -f ' + dir+'/'+output_name+'_mask.fits' ;old files
    
    if file_test(dir+'/'+output_name+raster_str+'_mask.fits') and not keyword_set(cleanmask) then begin
    
      str = ['Using the mask file: ' + dir+'/'+output_name+raster_str+'_mask.fits', 'Remove this file if you want the mask to be recalculated']
      write_message, str, lun
      mask_cube = readfits(dir+'/'+output_name+raster_str+'_mask.fits', /sil)
      
    endif else begin
    
      str = ['Building mask cube']
      write_message, str, lun
      
      ;mask is 1 everywhere by default
      mask_cube = nint(flux_cube) * 0 + 1
      
      if keyword_set(rebinned) then mask_names = ['flag'] else mask_names = ['blindpixels', 'saturation', 'rawsaturation', 'noisypixels', 'badpixels', 'gratmove', 'glitch', 'outliers'] ;, 'uncleanchop']
      
      for mm = 0, n_elements(mask_names) - 1 do begin
      
        str = 'EXTENSION: ' + strtrim(mask_names(mm))
        write_message, str, lun
        ;extension was found
        if keyword_set(combine_rasters) then begin
          ind_ingroup = where( cube_groups.group eq n, c_ingroup )
          filenames = cube_groups[ind_ingroup].filename
          filenames = filenames[uniq(filenames)]
          file = filenames[0]
        endif else file = dir+'/'+cube_name+a2str(n)+'.fits'
        
        next = getextension(file, strtrim(mask_names(mm)))
        if next ne -1 then begin
        
          str = '...' + file
          write_message, str, lun
          
          mask_i_cube = long(readfits(file, h, exten_no=next, /silent))
          mask_i_cube = convert2bool(mask_i_cube)
          
          
          ;CONCATENE SPATIALLY-COINCIDING RASTERS?
          if keyword_set(combine_rasters) then begin
            if n_elements(filenames) gt 1 then begin
              ;print, 'Concatenating masks at same position'
              for j = 1, n_elements(filenames) - 1 do begin
                file = filenames[j]
                str = '......' + file
                write_message, str, lun
                next = getextension(file, strtrim(mask_names(mm)))
                tmp = long(readfits(file, h, exten_no=next, /silent))
                tmp = convert2bool(tmp)
                mask_i_cube = [[[mask_i_cube]], [[tmp]]]
              endfor
            endif
          endif
          ind2 = where( mask_i_cube eq 1, c)
          if c gt 0 then mask_cube(ind2) = 0
          
          strcount = '         ' + strtrim(mask_names(mm))+' (all spaxels): '+a2str(c)
          if mask_names(mm) eq 'saturation' and c gt 0 then write_message, '/!\ WARNING: some ramps are saturated !', lun
          
          ;MASK OTHER NOD
          if keyword_set(cube_nod) then begin ;create mask for the other nod if needed. If the data is bad in the other nod we cannot use it
            ;CONCATENE SPATIALLY-COINCIDING RASTERS?
            if keyword_set(combine_rasters) then begin
              ind_ingroup = where( cube_groups.group eq n, c_ingroup )
              filenames = strtrans(cube_groups[ind_ingroup].filename, cube_name, cube_nod)
              filenames = filenames[uniq(filenames)]
              file = filenames[0]
            endif else file = dir+'/'+cube_nod+a2str(n)+'.fits'
            str = '...' + file
            write_message, str, lun
            mask_i_nod = long(readfits(file, h, exten_no=getextension(file, strtrim(mask_names(mm))), /silent))
            mask_i_nod = convert2bool(mask_i_nod)
            
            ;CONCATENE SPATIALLY-COINCIDING RASTERS?
            if keyword_set(combine_rasters) then begin
              if n_elements(filenames) gt 1 then begin
                ;print, 'Concatenating masks at same position [nod]'
                for j = 1, n_elements(filenames) - 1 do begin
                  file = filenames[j]
                  str = '......' + file
                  write_message, str, lun
                  next = getextension(file, strtrim(mask_names(mm)))
                  tmp = long(readfits(file, h, exten_no=next, /silent))
                  tmp = convert2bool(tmp)
                  mask_i_nod = [[[mask_i_nod]], [[tmp]]]
                endfor
              endif
            endif
            ind2 = where( mask_i_nod eq 1, c)
            if c gt 0 then mask_cube(ind2) = 0
            
            if size(tmp, /tname) ne 'UNDEFINED' then help, temporary(tmp)
            
            strcount = strcount + ' & '+a2str(c)
            if mask_names(mm) eq 'saturation' and c gt 0 then write_message, '/!\ WARNING: some ramps are saturated !', lun
          endif
          
          write_message, strcount, lun
          
          ;print, 'MASK check: ', n_elements(mask_i_cube[0,0,*]), n_elements(lambda_cube[0,0,*])
          if n_elements(mask_i_cube[0, 0, *]) ne n_elements(lambda_cube[0, 0, *]) then begin
            str = '/!\ The mask has a wrong number of readouts!'
            write_message, str, lun
            ;close, /all
            ;retall
          endif
        ;if strtrim(mask_names(mm)) eq 'uncleanchop' then stop
          
        endif ;extension not found
        
      endfor
      help, temporary(mask_i_cube)
      if keyword_set(cube_nod) then help, temporary(mask_i_nod)
      writefits, dir+'/'+output_name+raster_str+'_mask.fits', mask_cube
      str = ['Mask cube written in '+dir+'/'+output_name+raster_str+'_mask.fits']
      write_message, str, lun
      
    endelse
    ;==========================================================================================
    ;==========================================================================================
    ;==========================================================================================
    
    
    
    
    if keyword_set(savedatacloud) and not keyword_set(central3by3combine) and not keyword_set(central5by5combine) then begin
      datacloud = replicate({lambda: fltarr(n_elements(lambda_cube(0, 0, *))), flux: fltarr(n_elements(lambda_cube(0, 0, *))), error: fltarr(n_elements(lambda_cube(0, 0, *))), $
        mask: intarr(n_elements(lambda_cube(0, 0, *))), ra: fltarr(n_elements(lambda_cube(0, 0, *))), dec: fltarr(n_elements(lambda_cube(0, 0, *))), $
        spaxel_X: 0, spaxel_Y: 0, lambda_obs: fltarr(n_elements(lambda_cube(0, 0, *)))}, 5*5)
      k = 0
      for x = 0, 4 do for y = 0, 4 do begin
        datacloud(k).spaxel_X = x+1
        datacloud(k).spaxel_Y = y+1
        datacloud(k).lambda = lambda_cube(x, y, *)
        datacloud(k).lambda_obs = lambda_cube(x, y, *) * (1. + redshift) + lambda_obs*lineshift/cspeed
        datacloud(k).flux = flux_cube(x, y, *)
        datacloud(k).error = noise_cube(x, y, *)
        datacloud(k).mask = mask_cube(x, y, *)
        datacloud(k).ra = ra_cube(x, y, *)
        datacloud(k).dec = dec_cube(x, y, *)
        k += 1
      endfor
      write_message, 'Saving data cloud in ' + a2str(dir+'/'+output_name+raster_str+'_datacloud.sav'), lun
      save, filename = dir+'/'+output_name+raster_str+'_datacloud.sav', datacloud
    endif
    
    lambda_cube -= lambda_rest    ;lambda_rest centered, easier to read the plots
    ;lambda_raw -= lambda_rest
    
    ;save reference wavelength scale for common scale for smoothing later. Ref is central spaxel of raster 1
    if n eq 1 then lambda_ref = lambda_cube(2, 2, *)
    
    
    
    
    ;==========================================================================================
    ;==========================================================================================
    ;==========================================================================================
    ;COMBINE SPAXELS / FLUX CALIBRATION
    
    ;PLS correction factor
    if keyword_set(central3by3combine) or keyword_set(central5by5combine) then begin
      write_message, 'Combining spaxels...', lun
      
      if n_elements(combine) eq 2 then begin
        indx = combine(0)
        indy = combine(1)
      endif else begin
        if file_test(dir+'/footprint_results.sav') then begin
          write_message, 'Using results from footprint_results.sav', lun
          restore, dir + '/footprint_results.sav'
          indx = result.brightest(0)-1
          indy = result.brightest(1)-1
          write_message, 'Brightest spaxel: ' + a2str(result.brightest(0)) + ', ' + a2str(result.brightest(1)), lun
        endif else begin
          indx = 2
          indy = 2
        endelse
      endelse
      if keyword_set(central5by5combine) then begin
        indx = 2
        indy = 2
      endif
      ;indx=2
      ;indy=2
      print, 'Combination around spaxel: indx=',indx+1, ' - indy=',indy+1
      
      
      set_plot, 'z'
      ERASE
      device, set_resolution=reso, set_pixel_depth = 24, decomposed = 0
      
      ;following just to get lambda_smooth
      bin = fwhm_um / sampling
      
      ind = where( mask_cube(indx, indy, *) ge 0, c )
      spec_smooth, lambda_cube(indx, indy, ind), flux_cube(indx, indy, ind), lambda_smooth, flux_smooth, bin = bin, lambda_ref = lambda_ref, $
        line_range = constraints_continuum(0)*fwhm_um, getbin = getbin, lun = lun
      lambda_comb = lambda_smooth
      ;print, lambda_smooth
      flux_comb = flux_smooth * 0.
      error_comb = flux_smooth * 0.
      
      set_plot, 'PS'
      !p.font = 0
      !p.position = [0.12, 0.15, 0.95, 0.85]
      if keyword_set(central5by5combine) then begin
        spawn, '\rm -f '+dir+'/combinedspectra_5x5'+raster_str+'.png'
        device, filename=dir+'/combinedspectra_5x5'+raster_str+'.eps', /encaps, /color, /helvetica
      endif else begin
        spawn, '\rm -f '+dir+'/combinedspectra_3x3'+raster_str+'.png'
        device, filename=dir+'/combinedspectra_3x3'+raster_str+'.eps', /encaps, /color, /helvetica
      endelse
      loadct, 0, /sil
      plot, indgen(10), /nodata, yran = [-max(flux_smooth, /nan), 2.5*max(flux_smooth, /nan)], xran = minmax(lambda_smooth, /nan), $
        ytit = 'Signal [Jy]', xtit = textoidl('\lambda')+' ['+textoidl('\mu')+'m]'
      plots, !x.crange, [0., 0.], line = 1
      loadct, 38, /sil
      
      ;save the 3x3 combined cube results
      if keyword_set(align) then begin
        save, filename = dir+'/'+output_name+'.sav', cube_lres, cube_fitparams, cube_spectra, cube_params, cube_header
        if keyword_set(central5by5combine) then restore, dir+'/'+strtrans(output_name, '_5x5', '')+'.sav' $
        else restore, dir+'/'+strtrans(output_name, '_3x3', '')+'.sav'
      endif
      
      ;combine
      ;combcube = dblarr(n_elements(lambda_comb), 9)
      ;k = 0
      if keyword_set(central5by5combine) then nbcomb = 2 else nbcomb = 1
      for x = indx-nbcomb, indx+nbcomb do for y = indy-nbcomb, indy+nbcomb do begin
        if x lt 0 or y lt 0 or x gt 4 or y gt 4 then continue
        ind = where( mask_cube(x, y, *) ge 0, c )
        lambda_in = lambda_cube(x, y, ind)
        if keyword_set(align) then begin
          ind_cube = where( cube_lres.spaxel_X eq (x+1) and cube_lres.spaxel_Y eq (y+1) and cube_lres.raster eq n )
          lambda_in = (lambda_in+lambda_rest) * (1. - cube_lres(ind_cube).velocity / cspeed) - lambda_rest
        endif
        spec_smooth, lambda_in, flux_cube(x, y, ind), lambda_smooth, flux_smooth, error_smooth, bin = bin, lambda_ref = lambda_ref, line_range = constraints_continuum(0)*fwhm_um
        
        ;---   remove a continuum before adding up the rebinned spectra
        ;first guess of the polynomial around the line
        continuum_ind = where( (abs(lambda_smooth) gt constraints_continuum(0)*fwhm_um) and (lambda_smooth gt 0.-constraints_continuum(1)*fwhm_um) and (lambda_smooth lt constraints_continuum(2)*fwhm_um), c, complement = line_ind )
        if c eq 0 then stop
        
        ;polynomial
        parinfo = replicate({value:0.D, fixed:0, limited:[0, 0], limits:[0.D, 0]}, 3 + poly_degree+1)
        parinfo[3].value = median(flux_smooth(continuum_ind))
        parinfo[1].value = 1. ;cannot be zero
        parinfo[0:2].fixed = replicate(1, 3)
        parinfo[0].value = 0.
        parinfo[2].value = 0.
        
        ;make the fit
        weights = 1./error_smooth(continuum_ind)
        params_cont = mpfitfun('ContinuumModel', lambda_smooth(continuum_ind), flux_smooth(continuum_ind), weights = weights, parinfo = parinfo, /nan, /quiet, yfit = yfit_cont, STATUS=status, ERRMSG=errmsg, functargs={poly_degree: poly_degree})
        if status LE 0 then message, errmsg
        
        ;remove the continuum
        continuum_smooth = interpol(yfit_cont, lambda_smooth(continuum_ind), lambda_smooth)
        flux_smooth -= continuum_smooth
        ;---
        
        ;plot smoothed spectrum
        oplot, lambda_smooth, flux_smooth, color = 25 + x*25 + y*25
        flux_comb += interpol(flux_smooth, lambda_smooth, lambda_comb) ;should be on the same grid already
        error_comb += interpol(error_smooth^2, lambda_smooth, lambda_comb)
      endfor
      
      error_comb = sqrt(error_comb)
      oplot, lambda_comb, flux_comb, ps = 10, thick = 4
      ;print, flux_comb
      ;write_png, dir+'/combinedspectra_3x3.png', tvrd(/true)
      device, /close
      ind = where( finite(flux_comb), c )
      
      ;ind2 = where( (lambda_cube(indx, indy, *) gt min(lambda_comb(ind), /nan)) and (lambda_cube(indx, indy, *) lt max(lambda_comb(ind), /nan)), c, complement = noind2 )
      ;flux_cube(indx, indy, ind2) = interpol(flux_comb(ind), lambda_comb(ind), lambda_cube(indx, indy, ind2))
      ;noise_cube(indx, indy, ind2) = interpol(error_comb(ind), lambda_comb(ind), lambda_cube(indx, indy, ind2))
      ;if noind2(0) ne -1 then begin
      ;  flux_cube(indx, indy, noind2) = !values.f_nan
      ;  noise_cube(indx, indy, noind2) = !values.f_nan
      ;endif
      ;-- use directly the rebinned cube and not the reconstructed cloud!
      flux_cube = fltarr(5,5,n_elements(lambda_comb))
      noise_cube = fltarr(5,5,n_elements(lambda_comb))
      lambda_cube = fltarr(5,5,n_elements(lambda_comb))
      mask_cube = fltarr(5,5,n_elements(lambda_comb))+1 ;keep all datapoints
      flux_cube[indx, indy, *] = flux_comb
      noise_cube[indx, indy, *] = error_comb
      lambda_cube[indx, indy, *] = lambda_comb
      
      ;restore the 3x3 combined cube result
      if keyword_set(align) then restore, dir+'/'+output_name+'.sav' ;, cube_lres, cube_fitparams, cube_spectra, cube_params, cube_header
      
      if keyword_set(central5by5combine) then begin
        pls_corr = 1./interpol(pls_frac_in5x5, pls_lambda, (lambda_cube(indx, indy, *)+lambda_rest)*(1.+redshift))
        write_message, 'Point-like source correction (5x5)~ '+a2str(avg(pls_corr, /nan)), lun
      endif else begin
        pls_corr = 1./interpol(pls_frac_in3x3, pls_lambda, (lambda_cube(indx, indy, *)+lambda_rest)*(1.+redshift))
        write_message, 'Point-like source correction (3x3)~ '+a2str(avg(pls_corr, /nan)), lun
      endelse
      flux_cube(indx, indy, *) *= pls_corr
      noise_cube(indx, indy, *) *= pls_corr
      indx = [indx, indx]
      indy = [indy, indy]
    endif else begin ;pls_corr keyword was set? but not central3x3combine
      indx = [0, 4]
      indy = [0, 4]
      case pls_corr of
        1: begin
          pls_corr = 1./interpol(pls_frac_in1, pls_lambda, (lambda_cube+lambda_rest)*(1.+redshift))
          ;pls_corr = pls_corr[0]
          write_message, 'Point-like source correction (central spaxel)~ '+a2str(avg(pls_corr, /nan)), lun
          flux_cube *= pls_corr
          noise_cube *= pls_corr
        end
        2: begin
          pls_corr = 1./interpol(pls_frac_in3x3, pls_lambda, (lambda_cube+lambda_rest)*(1.+redshift))
          ;pls_corr = pls_corr[0]
          write_message, 'Point-like source correction (3x3)~ '+a2str(avg(pls_corr, /nan)), lun
          flux_cube *= pls_corr
          noise_cube *= pls_corr
        end
        3: begin
          pls_corr = 1./interpol(pls_frac_in5x5, pls_lambda, (lambda_cube+lambda_rest)*(1.+redshift))
          ;pls_corr = pls_corr[0]
          write_message, 'Point-like source correction (5x5)~ '+a2str(avg(pls_corr, /nan)), lun
          flux_cube *= pls_corr
          noise_cube *= pls_corr
        end
        else:
      endcase
    endelse
    ;==========================================================================================
    ;==========================================================================================
    ;==========================================================================================
    
    
    
    
    ;###########################################
    ;###########################################
    ;###########################################
    ;###########################################
    ;###########################################
    ;SPAXEL LOOP
    ;###########################################
    ;###########################################
    ;###########################################
    ;###########################################
    ;###########################################
    ;###########################################
    for x = indx(0), indx(1) do for y = indy(0), indy(1) do begin ;loop on spatial pixels, or 3x3 around brightest spaxel if spaxels are combined
    
      ;if x ne 1 or y ne 2 then continue
      ;if x ne 0 or y ne 1 then continue
    
      ind_cube = where( cube_lres.spaxel_X eq (x+1) and cube_lres.spaxel_Y eq (y+1) and cube_lres.raster eq n )
      
      str = ['___________________________________', 'Raster '+a2str(n)+', Pixel '+a2str(x+1)+' '+a2str(y+1)]
      write_message, str, lun
      
      if not keyword_set(skipfit) then begin
      
        flux = reform(flux_cube(x, y, *))
        lambda = reform(lambda_cube(x, y, *))
        noise = reform(noise_cube(x, y, *))
        mask = reform(mask_cube(x, y, *))
        
        ;sort by wavelengths
        ind = sort(lambda)
        lambda = lambda(ind)
        flux = flux(ind)
        noise = noise(ind)
        mask = mask(ind)
        
        ;keep all data points to plot flagged data afterwards
        lambda_raw = lambda
        flux_raw = flux
        
        
        ;==========================================================================================
        ;==========================================================================================
        ;==========================================================================================
        ;ADAPT CENTROID POSITION BASED ON VELOCITY MAP?
        ;centroid is around the max
        ;centroid = lambda( (where( flux_smooth eq max(flux_smooth) ))(0) )
        ;if abs( centroid ) gt 3.*fwhm then centroid = 0. ;takes 0 if guess is too far away
        centroid = 0.        ;we know where the line is supposed to be
        if keyword_set(use_velomap) then begin
          write_message, ['Using smoothed velocity map to refine 1st guess of line centroid'], lun
          if size(use_velomap, /tname) eq 'INT' then begin ;the user didn't provide a file
            spawn, '\ls ' + dir+'/*_LINE'+line+'_Flux.fits', files
            velofile = files(0)
            if n_elements(files) gt 1 then begin
              str = 'Several _Flux.fits were found, aborting...'
              write_message, str, lun
              close, /all
              retall
            endif
          endif else begin
            if file_test(use_velomap) then velofile = use_velomap else begin
              write_message, ['/!\ Velocity map could not be found'], lun
              close, /all
              retall
            endelse
          endelse
          tmp = readfits(velofile, h)
          velomap = reform(tmp(*, *, 2))
          velomap = filter_image(velomap, median = 9, /all) ;the median box should be changed depending on the map size and S/N
          extast, h, astr
          ad2xy, cube_lres(ind_cube).ra, cube_lres(ind_cube).dec, astr, velo_x, velo_y
          velo_x = floor(velo_x)
          velo_y = floor(velo_y)
          velo_x = max([velo_x, 0])
          velo_y = max([velo_y, 0])
          velo_x = min([velo_x, n_elements(velomap(*, 0))-1])
          velo_y = min([velo_y, n_elements(velomap(0, *))-1])
          centroid = (velomap(velo_x, velo_y)/cspeed) * lambda_rest
        ;print, velo_x, velo_y, centroid
        endif else centroid = 0.
        ;==========================================================================================
        ;==========================================================================================
        ;==========================================================================================
        
        
        
        
        ;select only good data from now on
        ind = where( finite(flux) and finite(lambda) and (mask eq 1), count, complement = flagged )
        if count eq 0 then begin
          str = ['No valid data point was found']
          write_message, str, lun
        endif
        lambda = lambda(ind)
        flux = flux(ind)
        noise = noise(ind)
        mask = mask(ind)
        
        ;###########################################
        ;keep points for continuum and make a
        ;smooth spectrum by rebinning
        
        bin = fwhm_um / sampling
        
        if keyword_set(rebinned) or keyword_set(central3by3combine) or keyword_set(central5by5combine) then begin
          lambda_smooth = lambda ;; already rebinned
          flux_smooth = flux
          error_smooth = noise
          error_indiv = noise
          goodinds = indgen(n_elements(lambda))
        ;ninbins= readouts
        endif else spec_smooth, lambda, flux, lambda_smooth, flux_smooth, error_smooth, bin = bin, lambda_ref = lambda_ref, error_indiv = error_indiv, $
          line_range = constraints_continuum(0)*fwhm_um, getbin = getbin, lun = lun
        ;if n eq 1 and x eq 3 and y eq 1 then stop
          
          
          
        ;============================================================================
        ;============================================================================
        ;============================================================================
        ;DO WE SUBTRACT SPECTRA OF RASTERS DERIVED FROM A PREVIOUS RUN?
        if keyword_set(subtract_rasters) then begin
          str = ['Subtracting some rasters: ', a2str(subtract_rasters)]
          ;npoints = 8000
          ;cube_spectra = replicate({spaxel_X: nint(0), spaxel_Y: nint(0), raster: nint(0), lambda: fltarr(npoints)*!values.f_nan, flux: fltarr(npoints)*!values.f_nan, $
          ;                          normalized_flux: fltarr(npoints)*!values.f_nan, cont_subtracted_flux: fltarr(npoints)*!values.f_nan, lambda_obs: fltarr(npoints)*!values.f_nan, error: fltarr(npoints)*!values.f_nan, ra: double(0.), dec: double(0.)}, 25*n_cubes)
          match2, cube_spectra_before_subtracting.raster, subtract_rasters, valrastersubtract ;tmp has the same number of elements as A (cube_spectra_before_subtracting.raster)
          indspaxelsubtract = where( valrastersubtract gt -1 and cube_spectra_before_subtracting.spaxel_X eq (x+1) and cube_spectra_before_subtracting.spaxel_Y eq (y+1), countrasters )
          lambda_subtract = cube_spectra_before_subtracting[indspaxelsubtract[0]].lambda - lambda_rest
          indok = where( finite(lambda_subtract) )
          lambda_subtract = lambda_subtract[indok]
          flux_subtract = fltarr(n_elements(lambda_subtract))
          for i = 0, n_elements(lambda_subtract) - 1 do begin
            arr = fltarr(countrasters)
            for j = 0, countrasters - 1 do begin
              tmpw = cube_spectra_before_subtracting[indspaxelsubtract[j]].lambda - lambda_rest
              tmpf = cube_spectra_before_subtracting[indspaxelsubtract[j]].flux / ( 1.e-6 / (9.4*!pi/180./3600.)^2. )
              indok = where( finite(tmpw) and finite(tmpf) )
              tmpf = interpol(tmpf[indok], tmpw[indok], lambda_subtract)
              arr[j] = tmpf[i]
            endfor
            flux_subtract[i] = median(arr, /even)
          endfor
          ;resample to keep only large-scale variations (>0.01um)
          bin_subtract = median( abs( lambda_subtract - shift(lambda_subtract, 1) ) ) ;that's the bin size in um ;we want 0.01um
          if bin_subtract lt 0.01 then begin ;if bin size is too small we smooth though a regridding
            newn = nint( n_elements(lambda_subtract) * (bin_subtract/0.01) )
            flux_subtract = congrid(flux_subtract, newn)
            lambda_subtract = congrid(lambda_subtract, newn)
          endif
          if file_test(dir+'/PLOTS_' + output_name + '/background/background_subtract_'+a2str(x+1)+a2str(y+1)+'.eps') eq 0 then begin
            set_plot, 'PS'
            device, filename = dir+'/PLOTS_' + output_name + '/background/background_subtract_'+a2str(x+1)+a2str(y+1)+'.eps', /encaps, /helvetica, xsize = 20, ysize = 15
            plot, lambda_subtract, flux_subtract, ps = 10, xsty = 1, ysty = 1, xtit = 'Wavelength [um]', ytit = 'Flux density [Jy]'
            device, /close
          endif
          flux -= interpol(flux_subtract, lambda_subtract, lambda)
          spec_smooth, lambda, flux, lambda_smooth, flux_smooth, error_smooth, bin = bin, lambda_ref = lambda_ref, error_indiv = error_indiv, $
            line_range = constraints_continuum(0)*fwhm_um, getbin = getbin, lun = lun
        endif
        ;============================================================================
        ;============================================================================
        ;============================================================================
        
        
        
        
        ;- comment because those quantities are already the same (albeit mask and finite not applied to comb)
        ;if keyword_set(central3by3combine) or keyword_set(central5by5combine) then begin
        ;  error_smooth = interpol(error_comb, lambda_comb, lambda_smooth) ;modify error array for spaxel combination
        ;endif
        
        
        cube_spectra(ind_cube).lambda(0:n_elements(lambda_smooth)-1) = lambda_smooth + lambda_rest
        cube_spectra(ind_cube).lambda_obs(0:n_elements(lambda_smooth)-1) = (lambda_smooth + lambda_rest) * (1. + redshift) + lambda_obs*lineshift/cspeed
        cube_spectra(ind_cube).flux(0:n_elements(lambda_smooth)-1) = flux_smooth * 1.e-6 / (9.4*!pi/180./3600.)^2. ;MJy sr-1
        cube_spectra(ind_cube).error(0:n_elements(lambda_smooth)-1) = error_smooth * 1.e-6 / (9.4*!pi/180./3600.)^2. ;MJy sr-1
        
        
        
        ;==========================================================================================
        ;==========================================================================================
        ;==========================================================================================
        ;sigma clipping
        if not keyword_set(central3by3combine) and not keyword_set(central5by5combine) then begin
          residual = flux - interpol(flux_smooth, lambda_smooth, lambda)
          indclip = where( abs(residual) gt n_sigma_clipping*robust_sigma(residual, /zero), c )
          ;indclip = where( abs(flux - interpol(flux_smooth, lambda_smooth, lambda)) gt n_sigma_clipping*interpol(error_smooth, lambda_smooth, lambda), c)
          if c gt 0 then begin
            str = ['Sigma-clipping ('+a2str(nint(n_sigma_clipping))+'-sigma): '+a2str(c)+' readouts (/'+a2str(n_elements(flux))+')']
            write_message, str, lun
            mask(indclip) = 0
          endif
        endif else begin
          goodinds = indgen(n_elements(lambda))
          error_indiv = error_smooth
        endelse
        lambda_flagged = lambda
        flux_flagged = flux
        ind = where( (mask eq 1), complement=flagged_sigmaclip )
        lambda = lambda(ind)
        flux = flux(ind)
        noise = noise(ind)
        mask = mask(ind)
        
        ;re-calculate the smoothed spectrum after sigma-clipping
        bin = fwhm_um / sampling
        if not keyword_set(rebinned) and not keyword_set(central3by3combine) and not keyword_set(central5by5combine) then $
          spec_smooth, lambda, flux, lambda_smooth, flux_smooth, error_smooth, bin = bin, lambda_ref = lambda_ref, $
          error_indiv = error_indiv, line_range = constraints_continuum(0)*fwhm_um, getbin = getbin, goodinds = goodinds
        ;test
        ;lambda = lambda[goodinds]
        ;flux = flux[goodinds]
        ;noise = noise[goodinds]
        ;==========================================================================================
        ;==========================================================================================
        ;==========================================================================================
          
        ;if keyword_set(rebinned) then error_indiv = error_smooth

        ;rbcube_name = dir+'/cube01.fits'
        ;if file_test(rbcube_name) then begin
        ;   stop
        ;   rbcube_f = readfits(rbcube_name, exten_no=get_extension(rbcube_name, 'image'))
        ;   rbcube_w = readfits(rbcube_name, exten_no=get_extension(rbcube_name, 'waveGrid'))
        ;   rbcube_ra = readfits(rbcube_name, exten_no='ra')
        ;   rbcube_dec = readfits(rbcube_name, exten_no='dec')
        ;endif

          
          
        ;sigma value
        ;sigma = fwhm_um / (2.*sqrt(2.*alog(2.)))
        ;fwhm_kms = cspeed * fwhm / lambda_rest
          
          
        ;==========================================================================================
        ;==========================================================================================
        ;==========================================================================================
        ;select only data around line from now on
        indfit = where( ((lambda-centroid) gt 0.-constraints_continuum(1)*fwhm_um) and ((lambda-centroid) lt constraints_continuum(2)*fwhm_um), count, complement = flagged_outsiderange )
        if count eq 0 then begin
          str = ['No data point was found around the centroid wavelength ' + a2str(centroid), 'min/max of the input wavelength array: ' + a2str(minmax(lambda, /nan))]
          write_message, str, lun
          stop
        endif
        
        lambda_fullrange = lambda
        flux_fullrange = flux
        noise_fullrange = noise
        
        lambda = lambda(indfit)
        flux = flux(indfit)
        noise = noise(indfit)
        mask = mask(indfit)
        ;==========================================================================================
        ;==========================================================================================
        ;==========================================================================================
        
        
        
        
        ;==========================================================================================
        ;==========================================================================================
        ;==========================================================================================
        ;CONTINUUM FIT
        indfit_smooth = where( ((lambda-centroid) gt 0.-constraints_continuum(1)*fwhm_um) and ((lambda-centroid) lt constraints_continuum(2)*fwhm_um), count, complement = flagged_outsiderange )
        
        ;###########################################
        ;first guess of the polynomial around the line
        continuum_ind = where( (abs(lambda-centroid) gt constraints_continuum(0)*fwhm_um) and ((lambda-centroid) gt 0.-constraints_continuum(1)*fwhm_um) and ((lambda-centroid) lt constraints_continuum(2)*fwhm_um), c, complement = line_ind )
        if c eq 0 then stop
        
        parinfo = replicate({value:0.D, fixed:0, limited:[0, 0], limits:[0.D, 0]}, 3 + poly_degree+1)
        
        ;polynomial
        parinfo[3].value = median(flux(continuum_ind))
        ;fringes ?
        if keyword_set(fringes) then begin
          parinfo[0:2].value = [1., 1.5, 0.] ;amplitude (Jy), period (um), phase
          parinfo[0].limited = [1, 1]
          parinfo[0].limits  = [0., 2.5]
          parinfo[1].limited = [1, 1]
          parinfo[1].limits  = [0.5, 2.5]
        endif else begin
          parinfo[0:2].fixed = replicate(1, 3)
          parinfo[0].value = 0.
          parinfo[1].value = 1. ;cannot be zero
          parinfo[2].value = 0.
        endelse
        
        ;parinfo[1:poly_degree].value = fltarr(poly_degree)
        ;if keyword_set(use_errors) then weights = 1./noise(continuum_ind) else weights = 1./interpol(sqrt(ninbins)*error_smooth, lambda_smooth, lambda(continuum_ind))^2. * max([1.-lambda(continuum_ind)/(3.*fwhm_um), 0]) ;+fltarr(n_elements(continuum))
        if keyword_set(use_errors) then weights = 1./noise(continuum_ind) else weights = 1./interpol(error_smooth, lambda_smooth, lambda(continuum_ind))^2. * exp(-0.01*lambda(continuum_ind)^2./fwhm_um^2.) ;* max([1.-lambda(continuum_ind)/(3.*fwhm_um), 0]) ;+fltarr(n_elements(continuum))
        
        ;noise_eff = noise(continuum)
        params_cont = mpfitfun('ContinuumModel', lambda(continuum_ind), flux(continuum_ind), weights = weights, parinfo = parinfo, /nan, /quiet, yfit = yfit_cont, STATUS=status, ERRMSG=errmsg, functargs={poly_degree: poly_degree})
        
        if status LE 0 then message, errmsg
        ;mean = moment(flux(continuum_ind)-ContinuumModel(lambda(continuum), params_cont, poly_degree = poly_degree), sdev = sdev) ;keep standard deviation
        tmp = flux(continuum_ind)-ContinuumModel(lambda(continuum_ind), params_cont, poly_degree = poly_degree)
        sdev = robust_sigma(tmp, /zero);keep standard deviation
        if sdev eq -1 then begin
          str = 'replacing by stdev value'
          write_message, str, lun
          sdev = stdev(tmp)
        endif
        ind = where( abs(flux(continuum_ind)-ContinuumModel(lambda(continuum_ind), params_cont, poly_degree = poly_degree)) gt 3.*sdev, count, complement = noind )
        if count gt 0 then weights(ind) = 0.
        if keyword_set(use_errors) then weights(noind) = 1./noise(continuum_ind(noind)) else weights(noind) = 1./interpol(error_smooth, lambda_smooth, lambda(continuum_ind(noind)))^2. * max([1.-lambda(continuum_ind(noind))/(3.*fwhm_um), 0])
        ;if count gt 0 then noise_eff(ind) = 1.e8
        parinfo[*].value = params_cont
        params_cont = mpfitfun('ContinuumModel', lambda(continuum_ind), flux(continuum_ind), weights = weights, parinfo = parinfo, /nan, /quiet, yfit = yfit_cont, STATUS=status, ERRMSG=errmsg, functargs={poly_degree: poly_degree})
        if status LE 0 then message, errmsg
        
        ;str = ['1st guess Continuum: degree 0, 1, 2, ...', $
        ;  a2str(params_cont(3:5))]
        ;write_message, str, lun
        
        continuum_smooth_ind = where( (abs(lambda_smooth-centroid) gt constraints_continuum(0)*fwhm_um) and ((lambda_smooth-centroid) gt 0.-constraints_continuum(1)*fwhm_um) and ((lambda_smooth-centroid) lt constraints_continuum(2)*fwhm_um), c, complement = line_smooth_ind )
        if c le 3 then begin
           str = 'Problem with continuum range,not enough points... use constraints_continuum option?'
           write_message, str, lun
           retall
        endif
        cube_lres(ind_cube).continuum = interpol(yfit_cont, lambda(continuum_ind), centroid)
        tmp = flux_smooth(continuum_smooth_ind) - ContinuumModel(lambda_smooth(continuum_smooth_ind), params_cont, poly_degree = poly_degree)
        sdev = robust_sigma(tmp, /zero)
        if sdev eq -1 then begin
          str = 'replacing by stdev value'
          write_message, str, lun
          sdev = stdev(tmp)
        endif
        cube_lres(ind_cube).continuum_error = sdev
        ;==========================================================================================
        ;==========================================================================================
        ;==========================================================================================
        
        
        
        
        
        ;==========================================================================================
        ;==========================================================================================
        ;==========================================================================================
        ;LINE FIT
        ;###########################################
        ;initial value of maximum from the smoothed spectrum
        maxi = 1.5 * max( flux_smooth - ContinuumModel(lambda_smooth, params_cont, poly_degree = poly_degree), /nan )
        ;print, 'maximum =', maxi
        ;if x eq 2 and y eq 2 then stop
        
        
        ;###########################################
        ;initial value of sigma
        sigma = (fwhm_kms/cspeed)*lambda_rest /(2.*sqrt(2.*alog(2.)))
        ;fwhm_kms = cspeed * fwhm / lambda_rest
        ;print, 'sigma =', sigma, fwhm
        
        ;###########################################
        ;redefine centroid here
        if keyword_set(find) then begin
          ;write_message, ['/!\ Warning: find option has been deactivated temporarily'], lun
          ;weights = 1./(2.*sqrt(ninbins)*error_smooth^2.) * max([1.-lambda_smooth/(3.*fwhm_um), 0]);1.+fltarr(n_elements(lambda))
          weights = 1./(error_indiv^2.) * max([1.-lambda_smooth/(3.*fwhm_um), 0]);1.+fltarr(n_elements(lambda))
          ind = where( finite(lambda_smooth) and finite(flux_smooth) and abs(lambda_smooth) lt find*fwhm_um )
          chisqarr = fltarr(n_elements(ind))
          for i = 0, n_elements(ind) - 1 do begin
            parinfo = replicate({value:0.D, fixed:0, limited:[0, 0], limits:[0.D, 0]}, 12+poly_degree+1)
            parinfo[0:2].value = [flux_smooth(ind(i))-median(flux_smooth(ind)), lambda_smooth(ind(i)), sqrt(sigma^2.+(broadening_sigma(1)*lambda_rest/cspeed)^2.)]
            parinfo[1:2].fixed = 1
            parinfo[indparams_contfringes].value = params_cont
            for ncomp = 0, n_elements(addcomponents) - 1 do begin
              parinfo[3+3*ncomp:3+3*ncomp+2].value = [0., 0., sigma] ;fix other components
              parinfo[3+3*ncomp:3+3*ncomp+2].fixed = replicate(1, 3)
            endfor
            params = mpfitfun('SpectrumModel', lambda_smooth, flux_smooth, weights = weights, parinfo = parinfo, /nan, /quiet, $
              yfit = yfit, STATUS=status, ERRMSG=errmsg, dof = dof, nfree = nfree, bestnorm = bestnorm, perror = perror, $
              functargs={poly_degree: poly_degree, ncomps: n_elements(addcomponents)})
            if status LE 0 then message, errmsg
            chisqarr(i) = bestnorm/dof
          endfor
          indbest = where( chisqarr eq min(chisqarr, /nan) )
          str = [fwhm_um, lambda_smooth(ind)]
          write_message, str, lun
          str = chisqarr
          write_message, str, lun
          if finite(chisqarr(indbest), /nan) then centroid = 0. else centroid = lambda_smooth(ind(indbest))
          write_message, ['Line centroid = ']+a2str(centroid), lun
        endif
        
        
        ;###########################################
        ;FIT
        
        parinfo = replicate({relstep: 0.5, value:0.D, fixed:0, limited:[0, 0], limits:[0.D, 0]}, nparams)
        ;0:2 => line
        ;3:5 => add locals
        ;6:9 => fringes
        ;9:9+d+1 => poly
        
        ;main compo
        ;parinfo[0:2].relstep = replicate(0.01, 3) ;line parameters
        
        parinfo[0:2].value = [max([maxi, 0.]), centroid, sqrt(sigma^2.+(broadening_sigma(1)*lambda_rest/cspeed)^2.)]
        ;parinfo(0).fixed = 1
        ;parinfo(0).value = 20.
        
        ;scaling factor
        if keyword_set(blend) then begin
          parinfo[0].limited[0] = 1
          parinfo[0].limits[0]  = 0.D
        endif
        ;;test
        ;parinfo[0].fixed = 1
        ;parinfo[0].value  = 2.
        ;wavelength
        parinfo[1].limited = [1, 1]
        parinfo[1].limits  = [centroid-constraints_wavelength*fwhm_um, centroid+constraints_wavelength*fwhm_um]
        ;sigma
        parinfo[2].limited = [1, 1]
        ;parinfo[2].limits  = [(1.-0.01*constraints_sigma(0))*sigma, (1.+0.01*constraints_sigma(1))*sigma]
        parinfo[2].limits  = [sqrt(sigma^2.+(broadening_sigma(0)*lambda_rest/cspeed)^2.), sqrt(sigma^2.+(broadening_sigma(1)*lambda_rest/cspeed)^2.)]
        ;print, sigma, parinfo[2].limits
        
        ;continuum from 1st guess
        ;parinfo[6:6+extra_degree+poly_degree].relstep = replicate(0.01, extra_degree+poly_degree+1)
        ;parinfo[6:6+extra_degree+poly_degree].value = params_cont
        
        ;continuum
        parinfo[indparams_contfringes].value = params_cont ;(with fringes)
        ;set fringes?
        if keyword_set(fringes) then begin
          parinfo[indparams_fringes].value = [1., 1.5, 0.] ;amplitude (Jy), period (um), phase
          parinfo[indparams_fringes(0)].limited = [1, 1]
          parinfo[indparams_fringes(0)].limits  = [0., 2.5]
          parinfo[indparams_fringes(1)].limited = [1, 1]
          parinfo[indparams_fringes(1)].limits  = [0.5, 2.5]
        endif else begin
          parinfo[indparams_fringes].fixed = replicate(1, 3)
          parinfo[indparams_fringes(0)].value = 0.
          parinfo[indparams_fringes(1)].value = 1. ;cannot be zero
          parinfo[indparams_fringes(2)].value = 0.
        endelse
        
        ;other components?
        if keyword_set(addcomponents) then begin
          for ncomp = 0, n_elements(addcomponents) - 1 do begin
          
            addsigma = components_fwhm(ncomp)/(2.*sqrt(2.*alog(2.)))*lambda_rest/cspeed
            ;print, sigma, addlocal, addlocalsigma, sqrt(sigma^2.+addlocalsigma^2.)
            parinfo[3+3*ncomp:3+3*ncomp+2].value = [max([0.5*maxi, 0.]), addcomponents(ncomp)/cspeed*lambda_rest, sqrt(sigma^2.+addsigma^2.)]
            
            if keyword_set(blend) then begin
              parinfo[3+3*ncomp].limited[0] = 1
              parinfo[3+3*ncomp].limits[0]  = 0.D
            endif
            
            parinfo[3+3*ncomp+1].fixed = 1
            parinfo[3+3*ncomp+2].fixed = 1
          ;parinfo[5].limited = [1, 1]
          ;parinfo[5].limits  = [sqrt(sigma^2.-(broadening(0)*lambda_rest/cspeed)^2.), sqrt(sigma^2.+(broadening(1)*lambda_rest/cspeed)^2.)]
          endfor
        endif else begin
          for ncomp = 0, n_elements(addcomponents) -1 do begin
            parinfo[3+3*ncomp:3+3*ncomp+2].value = [0., 0., sigma]
            parinfo[3+3*ncomp:3+3*ncomp+2].fixed = replicate(1, 3)
          endfor
        endelse
        
        ;if not keyword_set(asymetric) then parinfo[3].tied = 'P[2]'
        ;indfit = where( abs(lambda-centroid) lt constraints_continuum(1)*fwhm_um )
        
        noise_eff = noise
        ;if keyword_set(use_errors) then weights = 1./noise else weights = (1./interpol(2.*sqrt(ninbins)*error_smooth, lambda_smooth, lambda)^2.) * max([1.-lambda/(3.*fwhm_um), 0]);1.+fltarr(n_elements(lambda))
        if keyword_set(use_errors) or keyword_set(central3by3combine) or keyword_set(central5by5combine) then $
          weights = 1./noise^2. else weights = 1.*(1./interpol(error_indiv, lambda_smooth, lambda)^2.)
        weights2 = weights * exp(-0.01*lambda^2./fwhm_um^2.) ; max([1.-lambda/(3.*fwhm_um), 0]);1.+fltarr(n_elements(lambda))
        ;err = interpol(error_indiv, lambda_smooth, lambda)
        ;if x eq 2 and y eq 2 then stop
        params = mpfitfun('SpectrumModel', lambda, flux, weights = weights2, parinfo = parinfo, /nan, /quiet, $
          yfit = yfit, STATUS=status, ERRMSG=errmsg, dof = dof, nfree = nfree, bestnorm = bestnorm, perror = perror, $
          functargs={poly_degree: poly_degree, ncomps: n_elements(addcomponents)})
        if status LE 0 then message, errmsg
        
        ;NOW we redo the fit to flag outliers
        ;first calculate sdev
        ;tmp = flux-SpectrumModel(lambda, params, poly_degree = poly_degree)
        ;sdev = robust_sigma(tmp, /zero) ;keep standard deviation
        ;if sdev eq -1 then begin
        ;  print, 'replacing by stdev value'
        ;  sdev = stdev(tmp)
        ;endif
        ;flag and put zero weights to outliers ;not necessarilly necessary :)
        ;ind = where( abs(flux-SpectrumModel(lambda, params, poly_degree = poly_degree)) gt 10.*sdev, count )
        ;if count gt 0 then weights(ind) = 0.
        ;params = mpfitfun('SpectrumModel', lambda, flux, weights = weights, parinfo = parinfo, /nan, /quiet, $
        ;  yfit = yfit, STATUS=status, ERRMSG=errmsg, dof = dof, nfree = nfree, bestnorm = bestnorm, perror = perror, $
        ;  functargs={poly_degree: poly_degree})
        ;if status LE 0 then message, errmsg
        
        ;if keyword_set(use_errors) then weights(noind) = 1./noise(noind) else weights(noind) = (1./interpol(2.*sqrt(ninbins)*error_smooth, lambda_smooth, lambda(noind))^2.) * max([1.-lambda(noind)/(3.*fwhm_um), 0])
        ;if keyword_set(use_errors) then weights(noind) = 1./noise(noind) else weights(noind) = (1./interpol(error_smooth, lambda_smooth, lambda(noind))^2.) * max([1.-lambda(noind)/(3.*fwhm_um), 0])
        if keyword_set(montecarlo) then begin
          parinfo_orig = parinfo
          ;params_final = params
          ;perror_final = perror
          ;bestnorm_final = bestnorm
          tab = dblarr(n_iter_noise, 3)
          tab2 = dblarr(n_iter_noise)
          for ii = 0, n_iter_noise - 1 do begin
            str = 'Monte-Carlo iteration: ' + a2str(ii+1) + '/' + a2str(n_iter_noise)
            write_message, str, lun
            parinfo = parinfo_orig
            new_err = randomn(seed, n_elements(flux), /double) * interpol(error_indiv, lambda_smooth, lambda)
            fluxtmp = flux + new_err
            ;no need to change weights
            params_tmp = mpfitfun('SpectrumModel', lambda, fluxtmp, weights = weights2, parinfo = parinfo, /nan, /quiet, $
              STATUS=status, ERRMSG=errmsg, functargs={poly_degree: poly_degree, ncomps: n_elements(addcomponents)})
            if status LE 0 then message, errmsg
            tab(ii, *) = params_tmp(0:2)
            tab2(ii) = ((params_tmp[0]*params_tmp[2]*sqrt(2.*!PI)) / (((params_tmp[1]+lambda_rest)^2.)*1.e13))*29.98
          endfor
          taberr = dblarr(3)
          for ii = 0, 2 do begin
            ind = where( tab(*, ii) ne 0., c )
            if c eq 0 then ind = where( tab(*, ii) eq tab(*, ii) )
            params(ii) = median(tab(ind, ii))
            taberr(ii) = stdev(tab(ind, ii))
          endfor
          ;print, tab2
          tab2err = stdev(tab2)
          tab2 = median(tab2)
          ;help, tab2, tab2err
          ;params = params_final
          ;perror = perror_final
          ;bestnorm = bestnorm_final
          ind = where( taberr ne 0., c )
          if c gt 0 then perror(ind) = taberr(ind)
        endif
        
        ;recompute params_cont (with fringes)
        params_cont = params(indparams_contfringes)
        ;==========================================================================================
        ;==========================================================================================
        ;==========================================================================================
        
        
        
        
        ;==========================================================================================
        ;==========================================================================================
        ;==========================================================================================
        ;DISPLAY FIT RESULTS
        write_message, 'FIT RESULT: ', lun
        ind = lambda[UNIQ(lambda, SORT(lambda))] ;real dof is for unique x elements
        ;if keyword_set(central3by3combine) or keyword_set(central5by5combine) then dof /= n_elements(lambda)/n_elements(lambda_smooth) ;the combined version replicates values in bins.
        ;dof = n_elements(ind)-nfree
        ;recalculate bestborm (chi^2 using weights and not weights2)
        ;bestnorm = total((flux-SpectrumModel(lambda, params, poly_degree = poly_degree, ncomps = n_elements(addcomponents)))^2. * abs(weights), /nan)
        bestnorm = total((flux[goodinds]-SpectrumModel(lambda[goodinds], params, poly_degree = poly_degree, ncomps = n_elements(addcomponents)))^2. * abs(weights[goodinds]), /nan)
        cube_lres(ind_cube).chisq = bestnorm/dof ;that's the measured -reduced- chi square
        
        ;print, total((flux-SpectrumModel(lambda, params, poly_degree = poly_degree))^2. * abs(weights), /nan)/dof, bestnorm, sqrt(bestnorm/dof)
        ;stop
        
        ;help, dof, n_elements(flux)-nfree, ind-nfree
        ;bestnorm is the same as TOTAL( (flux-SpectrumModel(lambda, params, poly_degree = poly_degree))^2. * ABS(WEIGHTS), /nan )
        
        ;print, bestnorm, dof, nfree, sqrt(bestnorm/dof)
        ;perror is 0 if held fixed or touching a boundary!
        ;if keyword_set(montecarlo) then pcerror = perror else pcerror = perror * cube_lres(ind_cube).chisq ;1 sigma errors assuming dof is correct
        pcerror = perror * sqrt(cube_lres(ind_cube).chisq)
        ;;print results
        ;str = ['Spectral line: signal peak, central wavelength, sigma of gaussian', $
        ;  a2str(params(0:2)) + ' ('+a2str(perror(0:2))+')']
        ;write_message, str, lun
        
        ;process params
        if perror(1) eq 0. then begin ; centroid
          str = 'Forcing error on v to 0.1%. You might wanna loosen this parameter. Final error on flux is probably not correct, consider using MC. '
          write_message, str, lun
          pcerror(1) = 1.e-3 * abs(params(1))
        endif
        if perror(2) eq 0. then begin
          str = 'Forcing error on sigma to same relative uncertainty as peak. You might wanna loosen this parameter. Final error on flux is probably not correct, consider using MC. '
          write_message, str, lun
          pcerror(2) = abs( params(2) * pcerror(0)/abs(params(0)) ) ;* (2.*sqrt(2.*alog(2.))*params(2)/params(0))
        endif
        
        if keyword_set(addcomponents) then begin
          for ncomp = 0, n_elements(addcomponents) - 1 do begin
            if perror(3+3*ncomp+1) eq 0. then begin
              str = 'Component '+a2str(ncomp+2)+': Forcing error on v to 0.1%. You might wanna loosen this parameter. Final error on flux is probably not correct, consider using MC. '
              write_message, str, lun
              pcerror(3+3*ncomp+1) = 1.e-3 * abs(params(3+3*ncomp+1))
            endif
            if perror(3+3*ncomp+2) eq 0. then begin
              str = '2nd component: Forcing error on sigma to same relative uncertainty as peak. You might wanna loosen this parameter. Final error on flux is probably not correct, consider using MC. '
              write_message, str, lun
              pcerror(3+3*ncomp+2) = abs( params(3+3*ncomp+2) * pcerror(3+3*ncomp)/abs(params(3+3*ncomp)) )
            endif
          endfor
        endif
        
        ;print results
        str = ['Spectral line: signal peak, central wavelength, sigma of gaussian', $
          a2str(params(0:2)) + ' ('+a2str(pcerror(0:2))+')']
        write_message, str, lun
        
        if keyword_set(addcomponents) then begin
          for ncomp = 0, n_elements(addcomponents) - 1 do begin
            str = ['Component '+a2str(ncomp+2)+': signal peak, central wavelength, sigma of gaussian', $
              a2str(params(3+3*ncomp:3+3*ncomp+2)) + ' ('+a2str(pcerror(3+3*ncomp:3+3*ncomp+2))+')']
            write_message, str, lun
          endfor
        endif
        
        if keyword_set(fringes) then begin
          str = ['Fringes: amplitude, frequency, phase', $
            a2str(params(indparams_fringes)) + ' ('+a2str(pcerror(indparams_fringes))+')']
          write_message, str, lun
        endif
        
        str = ['Continuum: degree 0, 1, 2, ...', $
          a2str(params(indparams_cont)) + ' ('+a2str(pcerror(indparams_cont))+')']
        write_message, str, lun
        
        str = 'Chi^2: '+a2str(cube_lres(ind_cube).chisq)
        write_message, str, lun
        
        str = [' Continuum [Jy]        = ' + a2str(cube_lres(ind_cube).continuum) + ' +/-  ' + a2str(cube_lres(ind_cube).continuum_error)]
        write_message, str, lun
        ;==========================================================================================
        ;==========================================================================================
        ;==========================================================================================
        
        
        
        
        ;==========================================================================================
        ;==========================================================================================
        ;==========================================================================================
        ;FLUX + ERROR HANDLING
        
        ;fit flux and noise
        ;integral of a Gaussian = sigma * sqrt(2*pi)
        ;area (Jy um) = peak(Jy) * sigma(um) * sqrt(2*pi)
        ;area(W m-2 Hz-1 um) = peak(Jy) * sigma(um) * sqrt(2*pi) * 1e-26
        ;area(W m-2 um-1 um) = peak(Jy) * sigma(um) * sqrt(2*pi) * 1e-26 * 3e14 / lambda(um)^2.
        if keyword_set(montecarlo) then cube_lres(ind_cube).flux = tab2 else cube_lres(ind_cube).flux = ((params[0]*params[2]*sqrt(2.*!PI)) / (((params[1]+lambda_rest)^2.)*1.e13))*29.98
        cube_lres(ind_cube).flux /= (1.+redshift)
        
        flux_rms = ((pcerror[0]*(params[2]+pcerror[2])*sqrt(2.*!PI)) / ((((params[1]-pcerror[1])+lambda_rest)^2.)*1.e13))*29.98 / (1.+redshift)
        
        ;write_message, params[0] , lun
        ;write_message, params[1]+lambda_rest, lun
        ;write_message, params[2], lun
        ;write_message, cube_lres(ind_cube).flux, lun
        cube_lres(ind_cube).wave  = params[1] + lambda_rest
        cube_lres(ind_cube).velocity  = cspeed * params[1] / lambda_rest
        cube_lres(ind_cube).velocity_error = cspeed * pcerror[1] / lambda_rest
        cube_lres(ind_cube).fwhm = cspeed * 2.*sqrt(2.*alog(2.)) * params[2] / lambda_rest
        cube_lres(ind_cube).fwhm_error = cspeed * 2.*sqrt(2.*alog(2.)) * pcerror[2] / lambda_rest
        
        ;errors
        if params(0) eq 0. then pc0 = 0. else pc0 = pcerror(0)/abs(params(0)) ;to avoid NaNs (signal)
        if params(1) eq 0. then pc1 = 0. else pc1 = 2.*pcerror(1)/abs(params(1)+lambda_rest) ;to avoid NaNs (signal)
        if params(2) eq 0. then pc2 = 0. else pc2 = pcerror(2)/params(2) ;to avoid NaNs (sigma)
        ;sigp = 30.
        ;sigs = 10.
        ;p = randomn(seed,5000)*sigp+100
        ;s = randomn(seed,5000)*sigs+50
        ;print, robust_sigma(p*s), median((p*s)*(abs(sigp/p)+abs(sigs/s))), median((p*s)*sqrt((sigp/p)^2.+(sigs/s)^2.))
        ;help, flux_rms, cube_lres(ind_cube).flux,  cube_lres(ind_cube).flux * sqrt( pc0^2. + pc2^2. + pc1^2. )
        ;help, (((params(0)+pcerror[0])*(params[2]+pcerror[2])*sqrt(2.*!PI)) / ((((params[1]-pcerror[1])+lambda_rest)^2.)*1.e13))*29.98 - ((params[0]*params[2]*sqrt(2.*!PI)) / (((params[1]+lambda_rest)^2.)*1.e13))*29.98
        ;if keyword_set(montecarlo) then cube_lres(ind_cube).error = tab2err else cube_lres(ind_cube).error = max([flux_rms, cube_lres(ind_cube).flux]) * sqrt( pc0^2. + pc2^2. + pc1^2. )
        cube_lres(ind_cube).error = abs(cube_lres(ind_cube).flux) * sqrt( pc0^2. + pc2^2. + pc1^2. )
        cube_lres(ind_cube).error /= (1.+redshift)
        
        if cube_lres(ind_cube).flux le 0. and keyword_set(nonegativefluxes) then begin
          write_message, 'Flux is <=0! setting it to 0. ', lun
          cube_lres(ind_cube).flux = 0.
          params(0) = 0.
        endif
        
        ;print, pc0, pc1, pc2
        ;help, params
        ;print, params
        ;print, perror
        ;print, pcerror
        
        if keyword_set(addcomponents) then begin
          for ncomp = 0, n_elements(addcomponents) - 1 do begin
            cube_lres(ind_cube).components_flux(ncomp) = ((params[3+3*ncomp]*params[3+3*ncomp+2]*sqrt(2.*!PI)) / (((params[3+3*ncomp+1]+lambda_rest)^2.)*1.e13))*29.98
            cube_lres(ind_cube).components_flux(ncomp) /= (1.+redshift)
            
            if params(3) eq 0. then pc0 = 0. else pc0 = pcerror(3+3*ncomp)/params(3+3*ncomp) ;to avoid NaNs (signal)
            if params(4) eq 0. then pc1 = 0. else pc1 = 2.*pcerror(3+3*ncomp+1)/(params(3+3*ncomp+1)+lambda_rest) ;to avoid NaNs (sigma)
            if params(5) eq 0. then pc2 = 0. else pc2 = pcerror(3+3*ncomp+2)/params(3+3*ncomp+2) ;to avoid NaNs (sigma)
            cube_lres(ind_cube).components_error(ncomp) = abs(cube_lres(ind_cube).components_flux(ncomp)) * sqrt( pc0^2. + pc2^2. + pc1^2. )
            cube_lres(ind_cube).components_error(ncomp) /= (1.+redshift)
            
            if cube_lres(ind_cube).components_flux(ncomp) le 0. then begin
              write_message, '(Extra component'+a2str(ncomp+1)+') flux is <=0! setting it to 0. ', lun
              cube_lres(ind_cube).components_flux(ncomp) = 0.
              params(3+3*ncomp) = 0.
            endif
          endfor
        endif
        
        ;help, params, cube_fitparams(ind_cube).params(0:9+poly_degree)
        cube_fitparams(ind_cube).params(0:nparams-1) = params
        cube_fitparams(ind_cube).params_errors(0:nparams-1) = pcerror
        
        ;integrated flux
        ran = where( abs(lambda_smooth-params[1]) lt constraints_continuum(0)*fwhm_um, c )
        if c eq 0 then begin
          str = 'It looks like the wavelength range has some problems... no points for line. Min and max of lambda corrected for redshift: '
          write_message, str, lun
          str = minmax(lambda_cube, /nan)
          write_message, str, lun
          close, /all
          retall
        endif
        ;step = median(lambda_smooth-shift(lambda_smooth, 1))
        contmodel = ContinuumModel(lambda_smooth(ran), params_cont, poly_degree = poly_degree)
        ;calculate equivalent width
        ;eqw = 0.
        ;for l = 0, n_elements(ran) - 1 do eqw += intfluxdensity(l) / contmodel(l) * bin
        ;write_message, ['EQW_rest [um]          = '+a2str(eqw)], lun
        cube_lres(ind_cube).eqw = total(SpectrumModel(lambda_smooth(ran), params, poly_degree = poly_degree, ncomps = n_elements(addcomponents))/ContinuumModel(lambda_smooth(ran), params_cont, poly_degree = poly_degree) - 1., /nan) * bin
        
        intfluxdensity = flux_smooth(ran) - contmodel
        ;intflux = total((intfluxdensity*3.e-12/(params[1]+lambda_rest)^2.)*step, /nan)
        cube_lres(ind_cube).integrated_flux = tsum(lambda_smooth(ran), (intfluxdensity*3.e-12/(params[1]+lambda_rest)^2.))
        cube_lres(ind_cube).integrated_flux /= (1.+redshift)
        
        ran = where( abs(lambda_smooth-params[1]) lt fwhm_um, c )
        if c gt 1 then cube_lres(ind_cube).integrated_noise = tsum(lambda_smooth(ran), error_smooth(ran)*(3.e-12/(params[1]+lambda_rest)^2.) )
        cube_lres(ind_cube).integrated_noise /= (1.+redshift)
        ;intnoise = dblarr(n_iter_noise)
        ;for r = 0, n_iter_noise - 1 do begin
        ;  intfluxdensity = dblarr(n_elements(ran))
        ;  randomNumbers = !RNG -> GetRandomNumbers(n_elements(ran))
        ;  intfluxdensity = flux_smooth(ran) + error_smooth(ran)*4.*(randomNumbers-0.5) - ContinuumModel(lambda_smooth(ran), params_cont, poly_degree = poly_degree)
        ;  ;for ii = 0, n_elements(ran) - 1 do begin
        ;  ;  intfluxdensity(ii) = flux_smooth(ran(ii)) + error_smooth(ran(ii))*2.*(randomu(seed)-0.5) - PolyModel(lambda_smooth(ran(ii)), params_poly, poly_degree = poly_degree)
        ;  ;endfor
        ;  arr = (intfluxdensity*3.e-12/(params[1]+lambda_rest)^2.)
        ;  ind = where( finite(arr, /nan), c )
        ;  if c gt 0 then arr[ind] = 0.
        ;  intnoise(r) = tsum(lambda_smooth(ran), arr)
        ;;print, intflux, intnoise(r)
        ;endfor
        ;ind = where(finite(intnoise))
        ;intnoise = robust_sigma(intnoise(ind))
        
        ;if fitnoise eq 0. then begin
        ;  str = ['Error on flux is zero (happens when touching the boundary of one parameter)', 'Replacing error with error on the integrated flux.']
        ;  write_message, str, lun
        ;  fitnoise = intnoise
        ;endif
        ;m = moment(intnoise(ind), sdev = sdev)
        ;intnoise = sdev
        
        ;replace flux by noise if less (can be zero)
        ;if fitflux eq 0. then begin
        ;  str = ['Flux is zero (happens when touching the boundary of one parameter)']
        ;  write_message, str, lun
        ;  ;fitflux = fitnoise
        ;endif
        
        ;if keyword_set(integrate) then begin
        ;  cube_lres(x, y, n-1, 0) = intflux
        ;  cube_lres(x, y, n-1, 1) = intnoise
        ;endif else begin
        
        write_message, 'FLUX (GAUSSIAN FIT): ', lun
        str = [' Flux [W m-2]          = '+a2str(cube_lres(ind_cube).flux) + ' +/- ' +a2str(cube_lres(ind_cube).error)]
        write_message, str, lun
        ;if keyword_set(addlocal) then begin
        ;  str = ['Flux (1+2) [W m-2]    = '+a2str(cube_lres(ind_cube).flux+cube_lres(ind_cube).localcomponent_flux) + ' +/- ' +a2str(cube_lres(ind_cube).error+cube_lres(ind_cube).localcomponent_error)]
        ;  write_message, str, lun
        ;endif
        
        write_message, 'FLUX (INTEGRATED): ', lun
        str = [' Flux [W m-2]          = '+a2str(cube_lres(ind_cube).integrated_flux) + ' +/- ' +a2str(cube_lres(ind_cube).integrated_noise)]
        write_message, str, lun
        if finite(cube_lres(ind_cube).integrated_flux, /nan) then begin
          str = 'Integrated flux is a NaN, maybe the sampling is too high?'
          write_message, str, lun
        endif
        
        if keyword_set(addcomponents) then begin
          for ncomp = 0, n_elements(addcomponents) - 1 do begin
            write_message, 'EXTRA COMPONENT '+a2str(ncomp+1)+' FLUX (GAUSSIAN FIT): ', lun
            str = [' Flux [W m-2]          = '+a2str(cube_lres(ind_cube).components_flux(ncomp)) + ' +/- ' +a2str(cube_lres(ind_cube).components_error(ncomp))]
            write_message, str, lun
          endfor
        endif
        ;==========================================================================================
        ;==========================================================================================
        ;==========================================================================================
        
        
        
        
        
        ;==========================================================================================
        ;==========================================================================================
        ;==========================================================================================
        ;PLOTS
        if keyword_set(dontdo_plots_publi) then piarr = [0, 0] else piarr = [0, 1]
        for pi = piarr[0], piarr[1] do begin
        
          case pi of
            0: begin
              plotdir = dir+'/PLOTS_' + output_name
              plots_publi = 0
            end
            1: begin
              plotdir = dir+'/PLOTSpubli_' + output_name
              plots_publi = 1
            end
          endcase
          
          
          ;###########################################
          ;###########################################
          ;###########################################
          ;###########################################
          ;full plot
          
          if keyword_set(display) then set_plot, 'X' else begin
            ;set_plot, 'z'
            ;ERASE
            ;device, set_resolution=reso, set_pixel_depth = 24, decomposed = 0
            set_plot, 'PS'
            !p.font = 0
            device, filename=plotdir + '/raw/raster'+a2str(n)+'_'+a2str(x+1)+a2str(y+1)+'_raw.eps', /encaps, /color, /helvetica
          ;device,
          ;filename=dir+'/PS/cube'+a2str(n)+'_'+a2str(x)+a2str(y)+'.ps',
          ;/color
          endelse
          loadct, 0, /silent
          xtit = textoidl('\lambda_{rest}-\lambda_{0} [\mum]')
          xran = [max([0.-viewrange(0)*fwhm_um, 1.1*min(lambda_smooth, /nan)], /nan), min([viewrange(1)*fwhm_um, 1.1*max(lambda_smooth, /nan)], /nan)]
          ;ind = where( lambda_flagged gt xran(0) and lambda_flagged lt xran(1), c )
          ;if c gt 0 then yran = minmax(flux_flagged(ind), /nan) else yran = minmax(flux_smooth)
          ind = where( lambda gt xran(0) and lambda lt xran(1) )
          yran = [min([flux(ind), SpectrumModel(lambda(indfit_smooth), params, poly_degree = poly_degree, ncomps = n_elements(addcomponents))], /nan), max([flux(ind), SpectrumModel(lambda(indfit_smooth), params, poly_degree = poly_degree, ncomps = n_elements(addcomponents))], /nan)]
          plot, indgen(10), /nodata, yran=yran, xran=xran, ysty = 1, $
            xtitle = xtit, ytitle = 'Flux density [Jy]', background = 255, color = 0, xstyle = 1+8, ymargin = [5, 5];, xmargin = [11, 11]
          axis, xaxis = 1, xrange = velocity, color = 0, xstyle = 1, xtitle = 'Velocity [km s'+textoidl('^{-1}')+']'
          ;axis, yaxis = 1, ysty = 1, yran = yran * 1.e-6 / (9.4*!pi/180./3600.)^2., color = 0, ytitle = 'Flux density [MJy sr'+textoidl('^{-1}')+']'
          
          xt = 0.02
          yt = 0.95
          if plots_publi eq 1 then begin
            st = line_label(line_index) + ', RA='+a2str(0.0001*nint(10000.*cube_lres(ind_cube).ra))+', DEC='+a2str(0.0001*nint(10000.*cube_lres(ind_cube).dec))
            cs = 0.75
          endif else begin
            st = line_label(line_index) + ', Raster ' + a2str(n)+', spaxel ('+a2str(x+1)+','+a2str(y+1)+')'
            cs = 0.75
          endelse
          xyouts, xt, yt, st, /normal, alignment=0, charsize=cs, color = 0
          velocity = cspeed * xran / lambda_rest
          
          plots, [0., 0.], !y.crange, linestyle = 0, color = 180, noclip = 0
          if plots_publi eq 0 then begin
            plots, params[1]*[1.,1.], !y.crange, linestyle = 2, color = 180
            plots, (0.-constraints_continuum(0)*fwhm_um)*[1.,1.], !y.crange, linestyle = 0, color = 150, noclip = 0
            plots, (0.+constraints_continuum(0)*fwhm_um)*[1.,1.], !y.crange, linestyle = 0, color = 150, noclip = 0
            plots, (0.-constraints_continuum(1)*fwhm_um)*[1.,1.], !y.crange, linestyle = 0, color = 150, noclip = 0
            plots, (0.+constraints_continuum(2)*fwhm_um)*[1.,1.], !y.crange, linestyle = 0, color = 150, noclip = 0
          ;plots, [0.-params[1], 0.-params[1]], !y.crange, linestyle = 0, color = 200, noclip = 0
          endif
          
          nsigma = cube_lres(ind_cube).flux/cube_lres(ind_cube).error
          nsigma = max([nsigma, 0])
          if finite(nsigma, /nan) then nsigma = 0.
          
          ;if fitflux le fitnoise then fitstr = '<'+str_round(a2str(fitnoise))+' W m-2' else $
          fitstr = str_round(a2str(cube_lres(ind_cube).flux)) + '+/-'
          tmp = cube_lres(ind_cube).error
          ;if not finite(cube_lres(ind_cube).localcomponent_error, /nan) then tmp += cube_lres(ind_cube).localcomponent_error
          fitstr += str_round(a2str(tmp)) + ' W m'+textoidl('^{-2}')
          ;fitstr = str_round(a2str(cube_lres(ind_cube).flux+cube_lres(ind_cube).localcomponent_flux))+'+/-'+str_round(a2str(cube_lres(ind_cube).error+cube_lres(ind_cube).localcomponent_error))+
          ;if intflux lt intnoise then intstr = 'F_int < '+str_round(a2str(intnoise))+' W m-2' else $
          ;  intstr = 'F_int = '+str_round(a2str(intflux))+'+/-'+str_round(a2str(intnoise))+' W m-2'
          
          cube_lres(ind_cube).broadening = max([sqrt(cube_lres(ind_cube).fwhm^2.-fwhm_kms^2.), 0.], /nan)
          cube_lres(ind_cube).broadening_error = cube_lres(ind_cube).broadening * cube_lres(ind_cube).fwhm_error / cube_lres(ind_cube).fwhm
          
          if plots_publi eq 0 then al_legend, /top, /right, charsize = 0.7, [fitstr, str_round(a2str(nsigma))+textoidl('\sigma'), $
            'EQW'+textoidl('_{rest}')+' = '+str_round(a2str(cube_lres(ind_cube).eqw))+ ' '+textoidl('\mu')+'m', $
            'v = '+str_round(a2str(cube_lres(ind_cube).velocity))+'+/-'+str_round(a2str(cube_lres(ind_cube).velocity_error))+' km s'+textoidl('^{-1}'), $
            'Broadening = '+str_round(a2str(cube_lres(ind_cube).broadening))+'+/-'+str_round(a2str(cube_lres(ind_cube).broadening_error))+ ' km s'+textoidl('^{-1}')], $
            textcolor = [0, 0, 0, 0, 0], box = 0
            
          if flagged(0) ne -1 then oplot, lambda_raw(flagged), flux_raw(flagged), color=150, psym=3, symsize=0.1
          if plots_publi eq 0 then al_legend, /bottom, /left, ['flagged'], textcolor=[150], box = 0
          
          loadct, 39, /silent
          if flagged_sigmaclip(0) ne - 1 then begin
            oplot, lambda_flagged(flagged_sigmaclip), flux_flagged(flagged_sigmaclip), color=120, psym=3, symsize=0.2
            if plots_publi eq 0 then al_legend, /bottom, /left, ['         (>'+a2str(nint(n_sigma_clipping))+''+textoidl('\sigma')+')'], textcolor=[120], box = 0
          endif
          
          ;threshold = 2.e8     ;if number of points > threshold then use frebin for a cleaner plot
          ;if readouts gt threshold then begin
          ;  print, 'rebinning to 2e8 elements for a cleaner plot'
          ;  oplot, frebin(lambda, threshold), frebin(flux, threshold), psym=3, color=50, symsize=0.1
          ;  oplot, frebin(lambda(continuum), threshold), frebin(flux(continuum), threshold), color=80, psym=3, symsize=0.1
          ;endif else begin
          oplot, lambda_fullrange, flux_fullrange, psym=3, color=50, symsize=0.1
          oplot, lambda(continuum_ind), flux(continuum_ind), color=80, psym=3, symsize=0.1
          ;endelse
          
          oploterror, lambda_smooth, flux_smooth, error_indiv, color=200, psym=-3, symsize=0.1, errcolor = 200 ;error_smooth*sqrt(ninbins)
          oploterror, lambda_smooth, flux_smooth, error_smooth, color=200, psym=-3, symsize=0.1, errcolor = 200 ;error_smooth*sqrt(ninbins)
          ;oplot, lambda_smooth, flux_smooth+error_smooth, color=50, psym=-3, symsize=0.1;, linestyle = 2
          ;oplot, lambda_smooth, flux_smooth-error_smooth, color=50, psym=-3, symsize=0.1;, linestyle = 2
          oplot, lambda(indfit_smooth), SpectrumModel(lambda(indfit_smooth), params, poly_degree = poly_degree, ncomps = n_elements(addcomponents)), color=250, thick=3
          oplot, lambda(indfit_smooth), LineModel(lambda(indfit_smooth), params[0:2])+ContinuumModel(lambda(indfit_smooth), params_cont, poly_degree = poly_degree), color=250, thick=1, linestyle = 0
          if keyword_set(addcomponents) then begin
            for ncomp = 0, n_elements(addcomponents) - 1 do begin
              oplot, lambda(indfit_smooth), LineModel(lambda(indfit_smooth), params[3+3*ncomp:3+3*ncomp+2])+ContinuumModel(lambda(indfit_smooth), params_cont, poly_degree = poly_degree), color=250, thick=1, linestyle = 0
            endfor
          endif
          oplot, lambda(indfit_smooth), ContinuumModel(lambda(indfit_smooth), params_cont, poly_degree = poly_degree), color=250, thick=1, linestyle = 0
          
          if plots_publi eq 0 then al_legend, /top, /left, ['data', 'data (continuum)', 'rebinned data', 'fit'], textcolor=[50, 80, 200, 250], box = 0
          
          loadct, 0, /silent
          if not keyword_set(display) then begin
            ;  write_png, dir+'/PLOTS_' + output_name + '/raw/raster'+a2str(n)+'_'+a2str(x+1)+a2str(y+1)+'_raw.png', tvrd(/true)
            device, /close
            set_plot, 'X'
          endif
          
          ;###########################################
          ;###########################################
          ;###########################################
          ;###########################################
          ;rebinned plot
          if keyword_set(display) then set_plot, 'X' else begin
            ;set_plot, 'z'
            ;ERASE
            ;device, set_resolution=reso, set_pixel_depth = 24, decomposed = 0
            set_plot, 'PS'
            !p.font = 0
            device, filename=plotdir + '/raster'+a2str(n)+'_'+a2str(x+1)+a2str(y+1)+'.eps', /color, /encaps, /helvetica
          endelse
          loadct, 0, /silent
          xtit = textoidl('\lambda_{rest}-\lambda_{0} [\mum]')
          xran = [max([0.-viewrange(0)*fwhm_um, 1.1*min(lambda_smooth, /nan)], /nan), min([viewrange(1)*fwhm_um, 1.1*max(lambda_smooth, /nan)], /nan)]
          ind = where( lambda_smooth gt xran(0) and lambda_smooth lt xran(1) )
          yran = [min([flux_smooth(ind), SpectrumModel(lambda(indfit_smooth), params, poly_degree = poly_degree, ncomps = n_elements(addcomponents))], /nan), max([flux_smooth(ind), SpectrumModel(lambda(indfit_smooth), params, poly_degree = poly_degree, ncomps = n_elements(addcomponents))], /nan)]
          
          if plots_publi eq 1 then cs = 1. else cs = 1.
          
          plot, indgen(10), /nodata, yran=yran, xran=xran, charsize = cs, $
            xtitle = xtit, ytitle = 'Flux density [Jy]', background = 255, color = 0, xstyle = 1+8, ymargin = [5, 5], ystyle = 1;, xmargin = [11, 11]
          ;axis, yaxis = 1, ysty = 1, yran = yran * 1.e-6 / (9.4*!pi/180./3600.)^2., charsize = cs, color = 0, ytitle = 'Flux density [MJy sr'+textoidl('^{-1}')+']'
          axis, xaxis = 1, xrange = velocity, charsize = cs, color = 0, xstyle = 1, xtitle = 'Velocity [km s'+textoidl('^{-1}')+']'
          
          xt = 0.02
          yt = 0.95
          if plots_publi eq 1 then begin
            st = 'RA='+a2str(0.0001*nint(10000.*cube_lres(ind_cube).ra))+', DEC='+a2str(0.0001*nint(10000.*cube_lres(ind_cube).dec))
            cs = 0.75
            coldata = 0
            colerr = 100
            colfit = 250
          endif else begin
            st = line_label(line_index) + ', Raster ' + a2str(n)+', spaxel ('+a2str(x+1)+','+a2str(y+1)+')'
            cs = 0.75
            coldata = 100
            colerr = 200
            colfit = 250
          endelse
          xyouts, xt, yt, st, /normal, alignment=0., charsize=cs, color = 0
          velocity = cspeed * xran / lambda_rest
          
          plots, [0., 0.], !y.crange, linestyle = 0, color = 180, noclip = 0
          if plots_publi eq 0 then begin
            plots, [params[1], params[1]], !y.crange, linestyle = 2, color = 180, noclip = 0
            plots, (0.-constraints_continuum(0)*fwhm_um)*[1.,1.], !y.crange, linestyle = 0, color = 150, noclip = 0
            plots, (0.+constraints_continuum(0)*fwhm_um)*[1.,1.], !y.crange, linestyle = 0, color = 150, noclip = 0
            plots, (0.-constraints_continuum(1)*fwhm_um)*[1.,1.], !y.crange, linestyle = 0, color = 150, noclip = 0
            plots, (0.+constraints_continuum(2)*fwhm_um)*[1.,1.], !y.crange, linestyle = 0, color = 150, noclip = 0
          ;plots, [0.-params[1], 0.-params[1]], !y.crange, linestyle = 0, color = 200, noclip = 0
          endif
          
          nsigma = cube_lres(ind_cube).flux/cube_lres(ind_cube).error
          nsigma = max([nsigma, 0])
          if finite(nsigma, /nan) then nsigma = 0.
          
          ;if fitflux le fitnoise then fitstr = '<'+str_round(a2str(fitnoise))+' W m-2' else $
          
          ;fitstr = str_round(a2str(cube_lres(ind_cube).flux+cube_lres(ind_cube).localcomponent_flux))+'+/-'+str_round(a2str(cube_lres(ind_cube).error+cube_lres(ind_cube).localcomponent_error))+' W m'+textoidl('^{-2}')        ;if intflux lt intnoise then intstr = 'F_int < '+str_round(a2str(intnoise))+' W m-2' else $
          ;  intstr = 'F_int = '+str_round(a2str(intflux))+'+/-'+str_round(a2str(intnoise))+' W m-2'
          
          if plots_publi eq 1 then al_legend, /top, /left, box = 1, line_label(line_index), textcolor = 0 else al_legend, /top, /right, charsize = 0.7, [fitstr, str_round(a2str(nsigma))+textoidl('\sigma'), $
            'EQW'+textoidl('_{rest}')+' = '+str_round(a2str(cube_lres(ind_cube).eqw))+ ' '+textoidl('\mu')+'m', $
            'v = '+str_round(a2str(cube_lres(ind_cube).velocity))+'+/-'+str_round(a2str(cube_lres(ind_cube).velocity_error))+' km s'+textoidl('^{-1}'), $
            'Broadening = '+str_round(a2str(cube_lres(ind_cube).broadening))+'+/-'+str_round(a2str(cube_lres(ind_cube).broadening_error))+ ' km s'+textoidl('^{-1}')], $
            textcolor = [0, 0, 0, 0, 0], box = 0
            
          loadct, 0, /silent
          
          oploterror, lambda_smooth, flux_smooth, error_smooth, color=coldata, psym=10, symsize=0.1, errcolor = colerr
          oplot, lambda_smooth, flux_smooth, color=coldata, psym=10, symsize=0.1, thick = 2
          
          loadct, 39, /silent
          oplot, lambda(indfit_smooth), SpectrumModel(lambda(indfit_smooth), params, poly_degree = poly_degree, ncomps = n_elements(addcomponents)), color=colfit, thick=3, linestyle = 0
          oplot, lambda(indfit_smooth), LineModel(lambda(indfit_smooth), params[0:2])+ContinuumModel(lambda(indfit_smooth), params_cont, poly_degree = poly_degree), color=colfit, thick=1, linestyle = 0
          if keyword_set(addcomponents) then begin
            for ncomp = 0, n_elements(addcomponents) - 1 do begin
              oplot, lambda(indfit_smooth), LineModel(lambda(indfit_smooth), params[3+3*ncomp:3+3*ncomp+2])+ContinuumModel(lambda(indfit_smooth), params_cont, poly_degree = poly_degree), color=250, thick=1, linestyle = 0
            endfor
          endif
          
          ;plot with theoretical spectral resolution
          params_tmp = params
          params_tmp(2) = sigma
          ;params_tmp(5) = sigma
          oplot, lambda(indfit_smooth), SpectrumModel(lambda(indfit_smooth), params_tmp, poly_degree = poly_degree, ncomps = n_elements(addcomponents)), color=colfit, thick = 1, linestyle = 2
          
          ;plot 3-sigma upper limit
          if nsigma lt 1. then begin
            params_tmp = params
            params_tmp(0) = 1.*pcerror(0)
            oplot, lambda(indfit_smooth), SpectrumModel(lambda(indfit_smooth), params_tmp, poly_degree = poly_degree, ncomps = n_elements(addcomponents)), color=colfit, thick = 1, linestyle = 1
            params_tmp(0) = 4.*pcerror(0) ;that's 2 sigma, and to keep the line profile close to expected, I put 2x in the peak and 2x for sigma in the peak as well
            oplot, lambda(indfit_smooth), SpectrumModel(lambda(indfit_smooth), params_tmp, poly_degree = poly_degree, ncomps = n_elements(addcomponents)), color=colfit, thick = 1, linestyle = 1
            params_tmp(0) = 6.*pcerror(0)
            oplot, lambda(indfit_smooth), SpectrumModel(lambda(indfit_smooth), params_tmp, poly_degree = poly_degree, ncomps = n_elements(addcomponents)), color=colfit, thick = 1, linestyle = 1
          ;params_tmp(0) = 20.*pcerror(0)
          ;oplot, lambda(indfit_smooth), SpectrumModel(lambda(indfit_smooth), params_tmp, poly_degree = poly_degree), color=colfit, thick = 1, linestyle = 1
          endif
          
          oplot, lambda(indfit_smooth), ContinuumModel(lambda(indfit_smooth), params_cont, poly_degree = poly_degree), color=colfit, thick=1, linestyle = 0
          
          if plots_publi eq 0 then al_legend, /top, /left, charsize = 0.7, ['FWHM='+str_round(a2str(cube_lres(ind_cube).fwhm))+'+/-'+str_round(a2str(cube_lres(ind_cube).fwhm_error))+' km s'+textoidl('^{-1}'), 'FWHM_th='+str_round(a2str(fwhm_kms))+' km s'+textoidl('^{-1}')], textcolor=[250, 250], box = 0, linestyle=[0, 2], color = [250, 250], thick = [3, 1]
          
          loadct, 0, /silent
          if not keyword_set(display) then begin
            ;write_png, dir+'/PLOTS_' + output_name + '/raster'+a2str(n)+'_'+a2str(x+1)+a2str(y+1)+'.png', tvrd(/true)
            device, /close
            set_plot, 'X'
          endif
          
          ;###########################################
          ;###########################################
          ;###########################################
          ;###########################################
          ;normalized plot
          if keyword_set(display) then set_plot, 'X' else begin
            ;set_plot, 'z'
            ;ERASE
            ;device, set_resolution=reso, set_pixel_depth = 24, decomposed = 0
            set_plot, 'PS'
            !p.font = 0
            device, filename=plotdir + '/continuum_normalized/raster'+a2str(n)+'_'+a2str(x+1)+a2str(y+1)+'_normalized.eps', /encaps, /color, /helvetica, xsize = 15, ysize = 10
          endelse
          
          cube_spectra(ind_cube).normalized_flux(0:n_elements(lambda_smooth)-1) = flux_smooth / ContinuumModel(lambda_smooth, params_cont, poly_degree = poly_degree)
          
          loadct, 0, /silent
          xtit = textoidl('\lambda_{rest}-\lambda_{0} [\mum]')
          xran = [max([0.-viewrange(0)*fwhm_um, 1.1*min(lambda_smooth, /nan)], /nan), min([viewrange(1)*fwhm_um, 1.1*max(lambda_smooth, /nan)], /nan)]
          ind = where( lambda_smooth gt xran(0) and lambda_smooth lt xran(1) )
          plot, indgen(10), /nodata, yran=[-0.2, $
            max([flux_smooth(ind), SpectrumModel(lambda_smooth, params, poly_degree = poly_degree, ncomps = n_elements(addcomponents))]/ContinuumModel(lambda_smooth, params_cont, poly_degree = poly_degree), /nan)], xran=xran, $
            xtitle = xtit, ytitle = 'Flux density (continuum-normalized)', background = 255, color = 0, xstyle = 1+8, ymargin = [5, 5];, xmargin = [7, 7]
            
          xt = 0.02
          yt = 0.95
          if plots_publi eq 1 then begin
            st = ''; line_label(line_index) + ', RA='+a2str(0.0001*nint(10000.*cube_lres(ind_cube).ra))+', DEC='+a2str(0.0001*nint(10000.*cube_lres(ind_cube).dec))
            cs = 0.75
            coldata = 0
            colerr = 100
            colfit = 250
          endif else begin
            st = line_label(line_index) + ', Raster ' + a2str(n)+', spaxel ('+a2str(x+1)+','+a2str(y+1)+')'
            cs = 0.75
            coldata = 100
            colerr = 200
            colfit = 250
          endelse
          xyouts, xt, yt, st, /normal, alignment=0., charsize=1.25, color = 0
          velocity = cspeed * xran / lambda_rest
          
          plots, [0., 0.], !y.crange, linestyle = 0, color = 180, noclip = 0
          if plots_publi eq 0 then begin
            plots, [params[1], params[1]], !y.crange, linestyle = 2, color = 180, noclip = 0
            plots, (0.-constraints_continuum(0)*fwhm_um)*[1.,1.], !y.crange, linestyle = 0, color = 150, noclip = 0
            plots, (0.+constraints_continuum(0)*fwhm_um)*[1.,1.], !y.crange, linestyle = 0, color = 150, noclip = 0
            plots, (0.-constraints_continuum(1)*fwhm_um)*[1.,1.], !y.crange, linestyle = 0, color = 150, noclip = 0
            plots, (0.+constraints_continuum(2)*fwhm_um)*[1.,1.], !y.crange, linestyle = 0, color = 150, noclip = 0
          ;plots, [0.-params[1], 0.-params[1]], !y.crange, linestyle = 0, color = 200, noclip = 0
          endif
          
          nsigma = cube_lres(ind_cube).flux/cube_lres(ind_cube).error
          nsigma = max([nsigma, 0])
          if finite(nsigma, /nan) then nsigma = 0.
          
          ;if fitflux le fitnoise then fitstr = '<'+str_round(a2str(fitnoise))+' W m-2' else $
          
          ;fitstr = str_round(a2str(cube_lres(ind_cube).flux+cube_lres(ind_cube).localcomponent_flux))+'+/-'+str_round(a2str(cube_lres(ind_cube).error+cube_lres(ind_cube).localcomponent_error))+' W m'+textoidl('^{-2}')        ;if intflux lt intnoise then intstr = 'F_int < '+str_round(a2str(intnoise))+' W m-2' else $
          ;  intstr = 'F_int = '+str_round(a2str(intflux))+'+/-'+str_round(a2str(intnoise))+' W m-2'
          if plots_publi eq 1 then al_legend, /top, /left, box = 1, line_label(line_index), textcolor = 0 else al_legend, /top, /right, charsize = 0.7, [fitstr, str_round(a2str(nsigma))+textoidl('\sigma'), $
            'EQW'+textoidl('_{rest}')+' = '+str_round(a2str(cube_lres(ind_cube).eqw))+ ' '+textoidl('\mu')+'m', $
            'v = '+str_round(a2str(cube_lres(ind_cube).velocity))+'+/-'+str_round(a2str(cube_lres(ind_cube).velocity_error))+' km s'+textoidl('^{-1}'), $
            'Broadening = '+str_round(a2str(cube_lres(ind_cube).broadening))+'+/-'+str_round(a2str(cube_lres(ind_cube).broadening_error))+ ' km s'+textoidl('^{-1}')], $
            textcolor = [0, 0, 0, 0, 0], box = 0
            
          loadct, 0, /silent
          
          oplot, lambda_smooth, replicate(1., n_elements(lambda_smooth)), color=0, thick=1, linestyle = 0
          oplot, lambda_smooth, replicate(0., n_elements(lambda_smooth)), color=0, thick=2, linestyle = 0
          
          oploterror, lambda_smooth, flux_smooth/ContinuumModel(lambda_smooth, params_cont, poly_degree = poly_degree), error_smooth/ContinuumModel(lambda_smooth, params_cont, poly_degree = poly_degree), color=coldata, psym=10, symsize=0.1, errcolor = colerr
          oplot, lambda_smooth, flux_smooth/ContinuumModel(lambda_smooth, params_cont, poly_degree = poly_degree), color=coldata, psym=10, symsize=0.1, thick = 2
          ;print, flux_smooth/ContinuumModel(lambda_smooth, params_cont, poly_degree = poly_degree)
          
          loadct, 39, /silent
          oplot, lambda(indfit_smooth), SpectrumModel(lambda(indfit_smooth), params, poly_degree = poly_degree, ncomps = n_elements(addcomponents))/ContinuumModel(lambda(indfit_smooth), params_cont, poly_degree = poly_degree), color=colfit, thick=3, linestyle = 0
          oplot, lambda(indfit_smooth), (LineModel(lambda(indfit_smooth), params[0:2])+ContinuumModel(lambda(indfit_smooth), params_cont, poly_degree = poly_degree))/ContinuumModel(lambda(indfit_smooth), params_cont, poly_degree = poly_degree), color=colfit, thick=1, linestyle = 0
          if keyword_set(addcomponents) then begin
            for ncomp = 0, n_elements(addcomponents) - 1 do begin
              oplot, lambda(indfit_smooth), (LineModel(lambda(indfit_smooth), params[3+3*ncomp:3+3*ncomp+2])+ContinuumModel(lambda(indfit_smooth), params_cont, poly_degree = poly_degree))/ContinuumModel(lambda(indfit_smooth), params_cont, poly_degree = poly_degree), color=colfit, thick=1, linestyle = 0
            endfor
          endif
          
          ;plot with theoretical resolution
          params_tmp = params
          params_tmp(2) = sigma
          ;params_tmp(5) = sigma
          oplot, lambda(indfit_smooth), SpectrumModel(lambda(indfit_smooth), params_tmp, poly_degree = poly_degree, ncomps = n_elements(addcomponents))/ContinuumModel(lambda(indfit_smooth), params_cont, poly_degree = poly_degree), color=colfit, thick = 1, linestyle = 2
          
          ;plot 3-sigma upper limit
          if nsigma lt 1. then begin
            params_tmp = params
            params_tmp(0) = 1.*pcerror(0)
            oplot, lambda(indfit_smooth), SpectrumModel(lambda(indfit_smooth), params_tmp, poly_degree = poly_degree, ncomps = n_elements(addcomponents))/ContinuumModel(lambda(indfit_smooth), params_cont, poly_degree = poly_degree), color=colfit, thick = 1, linestyle = 1
            params_tmp(0) = 4.*pcerror(0) ;that's 2 sigma, and to keep the line profile close to expected, I put 2x in the peak and 2x for sigma in the peak as well
            oplot, lambda(indfit_smooth), SpectrumModel(lambda(indfit_smooth), params_tmp, poly_degree = poly_degree, ncomps = n_elements(addcomponents))/ContinuumModel(lambda(indfit_smooth), params_cont, poly_degree = poly_degree), color=colfit, thick = 1, linestyle = 1
            params_tmp(0) = 6.*pcerror(0)
            oplot, lambda(indfit_smooth), SpectrumModel(lambda(indfit_smooth), params_tmp, poly_degree = poly_degree, ncomps = n_elements(addcomponents))/ContinuumModel(lambda(indfit_smooth), params_cont, poly_degree = poly_degree), color=colfit, thick = 1, linestyle = 1
          ;params_tmp(0) = 20.*pcerror(0)
          ;oplot, lambda(indfit_smooth), SpectrumModel(lambda(indfit_smooth), params_tmp, poly_degree = poly_degree), color=colfit, thick = 1, linestyle = 1
          endif
          
          oplot, lambda(indfit_smooth), ContinuumModel(lambda(indfit_smooth), params_cont, poly_degree = poly_degree)/ContinuumModel(lambda(indfit_smooth), params_cont, poly_degree = poly_degree), color=colfit, thick=1, linestyle = 0
          
          if plots_publi eq 0 then al_legend, /top, /left, charsize = 0.7, ['FWHM='+str_round(a2str(cube_lres(ind_cube).fwhm))+'+/-'+str_round(a2str(cube_lres(ind_cube).fwhm_error))+' km s'+textoidl('^{-1}'), 'FWHM_th='+str_round(a2str(fwhm_kms))+' km s'+textoidl('^{-1}')], textcolor=[250, 250], box = 0, linestyle=[0, 2], color = [250, 250], thick = [3, 1]
          
          loadct, 0, /silent
          axis, xaxis=1, xrange = velocity, color = 0, xstyle = 1, xtitle = 'Velocity [km s'+textoidl('^{-1}')+']'
          if not keyword_set(display) then begin
            ;write_png, dir+'/PLOTS_' + output_name + '/continuum_normalized/raster'+a2str(n)+'_'+a2str(x+1)+a2str(y+1)+'_normalized.png', tvrd(/true)
            device, /close
            set_plot, 'X'
          endif
          
          
          ;###########################################
          ;###########################################
          ;###########################################
          ;###########################################
          ;normalized plot
          if keyword_set(display) then set_plot, 'X' else begin
            ;set_plot, 'z'
            ;ERASE
            ;device, set_resolution=reso, set_pixel_depth = 24, decomposed = 0
            set_plot, 'PS'
            !p.font = 0
            device, filename=plotdir + '/continuum_subtracted/raster'+a2str(n)+'_'+a2str(x+1)+a2str(y+1)+'_nocontinuum.eps', /encaps, /color, /helvetica, xsize = 15, ysize = 10
          endelse
          
          cube_spectra(ind_cube).cont_subtracted_flux(0:n_elements(lambda_smooth)-1) = flux_smooth - ContinuumModel(lambda_smooth, params_cont, poly_degree = poly_degree)
          
          loadct, 0, /silent
          xtit = textoidl('\lambda_{rest}-\lambda_{0} [\mum]')
          xran = [max([0.-viewrange(0)*fwhm_um, 1.1*min(lambda_smooth, /nan)], /nan), min([viewrange(1)*fwhm_um, 1.1*max(lambda_smooth, /nan)], /nan)]
          ind = where( lambda_smooth gt xran(0) and lambda_smooth lt xran(1) )
          plot, indgen(10), /nodata, yran=[-0.2, $
            max([flux_smooth(ind), SpectrumModel(lambda_smooth, params, poly_degree = poly_degree, ncomps = n_elements(addcomponents))]-ContinuumModel(lambda_smooth, params_cont, poly_degree = poly_degree), /nan)], xran=xran, $
            xtitle = xtit, ytitle = 'Flux density (continuum-subtracted)', background = 255, color = 0, xstyle = 1+8, ymargin = [5, 5];, xmargin = [7, 7]
            
          xt = 0.02
          yt = 0.95
          if plots_publi eq 1 then begin
            st = ''; line_label(line_index) + ', RA='+a2str(0.0001*nint(10000.*cube_lres(ind_cube).ra))+', DEC='+a2str(0.0001*nint(10000.*cube_lres(ind_cube).dec))
            cs = 0.75
            coldata = 0
            colerr = 100
            colfit = 200
          endif else begin
            st = line_label(line_index) + ', Raster ' + a2str(n)+', spaxel ('+a2str(x+1)+','+a2str(y+1)+')'
            cs = 0.75
            coldata = 100
            colerr = 200
            colfit = 250
          endelse
          xyouts, xt, yt, st, /normal, alignment=0., charsize=1.25, color = 0
          velocity = cspeed * xran / lambda_rest
          
          plots, [0., 0.], !y.crange, linestyle = 0, color = 180, noclip = 0
          if plots_publi eq 0 then begin
            plots, [params[1], params[1]], !y.crange, linestyle = 2, color = 180, noclip = 0
            plots, (0.-constraints_continuum(0)*fwhm_um)*[1.,1.], !y.crange, linestyle = 0, color = 150, noclip = 0
            plots, (0.+constraints_continuum(0)*fwhm_um)*[1.,1.], !y.crange, linestyle = 0, color = 150, noclip = 0
            plots, (0.-constraints_continuum(1)*fwhm_um)*[1.,1.], !y.crange, linestyle = 0, color = 150, noclip = 0
            plots, (0.+constraints_continuum(2)*fwhm_um)*[1.,1.], !y.crange, linestyle = 0, color = 150, noclip = 0
          ;plots, [0.-params[1], 0.-params[1]], !y.crange, linestyle = 0, color = 200, noclip = 0
          endif
          
          nsigma = cube_lres(ind_cube).flux/cube_lres(ind_cube).error
          nsigma = max([nsigma, 0])
          if finite(nsigma, /nan) then nsigma = 0.
          
          ;if fitflux le fitnoise then fitstr = '<'+str_round(a2str(fitnoise))+' W m-2' else $
          
          ;fitstr = str_round(a2str(cube_lres(ind_cube).flux+cube_lres(ind_cube).localcomponent_flux))+'+/-'+str_round(a2str(cube_lres(ind_cube).error+cube_lres(ind_cube).localcomponent_error))+' W m'+textoidl('^{-2}')        ;if intflux lt intnoise then intstr = 'F_int < '+str_round(a2str(intnoise))+' W m-2' else $
          ;  intstr = 'F_int = '+str_round(a2str(intflux))+'+/-'+str_round(a2str(intnoise))+' W m-2'
          if plots_publi eq 1 then al_legend, /top, /left, box = 0, line_label(line_index), textcolor = 0 else  al_legend, /top, /right, charsize = 0.7, [fitstr, str_round(a2str(nsigma))+textoidl('\sigma'), $
            'EQW'+textoidl('_{rest}')+' = '+str_round(a2str(cube_lres(ind_cube).eqw))+ ' '+textoidl('\mu')+'m', $
            'v = '+str_round(a2str(cube_lres(ind_cube).velocity))+'+/-'+str_round(a2str(cube_lres(ind_cube).velocity_error))+' km s'+textoidl('^{-1}'), $
            'Broadening = '+str_round(a2str(cube_lres(ind_cube).broadening))+'+/-'+str_round(a2str(cube_lres(ind_cube).broadening_error))+ ' km s'+textoidl('^{-1}')], $
            textcolor = [0, 0, 0, 0, 0], box = 0
            
          loadct, 0, /silent
          
          ;oplot, lambda_smooth, replicate(1., n_elements(lambda_smooth)), color=0, thick=1, linestyle = 0
          oplot, lambda_smooth, replicate(0., n_elements(lambda_smooth)), color=0, thick=2, linestyle = 0
          
          oploterror, lambda_smooth, flux_smooth-ContinuumModel(lambda_smooth, params_cont, poly_degree = poly_degree), error_smooth, color=coldata, psym=10, symsize=0.1, errcolor = colerr
          oplot, lambda_smooth, flux_smooth-ContinuumModel(lambda_smooth, params_cont, poly_degree = poly_degree), color=coldata, psym=10, symsize=0.1, thick = 2
          ;print, flux_smooth/ContinuumModel(lambda_smooth, params_cont, poly_degree = poly_degree)
          
          loadct, 39, /silent
          oplot, lambda(indfit_smooth), SpectrumModel(lambda(indfit_smooth), params, poly_degree = poly_degree, ncomps = n_elements(addcomponents))-ContinuumModel(lambda(indfit_smooth), params_cont, poly_degree = poly_degree), color=colfit, thick=3, linestyle = 0
          oplot, lambda(indfit_smooth), (LineModel(lambda(indfit_smooth), params[0:2])+ContinuumModel(lambda(indfit_smooth), params_cont, poly_degree = poly_degree))-ContinuumModel(lambda(indfit_smooth), params_cont, poly_degree = poly_degree), color=colfit, thick=1, linestyle = 0
          if keyword_set(addcomponents) then begin
            for ncomp = 0, n_elements(addcomponents) - 1 do begin
              oplot, lambda(indfit_smooth), (LineModel(lambda(indfit_smooth), params[3:5])+ContinuumModel(lambda(indfit_smooth), params_cont, poly_degree = poly_degree))-ContinuumModel(lambda(indfit_smooth), params_cont, poly_degree = poly_degree), color=colfit, thick=1, linestyle = 0
            endfor
          endif
          
          ;plot with theoretical spectral resolution
          params_tmp = params
          params_tmp(2) = sigma
          params_tmp(5) = sigma
          oplot, lambda(indfit_smooth), SpectrumModel(lambda(indfit_smooth), params_tmp, poly_degree = poly_degree, ncomps = n_elements(addcomponents))-ContinuumModel(lambda(indfit_smooth), params_cont, poly_degree = poly_degree), color=colfit, thick = 1, linestyle = 2
          
          ;plot 3-sigma upper limit
          if nsigma lt 1. then begin
            params_tmp = params
            params_tmp(0) = 1.*pcerror(0)
            oplot, lambda(indfit_smooth), SpectrumModel(lambda(indfit_smooth), params_tmp, poly_degree = poly_degree, ncomps = n_elements(addcomponents))-ContinuumModel(lambda(indfit_smooth), params_cont, poly_degree = poly_degree), color=colfit, thick = 1, linestyle = 1
            params_tmp(0) = 4.*pcerror(0) ;that's 2 sigma, and to keep the line profile close to expected, I put 2x in the peak and 2x for sigma in the peak as well
            oplot, lambda(indfit_smooth), SpectrumModel(lambda(indfit_smooth), params_tmp, poly_degree = poly_degree, ncomps = n_elements(addcomponents))-ContinuumModel(lambda(indfit_smooth), params_cont, poly_degree = poly_degree), color=colfit, thick = 1, linestyle = 1
            params_tmp(0) = 6.*pcerror(0)
            oplot, lambda(indfit_smooth), SpectrumModel(lambda(indfit_smooth), params_tmp, poly_degree = poly_degree, ncomps = n_elements(addcomponents))-ContinuumModel(lambda(indfit_smooth), params_cont, poly_degree = poly_degree), color=colfit, thick = 1, linestyle = 1
          ;params_tmp(0) = 20.*pcerror(0)
          ;oplot, lambda(indfit_smooth), SpectrumModel(lambda(indfit_smooth), params_tmp, poly_degree = poly_degree), color=colfit, thick = 1, linestyle = 1
          endif
          
          oplot, lambda(indfit_smooth), ContinuumModel(lambda(indfit_smooth), params_cont, poly_degree = poly_degree)-ContinuumModel(lambda(indfit_smooth), params_cont, poly_degree = poly_degree), color=colfit, thick=1, linestyle = 0
          
          
          if plots_publi eq 0 then al_legend, /top, /left, charsize = 0.7, ['FWHM='+str_round(a2str(cube_lres(ind_cube).fwhm))+'+/-'+str_round(a2str(cube_lres(ind_cube).fwhm_error))+' km s'+textoidl('^{-1}'), 'FWHM_th='+str_round(a2str(fwhm_kms))+' km s'+textoidl('^{-1}')], textcolor=[250, 250], box = 0, linestyle=[0, 2], color = [250, 250], thick = [3, 1]
          
          
          
          loadct, 0, /silent
          axis, xaxis=1, xrange = velocity, color = 0, xstyle = 1, xtitle = 'Velocity [km s'+textoidl('^{-1}')+']'
          if not keyword_set(display) then begin
            ;write_png, dir+'/PLOTS_' + output_name + '/continuum_normalized/raster'+a2str(n)+'_'+a2str(x+1)+a2str(y+1)+'_normalized.png', tvrd(/true)
            device, /close
            set_plot, 'X'
          endif
          
        endfor
      ;==========================================================================================
      ;==========================================================================================
      ;==========================================================================================
        
        
        
        
        
      endif                   ;skipfit
      
      
      ;if not keyword_set(skipfit) then begin
      ;spawn, '\rm ' + dir+'/'+output_name+'_before_subtracting.sav'
      save, filename = dir+'/'+output_name+'.sav', cube_lres, cube_fitparams, cube_spectra, cube_params, cube_header
    ;endif
      
    endfor                     ;spaxel loop
    
  endfor                        ;raster loop
  
  close, /all

  if keyword_set(central3by3combine) or keyword_set(central5by5combine) then begin
    ind = where( cube_lres.spaxel_X eq (indx(0)+1) and cube_lres.spaxel_Y eq (indy(0)+1) )
    ;ind = ind(0)
    cube_lres = cube_lres(ind)
    cube_spectra = cube_spectra(ind)
    cube_fitparams = cube_fitparams(ind)
    save, filename = dir+'/'+output_name+'.sav', cube_lres, cube_fitparams, cube_spectra, cube_params, cube_header, cube_groups
  endif
  
  
  ;change dims of cube_spectra
  ;cube_spectra = cube_spectra(*, *, *, 0:n_elements(lambda_smooth)-1, *)
  ;writefits, dir+'/'+output_name+'_spectra.fits', cube_spectra, h_ini
  
  print, '-----------------------------------------------------------------------------------------'
  print, '- Fit results saved in ' + dir + '/'+line_label(line_index)+'_fitresults_'+output_name+'.txt'
  print, '- Output cubes written in ' + dir+'/'+output_name+'.sav'
  print, '- Plots saved in ' + dir+'/PLOTS_' + output_name
  print, '- and ' + dir+'/PLOTSpubli_' + output_name
  print, ''
  print, '- Example of calling sequence for PACSman_map:'
  print, "  PACSman_map, [object='objectname'], [/makeall (for batch mode)]"
  print, ''
  
end
