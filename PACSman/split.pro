pro split, filename

  if not keyword_set(filename) then begin
    print, "Syntax: split, 'file_ortho.fits'"
    spawn, '\ls *_ortho.fits'
    retall
  endif
  if file_test(filename) eq 0 then begin
    print, 'File does not exist: ', filename
    retall
  endif
  
  test = readfits(filename, h, /silent)
  sxdelpar, h, 'NAXIS3'
  sxaddpar, h, 'NAXIS', 2
  writefits, strtrans(filename, '.fits', '_flux.fits'), reform(test(*, *, 0)), h
  writefits, strtrans(filename, '.fits', '_error.fits'), reform(test(*, *, 1)), h
  writefits, strtrans(filename, '.fits', '_velo.fits'), reform(test(*, *, 2)), h
  writefits, strtrans(filename, '.fits', '_velo_error.fits'), reform(test(*, *, 3)), h
  writefits, strtrans(filename, '.fits', '_fwhm.fits'), reform(test(*, *, 4)), h
  writefits, strtrans(filename, '.fits', '_fwhm_error.fits'), reform(test(*, *, 5)), h
  writefits, strtrans(filename, '.fits', '_continuum.fits'), reform(test(*, *, 6)), h
  writefits, strtrans(filename, '.fits', '_continuum_error.fits'), reform(test(*, *, 7)), h
  ;bigresult = [[[nresult]], [[nnoise]], [[velomap]], [[velomap_error]], [[fwhmmap]], [[fwhmmap_error]], [[contmap]], [[contmap_error]]]
  
  
  
  
  
end
