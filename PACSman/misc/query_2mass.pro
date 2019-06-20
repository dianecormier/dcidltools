PRO query_2mass, target, Image,  Header, IMSIZE=ImSIze, NED=ned, BAND = band
;+
; NAME: 
;   Query2MASS
;
; PURPOSE: 
;    Query the 2MASS Image Catalog on-line at IRSA server and return
;    the 2MASS Quicklook images which are NOT recommended for
;    photometric measurements.  The returned images are the full scan
;    data.
;
; EXPLANATION: 
;     The script can query the 2MASS survey and retrieve an image and FITS 
;     header from the Infrared Science Archive (IRSA) at the Infrared
;     Processing & Analysis Center (IPAC).
;     See http://irsa.ipac.caltech.edu for details.
;
; CALLING SEQUENCE: 
;      Query2MASS, targetname_or_coords, Im, Hdr, [IMSIZE= , BAND= , /NED]
;
; INPUTS:
;      TARGETNAME_OR_COORDS - Either a scalar string giving a target name, 
;          (with J2000 coordinates determined by SIMBAD), or a 2-element
;          numeric vector giving the J2000 right ascension in *degrees* and
;          the target declination in degrees.
;
; OPTIONAL INPUTS: None
;
;
; OPTIONAL KEYWORD PARAMETERS: 
;     ImSize - Numeric scalar giving size of the image to be retrieved in 
;                 arcminutes.  Maximum allowed is 300.
;
;     /NED - Query the Nasa Extragalactic Database (NED) for the
;            target's coordinates.  The default is to use Simbad for
;            the target search.
;
;     BAND - Scalar string specifying which survey to retrieve.  
;          Possible values are 
;          'J'  - J-band (1.24 microns), this is the default
;          'H'  - H-band (1.66 microns)
;          'K'  - K_s band (2.16 microns)
; 
; OUTPUTS: 
;       Im - The image returned by the server. If there is an error, this 
;             contains a single 0.
;
;       Hdr - The FITS header of the image. Empty string in case of errors.
;
; SIDE EFFECTS: 
;     If Im and Hdr exist in advance,  they are overwritten.
;
; RESTRICTIONS: 
;      Relies on a working network connection. 
;
; PROCEDURE: 
;      Construct a query-url,  call WEBGET() and sort out the server's 
;      answer.
;
; EXAMPLE:           
;      Retrieve an 10' image surrounding the ultracompact HII region
;       G45.45+0.06.   Obtain the J-band image.
;
;        > Query2MASS, 'GAL045.45+00.06', image, header, imsize=10, band = 'J'
;        > tvscl, image
;        > hprint, header
;        > writefits,'jband_image.fits', image, header
; Note that the coordinates could have been specified directly, rather than
; giving the target name.
;        > Query2MASS, [288.587, 11.1510], image, header,imsize=10., band='J'
;
; PROCEDURES CALLED:
;       QUERYSIMBAD, WEBGET()
; MODIFICATION HISTORY: 
;       Written by J. Brauher, IPAC Apr 2003
;         Used QueryDSS as example for Query2MASS
;
;-
  if N_params() LT 2 then begin
      print,'query2mass : Syntax - Query2MASS, TargetName_or_coords, image, header' ; DBL 7/14/2006
      print,"query2mass :            [Imsize= , /ned, band = ['J','H','K'] ]" ; DBL 7/14/2006
      return
   endif

  if N_elements(target) EQ 2 then begin
      ra = float(target[0])
      dec = float(target[1])
  endif else begin
    if not (keyword_set(ned)) then begin
       QuerySimbad, target, ra,dec, Found = Found
       if found EQ 0 then begin 
         message,'Target name ' + target + $
                 ' could not be translated by SIMBAD', /informational
         return
       endif
    endif else begin
       QuerySimbad, target, ra,dec, /ned, Found = Found
       if found EQ 0 then begin
         message,'Target name ' + target + $
                 ' could not be translated by NED', /informational
         return
       endif
    endelse
  endelse  

  IF NOT Keyword_Set(ImSize) THEN ImSize = 10.
  Equinox = 'j2000'

;- Convert ImSize to arcseconds
  ImSize = Imsize * 60.

  IF N_elements(band) EQ 0 then band = 'J'

  QueryURL='http://irsadev.ipac.caltech.edu/cgi-bin/Oasis/' + $
           '2MASSImg/nph-2massimg?objstr=' + $
           strcompress(string(ra),/remove_all) + 'd+' + $
           strcompress(string(dec),/remove_all) + 'd+' + $
           'eq+' + equinox + '&size=' + $
           strcompress(string(imsize),/remove_all) + $
           '&band=' + $
           strcompress(band,/remove_all) + '&zipped=0'

  Result = webget(QueryURL)
  Image = Result.Image
  Header = Result.ImageHeader

  IF N_Elements(Image) NE 1 THEN return
  message, 'Problem retrieving your image! The server answered:', /info
  print, 'query2mass : Result.Text = ',Result.Text ; DBL 7/14/2006

END 
