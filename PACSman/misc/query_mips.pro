;+
; NAME:
;   queryIRAC
;
; PURPOSE:
;    Query the IRSA server and return the first Spitzer/MIPS image found for a
;    given object
;
; EXPLANATION:
;     The script can query andd retrieve an image and FITS
;     header from the Infrared Science Archive (IRSA) at the Infrared
;     Processing & Analysis Center (IPAC).
;     See http://irsa.ipac.caltech.edu for details.
;
; CALLING SEQUENCE:
;      QueryMIPS, object_name, Im, Hdr, BAND=
;
; INPUTS:
;      object_name - object name
;
;      band - MIPS band: '24um', '70um', '160um'
;
; OPTIONAL INPUTS: None
;
; OUTPUTS:
;       Im - The image returned by the server. If there is an error, this
;             contains a single 0.
;
;       Hdr - The FITS header of the image. Empty string in case of errors.
;
; RESTRICTIONS:
;      Relies on a working network connection.
;
; PROCEDURE:
;      Construct a query-url,  call WEBGET() and sort out the server's
;      answer.
;
; EXAMPLE:
;
; PROCEDURES CALLED:
;       QUERYSIMBAD, WEBGET()
; MODIFICATION HISTORY:
;       Written by V. Lebouteiller, CEA/SAp, 2010
;
;-
;
;+
; NAME:
;    WEBGET()
;
; PURPOSE:
;    Use the IDL SOCKET procedure to get data from http servers
;
; EXPLANATION:
;     WEBGET() can access http servers - even from behind a firewall -
;     and perform simple downloads. Currently, text and FITS files can be
;     accessed.
;
; CALLING SEQUENCE:
;      a=webget(URL)
;
; INPUTS:
;      URL - scalar string giving a fully qualified url of the form
;          'http://server.eso.org/path/file.html'.    WEBGET() can
;          also use other valid URLs that contain 'GET' or 'POST' codes.
;
; OPTIONAL INPUT KEYWORD PARAMETERS:
;       COPYFILE - if set to a valid filename (file must have write permission),
;            the data contents of the web server's answer is copied to that
;            file.
;       HTTP10 - If set, then use the HTTP 1.0
;       POST - if set to a structure, the structure tags and values
;              will be used as post variables and POST'ed to the URL.
;              If POST is not set, the normal HTTP GET is used to
;              retrieve the URL.
;       /SILENT - If set, the information error messages are suppressed
; OUTPUTS: A structure with the following fields:
;
;            .Header - the HTTP header sent by the server
;
;            .Text   - The text part of the downloaded file. If the
;                     content type of the file was not of class
;                     'text',  this will be an empty string.
;
;            .ImageHeader - Header file of a FITS-image. FITS images
;                          are read when the content type is
;                          'image/fits' or 'application/octet-stream'
;                          (for dss-access). If the file is not a FITS
;                          image,  this will be an empty string.
;
;            .Image - The FITS image read from the server. If the file
;                    did not contain a FITS image,  this will be zero.
;
;
; RESTRICTIONS:
;     The mime-type recognition is extremely limited. Only the content-type is
;     determined. Any text-file  will be stored in out.Text. The only other
;     category which can be fetched is FITS files,  which will be stored in
;     out.Image and out.ImageHeader.
;
;     PROXY: If you are behind a firewall and have to access the net through a
;         Web proxy,  set the environment variable 'http_proxy' to point to
;         your proxy server and port, e.g.
;         'setenv http_proxy=http://web-proxy.mpia-hd.mpg.de:3128'
;
;               The URL *MUST* begin with "http://".
;
; PROCEDURE:
;     Open a socket to the webserver and download the header. After deciding
;     whether it is text or binary, either store the text or try to read a
;     FITS file.
;
; EXAMPLE:
;      IDL> a=webget('http://www.mpia.de/index.html')
;      IDL> print,a.Text
;      or
;
;          > PointingRA=0.0
;          > PointingDE=30.0
;          > QueryURL = strcompress("http://archive.eso.org/dss/dss/image?ra="+$
;          >                          string(PointingRA)+$
;          >                          "&dec="+$
;          >                          string(PointingDE)+$
;          >                          "&x=10&y=10&Sky-Survey=DSS1&mime-type=download-fits", $
;          >                          /remove)
;          > a=webget(QueryURL)
;          > tvscl,a.Image
;          > print,a.ImageHead
;
;
; MODIFICATION HISTORY:
;     Written by M. Feldt, Heidelberg, Oct 2001 <mfeldt@mpia.de>
;     Use /swap_if_little_endian keyword to SOCKET  W. Landsman August 2002
;     Less restrictive search on Content-Type   W. Landsman   April 2003
;     Modified to work with FIRST image server-  A. Barth, Nov 2006
;     Better recovery from errors  W. Landsman  April 2007
;     Add support for POST access               J.D. Smith    June 2007
;     Recognize "fits" image type used by SKYVIEW   W. Landsman  June 2007
;     Upgraded, partially, to HTTP 1.1				M. Perrin, July 2007
;       The HTTP 1.1 support is presently INCOMPLETE: virtual servers are
;       supported, but chunked transfer encoding is not yet supported, so
;       technically this is not fully HTTP 1.1 compliant.
;     Added http10 keyword  W. Landsman   August 2007
;     Assume since V5.6, sockets always available  W. Landsman Nov 2007
;-

PRO MimeType,  Header, Class, Type, Length
  ;;
  ;; MIME type recognition
  ;
  Class = 'text'
  Type = 'simple'               ; in case no information found...
  def = strupcase(strmid(header,0,13))
  g = where(def EQ 'CONTENT-TYPE:', Ng)
  if Ng GT 0 then begin
    ClassAndType = strmid(Header[g[0]], 14, strlen(Header[g[0]])-1)
    Class = (strsplit(ClassAndType, '/', /extract))[0]
    Type = (strsplit(ClassAndType, '/', /extract))[1]
  ENDIF
  def = strupcase(strmid(header,0,15))
  g = where(def EQ 'CONTENT-LENGTH:', Ng)
  if Ng GT 0 then $
    Length = long(strmid(Header[g[0]], 15, strlen(Header[g[0]])-1))
  return
END

FUNCTION sm_webget,  url,  SILENT=silent, COPYFILE=copyfile, POST=post, $
    HTTP10=http10, timeout = timeout, status = status, force_text = force_text
  compile_opt idl2
  ;;
  ;;
  ;; define the result fields
  ;;
  Header = strarr(256)
  Data = strarr(256)
  Image = 0
  ImageHeader = ''
  
  ;; Setup post variables
  if n_elements(post) ne 0 then begin
    method='POST'
    t=tag_names(post)
    post_vars=strarr(n_elements(t))
    for i=0,n_elements(t)-1 do $
      post_vars[i]=strlowcase(t[i])+'='+strtrim(post.(i),2)
    post_vars=strjoin(post_vars,'&')
    post_data=['Content-Type: application/x-www-form-urlencoded',$
      'Content-Length: '+strtrim(strlen(post_vars),2), $
      '', $
      post_vars]
  endif else method='GET'
  
  
  ;;
  ;; open the connection and request the file
  ;;
  ProtocolString = keyword_set(http10) ? "HTTP/1.0" : " HTTP/1.1"
  UserAgentString= "IDL "+!version.release+' on '+!VERSION.OS+'/'+!VERSION.ARCH
  Proxy = getenv('http_proxy')
  IF Proxy NE '' THEN BEGIN
    ;;
    ;; sort out proxy name
    ;;
    LastColon = StrPos(Proxy, ':', /Reverse_Search)
    ProxyPort = fix(StrMid(Proxy, LastColon+1))
    ProxyServer = StrMid(Proxy, 7, LastColon-7)
    ;; open the connection and send the 'GET' command
    socket, unit, ProxyServer,  ProxyPort, /get_lun, /swap_if_little_endian, read_timeout = timeout, connect_timeout = timeout, write_timeout = timeout, error = error
    status = error
    if status ne 0 then return, ''
    printf, unit, method+' '+url+ProtocolString
  ENDIF ELSE BEGIN
    ;;
    ;; same thing easier without proxy
    ;;
    slash1 = StrPos(strmid(url, 7), '/')
    Server = StrMid(url, 7, slash1 )
    purl = strmid(url,slash1+7)
    Port = 80
    socket, unit, Server,  Port, /get_lun,/swap_if_little_endian, read_timeout = timeout, connect_timeout = timeout, write_timeout = timeout, error = error
    status = error
    if status ne 0 then return, ''
    printf, unit, method+' '+purl + ProtocolString
  ENDELSE
  ;; These lines are the same for either with or without proxy.
  ;; in HTTP 1.1 we MUST include the Host: line to allow requests
  ;; from co-hosted virtual servers to operate properly.
  printf, unit, "Host: "+Server
  printf, unit, 'User-Agent: '+ UserAgentString
  ;; HTTP 1.1 clients must either support persistent connections, or indicate
  ;; they do not by stating Connection: close
  printf, unit, "Connection: close"
  
  ;; Add the POST data, if requested
  if n_elements(post) ne 0 then printf,unit,transpose(post_data)
  ;; Blank line required to terminate HTTP request.
  printf, unit, ''
  
  LinesRead = 0
  text = 'xxx'
  ;;
  ;; now read the header
  ;;
  On_IOERROR, done
  WHILE  text NE '' do begin
    readf, unit, text
    Header[LinesRead] = text
    LinesRead = LinesRead+1
    IF LinesRead MOD 256 EQ 0 THEN $
      Header=[Header, StrArr(256)]
  ENDWHILE
  DONE: On_IOERROR, NULL
  ;;
  if LinesRead EQ 0 then begin
    message,'Unable to read HTTP server - Check timeout',/CON
    free_lun,unit
    return,{Header:'', Text:'', ImageHeader:ImageHeader,  Image: Image}
  endif
  Header = Header[0:LinesRead-1]
  MimeType, Header, Class,  Type, Length; analyze the header
  ;;
  IF Keyword_Set(CopyFile) THEN BEGIN
    openw, wunit, CopyFile, /get_lun
    aaa = bytarr(Length,/nozero)
    readu, unit, aaa
    writeu, wunit, aaa
    free_lun, wunit
    free_lun, unit
    return, 1
  ENDIF
  ;;
  text = '' ;initialize text fields
  LinesRead = 0l
  ;;
  
  if keyword_set(force_text) then class = 'text'
  
  CASE Class OF
    'text': BEGIN
      ;;
      ;; read anything of class 'text'
      WHILE  eof(unit) EQ 0 do begin
        readf, unit, text
        Data[LinesRead] = text
        LinesRead = LinesRead+1
        IF LinesRead MOD 256 EQ 0 THEN $
          Data=[Data, StrArr(256)]
      ENDWHILE
      if LinesRead EQ 0 then message,'ERROR - no lines of text read',/CON
      Data = Data[0:(LinesRead-1) > 0 ]
    END
    'image': BEGIN
      CASE Type OF
        'x-fits': Image = readfits(unit, ImageHeader)
        'fits': Image = readfits(unit, ImageHeader)
        else: message,'Unrecognized  image type of ' + type
      ENDCASE
    END
    'application':BEGIN
    CASE Type OF
      'octet-stream': BEGIN ; try reading a FITS file because ESO
        ; answers this way
        Image = readfits(unit, ImageHeader)
      END
      'octet-stream': BEGIN ; try reading a FITS file because ESO
        ; answers this way
        Image = readfits(unit, ImageHeader)
      END
      'fits': BEGIN     ; need this for FIRST survey
        Image = readfits(unit, imageheader)
      END
      
    ENDCASE
  END
ENDCASE

IF LinesRead EQ 0 THEN Data = ''
free_lun, unit
return, {Header:Header, Text:Data, ImageHeader:ImageHeader,  Image: Image}
END

;---------------------------------------------

pro query_mips, object, image, header, band = band, status = status

  status = 1
  
  case band of
    '24um': mipsim = '1'
    '70um': mipsim = '2'
    '160um': mipsim = '3'
  endcase
  
  print, 'Contacting Spitzer@IRSA... Please wait'
  q = sm_webget("http://irsa.ipac.caltech.edu/cgi-bin/Spitzer/Spitzer/nph-spitzer?mipsim="+mipsim+"&locstr="+object+"&mode=prog")
  r = stregex(q.text, 'http.*mipsim.*bat', /extract)
  ind = where(r ne '', c)
  if c ge 1 then begin
    r = (r[ind])[0]
    q = sm_webget(r, /force_text)
    print, 'Query result:'
    print, q
    r = stregex(q.text, 'http.*"', /extract)
    ind = where(r ne '', c)
    if c ge 1 then begin
      r = r[ind]
      valid = 0
      for i = 0, n_elements(r) - 1 do begin
        r2 = strmid(r[i], 0, strlen(r[i])-1)
        r2 = strtrans(r2, ':80', '')
        print, r2
        result = sm_webget(r2)
        if strpos(result.text[0], 'NO DATA') eq -1 then begin
          valid = 1
          break
        endif
      endfor
      if valid eq 0 then begin
        print, 'No result'
      endif else begin
        print, 'MIPS image retrieval:'
        print, result.header
        image = result.image
        header = result.imageheader
      endelse
    endif else begin
      print, 'no result'
      status = 0
    endelse
  endif else  begin
    print, 'no result'
    status = 0
  endelse
end
