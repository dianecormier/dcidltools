;+
; Project     : SOHO - LASCO/EIT
;
; Name        : WMESSAGE
;
; Purpose     :
;
; Category    : Widgets
;
; Explanation :
;
; Syntax      :
;
; Examples    :
;
; Inputs      :
;
; Opt. Inputs :
;
; Outputs     :
;
; Opt. Outputs:
;
; Keywords    : TITLE (a string containing the title to be used for the widget)
;               LABEL (a string containing a message in WIDGET_LABEL)
;               TEXT  (a string to display in WIDGET_TEXT)
; Common      :
;
; Restrictions: None.
;
; Side effects: None.
;
; History     :  15-nov-1995,Borut Podlipnik, MPAe,Written
;
; Contact     : BP, borut@lasco1.mpae.gwdg.de
;-
;

PRO wmessage, title=title, text=text, label=label, xsize=xsize, ysize=ysize ,scroll=scroll


  IF (XRegistered("wmessage") NE 0) THEN RETURN
  
  IF NOT KEYWORD_SET(title) THEN title = ""
  IF NOT KEYWORD_SET(xsize) THEN xsize=35
  IF NOT KEYWORD_SET(ysize) THEN ysize=4
  IF NOT KEYWORD_SET(scroll) THEN scroll=1
  
  bv1 = "Dismiss"
  
  WYNBase = WIDGET_BASE( TITLE=title, /COLUMN )
  
  IF KEYWORD_SET(label) THEN BEGIN
    wb0 = WIDGET_BASE( WYNBase, /COLUMN )
    wl0 = WIDGET_LABEL( wb0, VALUE = label,/DYNAMIC_RESIZE )
  ENDIF
  
  IF KEYWORD_SET(text) THEN BEGIN
    wb1 = WIDGET_BASE( WYNBase, /COLUMN)
    wt1 = WIDGET_TEXT( wb1, XSIZE=xsize, YSIZE=ysize, VALUE=text,/SCROLL)
  ENDIF
  
  wb2 = WIDGET_BASE( WYNBase, /ROW)
  b1 = WIDGET_BUTTON( wb2, VALUE=bv1, XSIZE=250 )
  
  WIDGET_CONTROL, WYNBase, /REALIZE
  
  event = WIDGET_EVENT(WYNBase)
    
    WIDGET_CONTROL,event.id, GET_VALUE=yn
    WIDGET_CONTROL,WYNBase, /DESTROY
  RETURN
END