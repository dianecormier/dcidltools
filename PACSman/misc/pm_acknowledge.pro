;#> acknowledge.dc3
;
; Identifier    acknowledge 
;
; Purpose      A function to pop up a widget with error or warning.
;
; Synopsis     acknowledge, text=text 
;
; Arguments     Name      I/O  Type        Description
;               --------------------------------------------------------
;               text      I    strarr     message to be displayed
;           
; Returns          
;
; Description 
;            
;
; Dependencies  Calls: MODAL_BASE
;               Called By: ISAP GUI MODULES, SAP_ERROR
;                 
; Comment      
;
; Example  
;
; Category   ISAP 
;
; Filename   acknowledge.pro
;
; Author     M. Buckley 
;
; Version     1.4 
;
; History     1.0   1996      initial design 
;             1.1  09-JUL-96  Changed so that hourglass is put up after operation
;             1.3  13-AUG-96  changed title to acknowledge (SJU) 
;             1.4   3-MAR-98  removed map option for IDL 5 modal problem (RN) 
;             1.5   6-May-05  Peter Hall, Cornell University.
;                             Added "parent_group" and "modal" processing.
;
; Copyright (C) 1997, California Institute of Technology.
;*
;#<

;******************************************************************************

pro pm_acknowledge_event,ev
widget_control,ev.top,/destroy				; destroy the modal widget
end

;*****************************************************************************

pro pm_acknowledge,text = text,$							; text to display
		tlb = tlb,$							; top level base to orient to
		ack_string = ack_string,$					; string to display on the accept button
		position = position,$						; string to give placement guide
		title = title,$
;		nowid = nowid							; flag, if set, just print message/s
		nowid = nowid, $	                        		; flag, if set, just print message/s ; PH 26 April 2005
                modal=modal, $                                                  ; PH 26 April 2005
                parent_group=parent_group                                       ; PH 26 April 2005


text = string(text)
if not keyword_set(ack_string) then ack_string = '    ok    '		;set up the button strings, if none s

if keyword_set(nowid) then begin
   print,'ACKNOWLEDGE:',text
   return								; return
endif

if not keyword_set(title) then title = 'Acknowledge'			; set up title, if none given

;ack_base = modal_base(/column,title = title); create the base; PH 26 April 2005

    if keyword_set(parent_group) then begin                                                  ; PH 26 April 2005
        modal = 1                                                                            ; PH 26 April 2005
        ack_base = widget_base(/column,title = title, modal=modal, group_leader=parent_group) ; PH 26 April 2005
    endif else begin                                                                         ; PH 26 April 2005
        ack_base = widget_base(/column,title = title)                                         ; PH 26 April 2005
    endelse                                                                                  ; PH 26 April 2005

   if keyword_set(text) then begin						; if we have been given some text
      xsiz = max(strlen(text))							; find the longest string in the input
      siz = size(text)								; find the SIZE info for the input
      if siz(0) gt 0 then ysiz = siz(1) else ysiz = 1				; if it is an array, use num of dims, else use 1
;     junk = widget_text(ack_base,value = text,xsize = xsiz,ysize = ysiz)	; set up the text wid
      junk = widget_text(ack_base,xsize = xsiz,ysize = ysiz) ; KS 15.7.98
   endif
   horiz_base = widget_base(ack_base,/row)					; make a nice base for the buttons
   acceptwid = widget_button(horiz_base,value = ack_string,uval = 'accept')	; set up the accept button

; (out in 1.5) widget_control,ack_base,map = 0							; make invisible
widget_control,ack_base,/realize						; realize it
position_base,ack_base,tlb=tlb,xpos,ypos,position = position			; position it on the tlb, if we were given one
widget_control,ack_base,tlb_set_xoffset = xpos,tlb_set_yoffset = ypos		; set its position
widget_control,junk,set_value=text ; KS 15.7.98
; (out in 1.5) widget_control,ack_base,map = 1							; make it visible
widget_control,acceptwid,/input_focus						; put the cursor on the text widget
xmanager,'pm_acknowledge',ack_base
widget_control,/hourglass							; put the hourglass up
end
