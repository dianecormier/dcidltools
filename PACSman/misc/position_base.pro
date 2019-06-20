;#> position_base.dc2
; Identifier	POSITION_BASE
;
; Purpose	to calculate the on-screen position of a given realized
;		widget
;
; Synopsis	pro position_base,
;			wid,
;			x,
;			y,
;			tlb = tlb,
;			position = position
;
;
; Arguments	Name		i/o type	description
;               --------------------------------------------------------
;		wid		in		widget id
;		x		out		outgoing x position
;		y		out		outgoing y position
;		position	in,optional	string to set position type
;
; Returns	none
;
; Dependencies   Calls:
;                Called By: ACKNOWLEDGE, QUERY_OK, STR_QUERY
;
; Category       ISAP
;
; Filename       position_base.pro 
;
; Author	Mark Buckley, RAL
;
; Version	1.0
;
; History
;
;*************************************************************
; Copyright 1996, Rutherford-Appleton Laboratory, UK
;*************************************************************
;#<

pro position_base,wid,$				; id of the widget base to centre (IN)
		  x,$				; return the x coord to place the widget (OUT)
		  y,$				; return the y coord to place the widget (OUT)
		  tlb = tlb,$			; id of the top level widget in which to centre it (IN,OPTIONAL)
		  position = position		; string to hold position type (IN,OPTIONAL)


if not keyword_set(position) then position = 'centre'

if keyword_set(tlb) then widget_control,tlb,tlb_get_offset = tlb_off,tlb_get_size = tlb_siz $
else device,get_screen_size = scrsiz

widget_control,wid,tlb_get_size = wid_siz				; get the size of the modal widget

x = 0
y = 0

case position of
'centre'   : if keyword_set(tlb) then begin
	       x = tlb_off(0) + tlb_siz(0)/2 - wid_siz(0)/2
	       y = tlb_off(1) + tlb_siz(1)/2 - wid_siz(1)/2
	     endif else begin
	        x = scrsiz(0)/2 - wid_siz(0)/2
		y = scrsiz(1)/2 - wid_siz(1)/2
	     endelse

'top-centre': if keyword_set(tlb) then begin
	       x = tlb_off(0) + tlb_siz(0)/2 - wid_siz(0)/2
	       y = tlb_off(1)
	     endif else begin
	        x = scrsiz(0)/2 - wid_siz(0)/2
		y = 0
	     endelse

'top-left' : if keyword_set(tlb) then begin
	        x = tlb_off(0)					; take the topleft x position
	        y = tlb_off(1)					; take the topleft y position
	     endif else begin
		x = 0
		y = 0
	    endelse

'bottom-left' : if keyword_set(tlb) then begin
		   x = tlb_off(0)
		   y = tlb_off(1) + tlb_siz(1) - wid_siz(1)
                endif else begin
		   x = 0
		   y = scrsiz(1) - wid_siz(1)
                endelse

'top-right' : if keyword_set(tlb) then begin
		 x = tlb_off(0) + tlb_siz(0)				; find top right x position
		 y = tlb_off(1)						; find top right y position
		 x = x - wid_siz(0)
	      endif else begin
		 x = scrsiz(0) - wid_siz(0)
		 y = 0
	      endelse

'bottom-right' : if keyword_set(tlb) then begin
		    x = tlb_off(0) + tlb_siz(0)
		    y = tlb_off(1) + tlb_siz(1)
		    x = x - wid_siz(0)
		    y = y - wid_siz(1)
		 endif else begin
		    x = scrsiz(0) - wid_siz(0)
		    y = scrsiz(1) - wid_siz(1)
		 endelse

else: print,'POSITION BASE : this positioning option not available : ',position
endcase
end
