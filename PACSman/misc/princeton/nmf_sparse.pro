;+
; NAME:
;   nmf_sparse
; PURPOSE:
;   non-negative PCA routine, sparse version
; CALLING SEQUENCE:
;   nmf_sparse, data, data_ivar, ncomp, mmatrix, tol, templates=,
;   coeffs=
; INPUTS:
;   data - sparse matrix struct (see below) [Nd, Nk] of data
;   data_ivar - sparse matrix struct (see below) [Nd, Nk] of inverse var
;   ncomp - number of basis vectors sought
;   mmatrix - [Nd, Nf] transformation from desired basis to observed coords
;             (observed coords are in Nd dimensions, desired basis is
;             in Nf dimensions)
;   tol - tolerance of fit
; INPUT/OUTPUTS:
;   templates - [nf, ncomp] initial guess at templates
;   coeffs - [ncomp, nk] initial guess at coefficients
; COMMENTS:
;   Nk are the number of data points
;   Nd are the number of types of measurements
;   Nf are the number of model components
;   
;   The sparse matrix structure referred to above is:
;       .VAL[NVAL]      - actual values in matrix
;       .X[NVAL]        - columns for each value in matrix
;       .NX             - number of columns
;       .NY             - number of rows
;       .ROWSTART[NY]   - starting position of each row in VAL, X
;       .NXROW[NY]      - number of columns in each now 
; REVISION HISTORY:
;   2005-Feb-5  Written by Mike Blanton, NYU
;               Adapted from Matlab code of Sam Roweis
;
;----------------------------------------------------------------------
pro nmf_sparse, data, data_ivar, ncomp, mmatrix, tol, coeffs=coeffs, $
                templates=templates

if(n_elements(tol) eq 0) then tol=100

;; check nonneg
ilezero=where(data.val le 0., nlezero)
if(nlezero gt 0) then begin
    message, 'nmf_sparse requires data to be all-positive'
endif

;; check dims
nd=data.nx
nk=data.ny
nd2=(size(mmatrix,/dim))[0]
nf=(size(mmatrix,/dim))[1]
splog, 'number of dimensions: '+strtrim(string(nd),2)
splog, 'number of data points: '+strtrim(string(nk),2)
splog, 'number of basis vectors sought: '+strtrim(string(ncomp),2)
if(nd ne nd2) then begin
    message, 'mmatrix dimensions do not match data dimensions'
endif

;; create datap=data*data_ivar
datap=data
datahatp=data
datahat=data
datap.val=data.val*data_ivar.val

;; initialize templates and coeffs
if(keyword_set(templates) eq 0) then $
  templates=randomu(seed, nf, ncomp)+0.5  ;; HACK: scale hardwired 
templates=templates/((dblarr(nf,1)+1.)#total(templates,1))
if(keyword_set(coeffs) eq 0) then begin
;;    mdata=transpose(mmatrix)#data
    mmsparse, mdata, mmatrix, data
    coeffs=transpose(templates)#mdata
endif
help,templates,coeffs

;; make first guess model
mcoeffs=templates#coeffs
mmeval, datahat, transpose(mmatrix), mcoeffs
err=total((datahat.val-data.val)^2*data_ivar.val,/double)
splog, 'initial error= '+strtrim(string(sqrt(err/(nd*nk))),2)

err=1.d+99
eold=1.d+100
iters=1

if(tol gt 1.) then begin
    maxiters=tol
    tol=0 
endif else begin
    maxiters=1000000000L
endelse

while(iters le maxiters and abs(err-eold)/err gt tol) do begin
    datahatp.val=datahat.val*data_ivar.val
;;    mdatap=transpose(mmatrix)#datap
    mmsparse, mdatap, mmatrix, datap
;;    mdatahatp=transpose(mmatrix)#datahatp
    mmsparse, mdatahatp, mmatrix, datahatp
    coeffs=coeffs*((transpose(templates)#mdatap)/ $
                   (transpose(templates)#mdatahatp))
    templates=templates*((mdatap#transpose(coeffs))/ $
                         (mdatahatp#transpose(coeffs)))
    templates=templates/((dblarr(nf,1)+1.)#total(templates,1))
    mcoeffs=templates#coeffs
;;  datahat=mmatrix#mcoeffs
    mmeval, datahat, transpose(mmatrix), mcoeffs

    eold=err
    err=total((datahat.val-data.val)^2*data_ivar.val,/double)
    splog, 'error at iter '+strtrim(string(iters),2)+' = '+ $
      strtrim(string(sqrt(err/(nd*nk))),2)
    iters=iters+1
endwhile

end
