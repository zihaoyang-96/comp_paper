PRO fg_rigidalign2, index=index, data, index_outo=index_out, data_out,          $
                   DX=dx,    DY=dy,           X0=x0,    Y0=y0,$
                   NT=nt,    I0=i0,     I1=i1,                $
                   DISPVECS=dispvecs,         USERVEC=uservec,$
                   NO_APPLY=no_apply,           SMOOTH=smooth,$
                   QUIET=quiet,               VERBOSE=verbose,$
                   DISPLAY=display

;+
; NAME:
;	FG_RIGIDALIGN
; PURPOSE:
;	Calculates and applies 2-D shifts to a series of images to bring
;       them into relative alignment with the first image in the series.
;
; CALLING SEQUENCE:
;	FG_RIGIDALIGN, index, data, index_out, data_out
;
; INPUT:
;     INDEX: the SSW image header structure array for the images in DATA.
;     DATA: [nx,ny,nfiles] image time series data cube.
;
; OUTPUT:
;     [optional] INDEX_OUT: modified SSW image header structure.
;     [optional] DATA_OUT:  shifted image time series data cube.
;
; OPTIONAL KEYWORD INPUT:
;     DX: the x-dimension of the box in which image-to-image correlations are calculated.
;         Default = 256.
;     DY: the y-dimension of the box in which image-to-image correlations are calculated.
;         Default = 256.
;     X0: the x-coordinate of the lower-left corner of the correlation box.
;         Default = image center.
;     Y0: the y-coordinate of  the lower-left corner of the correlation box.
;         Default = image center.
;     NT: the number of images correlated at one time, i.e. the number of images
;         assembled into a subgroup and aligned with respect to the first image in the group.
;         Default = 4. Set this parameter to be less than the number of images over which the 
;         correlated structures in the images change significantly. 
;     I0: the series number of the image to be used as the anchor, or reference, image.
;         Default = 0.
;     I1: the series number of the last image in the set to be aligned to the reference image.
;         Default = nfiles-1.
;     USERVEC: a 2xNFILES floating point vector of user-supplied displacement vectors for the 
;         image time series. If supplied, this routine does not calculate displacements and
;         simply applies the shifts in USERVEC to the data cube.
;     SMOOTH: if supplied, smooths the displacement vectors with a boxcar
;         filter of width SMOOTH.
;     NO_APPLY: if set, suppresses application of POLY_2D shifts. In this case the routine 
;         only calculates displacements.
;     QUIET: if set, all output is suppressed.
;     VERBOSE: if set, extra diagnostic output may be printed to STDOUT.
;     DISPLAY: if set, the sub-image correlation box images are TVSCL'd to an open window.
;
; OPTIONAL KEYWORD OUTPUT:
;     DISPVECS: the 2xNFILES floating point displacement vector array. DISPVECS[0,*] contains
;         the x-displacements and DISPVECS[1,*] contains the y-displacements.
; METHOD:
;	Uses the cross-correlation reoutine TR_GET_DISP to calculate displacments of images 
;       relative to the first image in a set of size NT. The default value of NT is 4. The 
;       full dataset is broken into NFILES/NT sets (plus a remainder set). Each set is internally
;       aligned to the first image in a set by tr_get_disp. For all sets after the first one,
;       uniform displacements are added to bring the set into alignement with the first set.
;;         
; REVISION HISTORY
;	1.0 Written, T. Berger, LMSAL, 14-December-2006.
;	1.1 AdW 20081026: use /nowindow option for tr_get_disp and center
;	cut-out in the FOV by default.
;       1.2 AdW 20091006: copy index to index_out and update history.
;       - RJM commented out index update
;
;-

t0 = SYSTIME(1)
crlf = STRING(10B)+STRING(13B)

;Name, version, and timing information
prognam = 'FG_RIGIDALIGN.PRO'
progver = 'V1.2'

loud = 1 - KEYWORD_SET(quiet)
verbose = KEYWORD_SET(verbose)
if (verbose eq 1) then loud = 1
if (loud) then PRINT, 'Running ', prognam, ' ', progver
if KEYWORD_SET(display) then display = 1 else display = 0
if KEYWORD_SET(no_apply) then apply=0 else apply=1

;Determine image size and number of images
sz = SIZE(data[*,*,*])
xs = sz[1]
ys = sz[2]
nfiles=sz[3]
;nfiles = (SIZE(index))[1]

if N_ELEMENTS(i0) ne 0 then i0=i0 else i0=0
if N_ELEMENTS(i1) ne 0 then i1=i1 else i1=nfiles-1
if KEYWORD_SET(dx) then dx=dx else dx=256
if KEYWORD_SET(dy) then dy=dy else dy=256
if N_ELEMENTS(x0) ne 0 then x0=x0 else x0=xs/2-dx/2
if N_ELEMENTS(y0) ne 0 then y0=y0 else y0=ys/2-dy/2

if KEYWORD_SET(nt) then nt=nt else nt=4
ng = nfiles/nt     ;number of groups of files 
nr = nfiles mod nt ;remainder at end of data series

docalcs=1
if KEYWORD_SET(uservec) then begin
   docalcs = 0
   nvec=N_ELEMENTS(uservec)
   szuv=SIZE(uservec,/dimen)
   if szuv[0] ne 2 or szuv[1] ne nfiles then begin
      MESSAGE,'User supplied displacement vector must be 2xN where N = number of images'
      RETURN
   end
   dispvecs = uservec
end

;Calculate displacements
if (docalcs) then begin
   t0c = SYSTIME(1)
   if (loud) then MESSAGE,/info,'Calculating rigid alignments...'
   dispvecs = FLTARR(2,nfiles)
   xlo = (x0) > 0
   xhi = (x0 + dx - 1) < (xs-1)
   ylo = (y0) > 0
   yhi = (y0 + dy - 1) < (ys-1)
   dx = xhi - xlo + 1
   dy = yhi - ylo + 1
   ntfac = (nt-1) > 1           ;nt-1 does the overlap of series

   imx = 0
   for jj=i0,i1,ntfac do begin   
      if MAX(imx) eq i1 then break
      tj0 = SYSTIME(1)     
      imx = (INDGEN(nt) + jj) < i1 ;file numbers for the group data cube
      imx = imx[UNIQ(imx)]
      nf = N_ELEMENTS(imx)

      datatmp = INTARR(dx,dy,nf)
      for i=0,nf-1 do begin
         datatmp[*,*,i] = data[xlo:xhi,ylo:yhi,imx[i]]
         if (display) then begin
            TVSCL,datatmp[*,*,i]
            WAIT,0.1
         end
      end

      ;calculate displacements for the set
      tmp = tr_get_disp_2b(datatmp,/nowindow,/quiet)

      ;now add the displacement of the first image in the set
      if jj gt i0 then for k=0,nf-1 do begin

         tmp[0,k] += dispvecs[0,imx[0]]
         tmp[1,k] += dispvecs[1,imx[0]]
                                ;print,tmp[0,k],tmp[1,k]
      end
      dispvecs[0,imx] = tmp[0,*]
      dispvecs[1,imx] = tmp[1,*]     
      dt = SYSTIME(1) - tj0
      if (loud) then begin
         MESSAGE,/info,'Number of images in displacement set:'+STRTRIM(nf,2)
         MESSAGE,/info,'Time for current displacement set: '+STRTRIM(dt,2)+' seconds.'
      end

   end  ;of jj loop

   ;for i=0,nfiles-1 do print,dispvecs[0,i],dispvecs[1,i]
   
   datatmp=0
   dt = (SYSTIME(1)-t0c)
   if (loud) then MESSAGE,/info,crlf+'     Total time for rigid alignment calculation: '+STRTRIM(dt,2)+' seconds.'

   if KEYWORD_SET(smooth) then begin
      sfac = ABS(smooth)
      dispvecs[0,*] = SMOOTH(dispvecs[0,*],sfac)
      dispvecs[1,*] = SMOOTH(dispvecs[1,*],sfac) 
   end

end  ;of docalcs block


;Apply shifts:
if (apply) then begin 
   ;index_out = index
   ;update_history, index_out, 'Rigid alignment applied', caller=prognam, version=version
   ts0 = SYSTIME(1)
   data_out = data[*,*,0]
   for kk=0,nfiles-1 do begin
      im = data[*,*,kk]
      im = POLY_2D(im,[-dispvecs(0,kk),0.,1.,0.],$
                   [-dispvecs(1,kk),1.,0.,0.],cubic=-0.5)
    
    ;Alternative shift function - RJM
    ; im=shift_frac2d(im,dispvecs(0,kk),dispvecs(1,kk))
      data_out = [ [[[data_out]]],[[im]] ]
   end
   data_out = data_out[*,*,1:*]
   dt = SYSTIME(1)-ts0
   if (loud) then MESSAGE,/INFO,'Total time for rigid alignment application = '+STRTRIM(dt,2)+' seconds'
end

dt = +SYSTIME(1)-t0
if (loud) then MESSAGE,/info,'Total time for rigid alignment = '+STRTRIM(dt,2)+' seconds'


RETURN
END
