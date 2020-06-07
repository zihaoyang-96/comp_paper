;+
; NAME:
;     TR_GET_DISP_2b.PRO
;
; PURPOSE:
; Given a datacube of images, measure the rigid displacement
; between each image and the first in the cube and optionally
; shift each image to coalign the entire cube.  Uses a centered
; square window for coalignment, side=2^n or 3*2^n.  Finds the whole-pixel
; shift of maximum cross-correlation, then interpolates for
; fractional pixel part either on cross-correlation function or
; optionally on an array of squared mean absolute deviations (MAD).  
; With MAD option, shifts are accurate to better than 0.1 pixel for 
; test images which are truly identical except for displacement, as long
;   as displacements are ~< 25 % of the square window
;
; CALLING SEQUENCE:
; disp = tr_get_disp(data [,/shift][,mad=mad][,/debug])
;
; INPUTS:
; data -- data cube of images, size (NX,NY,NT); 
;         data(*,*,0)= reference image to which others are aligned
;
; KEYWORD PARAMETERS:
; shift -- if set, data(*,*,1:*) are shifted to match the reference
; mad -- if non-zero, use mad x mad array of mean absolute difference (MAD) residuals to
;         find final displacement; if 1 < mad < 5, uses 5x5 array
;         mad has to have integer value!!!
;         USE ONLY IF IMAGES ARE SAME WAVELENGTHS AND EXPOSURE LEVELS
; debug  -- prints out parts of cross-correlation & MAD arrays
; themax - returns values maximum values of coherence
; quiet -- suppresses error output when MAD fails
; nomad -- counts the number of MAD fails
;
; OUTPUTS:
; Returns array of displacements disp = fltarr(2,NT)
; The sense is that data(i,j,0) <==> data(i-disp(0,k),j-disp(1,k),k)
;
; You can shift the images afterwards using data = shift_img(data, disp)
;
; TO DO:
; Needs better handling when mad x mad doesn't include the minimum - Partially fixed by RJM (program doesn't
; crash anyway!)
;
; MODIFICATION HISTORY:
; 29-Jan-98 (TDT) - adapted get_disp from H. Lin's flat fielding package ccdcal5.pro
;  1-Jul-98 (TDT) - added shift keyword, shift_img using poly_2d
; 11-Sep-98 (TDT) - variable MAD search area, fixed bug in subarea size 
; 22-Sep-98 (TDT) - added cross-correlation only feature
;  1-Oct-98 (TDT) - renamed tr_get_disp and put on-line
; 12-Mar-02 (TDT) - fixed bug so MAD now uses approx. same FOV as CC
; 12-Mar-02 (TDT) - changed to 2-D least square fitting to find minimum of CC and MAD
;  10-Jun-12 (RJM) - Edited to include change made by 26-Oct-08 (AdW) - added nowindow option to align on the whole frame
;  10-Jun-12 (RJM) - Added if statement to stop routine crashing during mad calculation, 
;                    i.e., if MAD cannot be performed takes whole pixel shift value
;  13-Mar-13 (RJM) - Added keyword themax which will return the maximum coherence value
;  Nov-13 (RJM) - Added keywords /nomad to return number of failed mads and /quiet to supress error string 
;  Jun-14 (RJM) - updated the max to give normalised values of CC -CAUTION -may not yet work for apodized images

;  tr_get_disp  get the image displacements
;
;  Method: Correlation tracks the image sequence using a 
;  square area centered on the image(s).  First image of sequence
;  is the reference.  Returns array of pixel displacements of images
;  with respect to first image.  May fail if the displacement > half of
;  the square area used for tracking
; The sense is that data(i,j,0) <==> data(i-disp(0,k),j-disp(1,k),k)
; Changed by TDT to return fractional pixel offsets,
;   added MAD algorithm  29-Jan-98
;   added shift keyword, 1-Jul-98:  if set, shifts all images to match img(*,*,0)
;   variable MAD search area (def = 5x5), fixed bug in subarea size, 11-Sep-98





;  is_in_range    true where x is inside the interval [lo,hi]
;
FUNCTION is_in_range, x, lo, hi
   RETURN, (x ge lo) and (x le hi)
END

;----------------------------------------------------------------------------
;  hanning    alternative (flatter than one in ~idl/lib)
;     Hanning function.  This one always square.
;
;From H. Lin's file ccdcal5.pro, 29-Jan-98
FUNCTION hanning_alt, n, m

   k = min([n,m])
   x = fltarr (k) + 1.0
   tenth =  long (k*.2)
   cons = !pi/tenth
   FOR i = 0,tenth DO BEGIN
       x(i) = (1.0 - cos (i*cons))/2.0
       x(k-i-1) = x(i)
   ENDFOR

   RETURN, x # x
END

;----------------------------------------------------------------------------
function apod_image,nx,ny
    ; apodizes image - alternative to others
    ; Based on RR & ADW routine from plotko_power
  
  apocube = fltarr(nx,ny)
  apod=0.1
  
 ; define spatial apodization
  apodx = fltarr(nx)+1
  apody = fltarr(ny)+1
  if (apod ne 0) then begin
    apodrimx=apod*nx
    apodrimy=apod*ny
    apodx[0] = (sin(!pi/2.*findgen(apodrimx)/apodrimx))^2
    apody[0] = (sin(!pi/2.*findgen(apodrimy)/apodrimy))^2
    apodx = apodx*shift(rotate(apodx,2),1)
    apody = apody*shift(rotate(apody,2),1)
    apodxy = apodx # apody
 endif

  ; spatial gradient removal + apodizing per image
  ;xf = fltarr(nx,ny)+1.
  ;yf = xf
  
 ; img = apocube[*,*,it]
 ; avg = avgstd(img)
 ; fitp = mpfit2dfun('gradient',xf,yf,img,fltarr(nx,ny)+1,[1000.,0.,0.],$
 ;  /quiet)
 ; fit = fitp[0]+xf*fitp[1]+yf*fitp[2]
 ; apocube[*,*,it] = (img-fit)*apodxy + avg
  

  ; done
  return,apodxy;apocube
end

;----------------------------------------------------------------------------
; shift_img   shift image by given offsets, using IDL 
;     routine poly_2d with cubic interpolation
;     Works on single image or datacube

FUNCTION shift_img,img,offsets
     sz = size(img)
    IF (sz(0) eq 2) THEN $
      RETURN, poly_2d(img,[-offsets(0),0.,1.,0.],[-offsets(1),1.,0.,0.],cubic=-0.5)
    IF (sz(0) eq 3) THEN FOR i=0,sz(3)-1 DO $
      img(*,*,i) = poly_2d(img(*,*,i),[-offsets(0,i),0.,1.,0.],[-offsets(1,i),1.,0.,0.],cubic=-0.5)
    RETURN,img
END

;----------------------------------------------------------------------------

FUNCTION tr_get_disp_2b, data, shift=shift, mad=mad, debug=debug,nowindow=nowindow,$   
                         themax=themax, quiet=quiet, nomad=nomad

IF NOT keyword_set(mad) THEN mad=0
nmad = mad > 5
nmad = 2*(nmad/2)+1
nmad2 = (nmad-1)/2
dxmad = [-1.,0.,1.]
dymad = dxmad
nomad = 0   ;added by RJM - allows for tracking of No. of failed MADs
;  TDT 11-Mar-02    arrays for 2-D least square fit
m = fltarr(5,3,3)
FOR i=0,2 DO BEGIN
  m(0,*,i) = dxmad 
  m(1,*,i) = dymad(i)
  m(2,*,i) = dxmad^2
  m(3,*,i) = dxmad*dymad(i)
  m(4,*,i) = dymad(i)^2
ENDFOR
  
m = reform(m,5,9)


w = 1. + 0.*fltarr(9)

IF NOT keyword_set(shift) THEN shift=0
IF NOT keyword_set(debug) THEN debug=0
errorstring = 'Minimum MAD not in '+string(nmad,format='(I2)')+'^2 area--image #, xmin, ymin:'

sz = size(data)
nx = sz(1) & ny = sz(2) & nz = sz(3)
disp = fltarr (2,nz)
themax=fltarr(nz)

;#####################################################################
;Sets size of region to carry our cross-correlation on
;#####################################################################

; ADW  20081025  option to use full frame
IF keyword_set(nowindow) THEN BEGIN
  nnx = nx
  nny = ny
ENDIF ELSE BEGIN
  ; TDT  11-Sep-98  added + 1.e-5 to make this work right!
  nnx = 2^long (alog10 (min ([nx, ny]))/.30103 + 1.e-5)
  nny = nnx
  
ENDELSE


; TDT 29-Jan-98  added float to this next statement
nnsqd = float(nnx)*float(nny)
appodize = hanning_alt (nnx, nny)



;#####################################################################
;Selects reference region from data, applies hanning funtion to reduce edge effects
;followed by a fourier transform
;#####################################################################

ref = data ((nx-nnx)/2:(nx+nnx)/2-1, (ny-nny)/2:(ny+nny)/2-1, 0)
tref = conj (fft ((ref-total(ref)/nnsqd)*appodize, -1))

var_ref=(moment(ref))[1] ;variance of ref
;tref_check

;#####################################################################
;Loop compares reference image to other images in array via cross-correlation
;#####################################################################

FOR i = 1, nz-1 DO BEGIN
   scene = data ((nx-nnx)/2:(nx+nnx)/2-1,(ny-nny)/2:(ny+nny)/2-1, i)
   tscene = fft ((scene-total(scene)/nnsqd)*appodize, -1)
   cc = shift (abs (fft (tref*tscene, 1)), nnx/2, nny/2)
   printerror = 1

   mx = max (cc, loc)   ; locate peak of Cross Correlation

   var_scene=(moment(scene))[1] ;variance of scene
   tscene_check = fft (scene, -1)
   
   ; normalised cc coefficient? What about hanning?
   themax[i]=mx/sqrt(var_ref*var_scene) 
   ;Alternative?
   ;themax[i]=max(abs (fft (tref_check*tscene_check, 1))) 

   xmax0 = loc mod nnx
   ymax0 = loc/nnx
   xmax = ( (xmax0 > nmad2) < (nnx-nmad2-1) )
   ymax = ( (ymax0 > nmad2) < (nny-nmad2-1) )
   
   IF debug THEN BEGIN 
      print,'Fourier Cross-correlation Peak: ',xmax0,ymax0
      print,cc(xmax-2:xmax+2,ymax-2:ymax+2), format='(5F8.1)'
   ENDIF
   
   ;#####################################################################  
   ;Creates an array (mad,mad) around maximum of Cross-correlation 
   ;##################################################################### 
   cc = -cc(xmax-nmad2:xmax+nmad2,ymax-nmad2:ymax+nmad2)
   

   ;#####################################################################   
   ; Mean Absolute Difference algorithm centered on xmax & ymax 
   ;##################################################################### 
   IF (mad) THEN BEGIN
  
      dx = nnx/2-xmax
      dy = nny/2-ymax
      ;  TDT  11-Mar-02  bug here -- shouldn't divide this by 2
      ; nnx2 = (nn/2-abs(dx)-nmad2-1)/2
      nnx2 = (nnx/2-abs(dx)-nmad2-1)
      nxl = nnx/2-nnx2
      nxh = nnx/2+nnx2
      ; nny2 = (nn/2-abs(dy)-nmad2-1)/2
      nny2 = (nny/2-abs(dy)-nmad2-1)
      nyl = nny/2-nny2
      nyh = nny/2+nny2
      area = float(nxh-nxl+1)*float(nyh-nyl+1)
 
      ;if statement added by RJM 12-06-12 - stops the routine crashing - if MAD cannot be
      ; performed CC'd array will be used even if MAD is set. Need to add feature that 
      ; distinguishes between shifts calculated by MAD or CC
      IF (nxh gt nxl and nyh gt nyl) THEN BEGIN
         cc = fltarr(nmad,nmad) ; only create new cc array if MAD can be performed
  
  
         FOR idx=-nmad2,nmad2 DO BEGIN
             FOR idy=-nmad2,nmad2 DO BEGIN

             cc(idx+nmad2,idy+nmad2)=total(appodize(nxl:nxh,nyl:nyh)*abs(ref(nxl:nxh,nyl:nyh) - $
                scene(nxl-dx+idx:nxh-dx+idx,nyl-dy+idy:nyh-dy+idy)))/area

             ENDFOR
         ENDFOR
      ENDIF
  
      cc = cc^2
      
      IF debug THEN BEGIN
         print,'Squared MAD array:'
         print,cc, format='('+string(nmad,format='(i2)')+'F8.1)'
      ENDIF

   ENDIF
     
  
   ;#####################################################################     
   ; Locate minimum of MAD^2 or -Cross-correlation function
   ;   hope nmad x nmad is big enough to include minimum
   ;#####################################################################

   mx = min (cc, loc)
  
   xmax7 = loc mod nmad
   ymax7 = loc/nmad
   IF (xmax7 gt 0 and xmax7 lt (nmad-1) and ymax7 gt 0 and ymax7 lt (nmad-1)) THEN BEGIN

       ;TDT  11-March-2002    2-D least square fit for minimum
       IF (mx gt 0) THEN w = replicate(1./float(mx),9) else w=replicate(1.,9) ;Creates fltarr[9] with components equal to min
       ccmin = cc(xmax7-1:xmax7+1,ymax7-1:ymax7+1) ;always 3 by 3 array
    
    
       ;/////////////////////////
       IF (debug) THEN $
       ;a=regress(X,Y,weights,yfit,const,sigma,ftest.....)
       a = regress(m,reform(ccmin(*,*),9),w,dum,a0,sigma,Ftest,R,Rmul,Chisq) else $
       a = regress(m,reform(ccmin(*,*),9),w,a0)  
       ;reform turns ccmin to fltarr[9] and is the dependent variable
       ; m is fltarr(5,9) and is the independent variable
       ; W is the weights (same weight for each value!)
       ; a is coefficients
       ;The fit is a quadratic polynomial regression equation of the form Z=a0+a(0)*x+a(1)*y+a(2)*x^2+a(3)*xy+a(4)*y^2
       ;i.e. fits a 3-D parabola to the local area around the CC maximum  (see, e.g. Jeffrey Edwards 2007 -Polynomial Regression and Response Surface Methodology)    
                                  
       ;/////////////////////////


       ;The location of the stationary points of the fitted 3-D parabola are given by the following
       det = abs(4.*a(2)*a(4)-a(3)*a(3))
       xfra = -(2.*a(0)*a(4)-a(1)*a(3))/det
       yfra = -(2.*a(2)*a(1)-a(0)*a(3))/det 
       
   ENDIF ELSE BEGIN 
       xfra = 0.
       yfra = 0.
    
       IF NOT keyword_set(quiet) THEN BEGIN ;added by RJM to suppress text
           if (printerror) then print,errorstring,i,xmax7-nmad2,ymax7-nmad2
       ENDIF
       IF (printerror) THEN nomad=nomad+1
       printerror=0
   ENDELSE




   xfra = xfra + xmax7 - nmad2
   yfra = yfra + ymax7 - nmad2
   IF debug THEN print,xfra,yfra,format='("Fractional dx, dy: ",2F10.3)'
   IF debug THEN print,det,Chisq,format='("Det, Chisq: ",2E12.4)'
   xmax = xfra + xmax 
   ymax = yfra + ymax 



   disp(0,i) = (nnx/2-xmax)
   disp(1,i) = (nny/2-ymax)
   
   if debug then print, i, disp(0,i), disp(1,i), $
     format='("Image ",I4, "    Final offsets ",2F10.2,/)'
   if (shift) then data(*,*,i) = shift_img(data(*,*,i),disp(*,i))

ENDFOR
return, disp
END
