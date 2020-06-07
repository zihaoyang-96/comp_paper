;PURPOSE: Calculates the angle of wave propagation from a data cube.
;         Current version is designed for CoMP velocity data.
;
;
;
;INPUTS: spec - cube of FFT'd data that runs from f=nt/2 to f=0 (i.e., half cube)
;
;OPTIONAL INPUTS: /save_coh - saves coherence values
;        
;OUTPUTS: wave_angle - array of wave angles
;         angle_error - error on wave angles
;         coh_measure - measure of coherence, only created if /save_coh is set
;
;
;NOTES: Array defined call mask. If STEREO idl files are installed then routine can get confused
;       and tries calling the STEREO routine. Maybe should change the name of mask array RJM 06/2015
;
;HISTORY:- Created by S Tomczyk
;        - Modularised by R J Morton 2015
;        - Created index to keep together useful quantities from wave tracking RJM 04/2016
;        - Edited as part of Python port - RJM 03/2019
;        - v0.3 (RJM 05/2019) Introducing versioning to keep track of code development. Given past major changes, 
;          the first numbered version is as version 0.3.
;
;TO DO: Is the current method of calculating coherence correct? Should summing be used instead of smoothing? - RJM
;


PRO wave_angle_calc,spec,index,config,wave_angle=wave_angle,$
                    angle_error=angle_error,coh_measure=coh_measure,save_coh=save_coh,debug=debug


nx=index.NAXIS1
ny=index.NAXIS2
nt=index.NFRAMES

;Unpack dictionaries here - slow to access so avoid in loops
mask=index.mask
coherenceBoxLength=2.*(config['coherence'].boxHalfLength)+1
coherenceSmoothing=(config['coherence'].smoothing)
coherenceBoxHalfLength=(config['coherence'].boxHalfLength)
minNumberCoherentPixels=(config['coherence'].minNumberPixels)
coherenceLimit =(config['coherence'].limit)



wave_angle=fltarr(nx,ny)
angle_error=fltarr(nx,ny)
angle_error2=fltarr(nx,ny) ; For testing RJM
coh_measure = fltarr(nx,ny,2)

xbox=rebin(findgen(coherenceBoxLength)-coherenceBoxHalfLength,coherenceBoxLength,coherenceBoxLength)  ;coordinates of points in box
ybox=transpose(xbox)

;Set up filter
nspec=nt/2
freq   = findgen(nspec)/(float(nspec*2)*index.normalCadence)
filter = exp( -(freq-(config['filter'].centralFrequency) )^2/(config['filter'].width)^2 )
filter(0) = 0.                      ; set dc to zero
filter(where (freq lt .001)) = 0.   ; set low frequencies to zero
filter=filter/total(filter)


; open windows
if keyword_set(debug) then begin
  loadct, 39, /silent
  f=8
  window,1,xs=coherenceBoxLength*f,ys=coherenceBoxLength*f,xpos=0,ypos=0
  window,4,xs=nx,ys=ny,retain=2,xpos=400+nx+20,ypos=0, title=date
  window,5,xs=ny,ys=ny,retain=2,xpos=400,ypos=ny+50
  window,6,xs=512,ys=512,xpos=400,ypos=800
  device,decomposed=0
endif

if not keyword_set(debug) then begin
 pixcounter = 0L
 old_perc_proc = 0
 tot_pix = long(n_elements(where(index.mask eq 1)))
endif

conjspec=conj(spec)

;Although this is quicker - there is some funny behaviour going on with smooth
;can't work out what!
;lar_g=smooth(real_part(spec*conjspec),[1,1,coherenceSmoothing])
lar_g=fltarr(nx,ny,nspec)
FOR iy=0,ny-1 DO FOR ix=0,nx-1 DO $
    IF (mask[ix,iy] EQ 1) THEN lar_g[ix,iy,*]=smooth(real_part(reform(spec[ix,iy,*]*conjspec[ix,iy,*])),coherenceSmoothing)

;calculates rms value of spectrum
;Used to provide secondary condition on for loops
;and exclude pixels with low signal
rms_spec=total(abs(spec),3)

message, /cont,'Starting processing...'
FOR iy=0,ny-1 DO FOR ix=0,nx-1 DO $
    IF (mask[ix,iy] EQ 1) and (rms_spec[ix,iy] GT 1e-7) THEN BEGIN
  
  IF NOT keyword_set(debug) THEN BEGIN
     perc_proc = round((float(pixcounter)/float(tot_pix))*100.)
     pixcounter = pixcounter+1
     IF perc_proc GT old_perc_proc THEN BEGIN
            IF strmid(getenv('TERM'),0,5) EQ 'xterm' THEN BEGIN
     		progress_comp, perc_proc, 100, msg='Computing wave propagation angles'
    	    ENDIF ELSE message, /cont, strcompress(string(perc_proc,format='(I3)'),/rem)+'% processed.'
    	    old_perc_proc = perc_proc
     ENDIF
  ENDIF
 
  ;================================================================
  ;  compute coherence for contiguous pixels with coherence greater than limit
  ;================================================================

  coh=fltarr(coherenceBoxLength,coherenceBoxLength)           ;initialize coherence array to zero
  coh_mask=fltarr(coherenceBoxLength,coherenceBoxLength)      ;initialize mask of good points in coherence
  coh(coherenceBoxLength/2,coherenceBoxLength/2)=1.
  coh_mask(coherenceBoxLength/2,coherenceBoxLength/2)=1.

  spec1=reform(spec[ix,iy,*])
  g1=reform(lar_g[ix,iy,*])

  count=1
  ncoh=1
  WHILE count GT 0 AND ncoh LT 100 DO BEGIN
        cont=where(coh_mask eq 1.)
        cnew=[cont-1,cont+1,cont+coherenceBoxLength,cont-coherenceBoxLength]  ;test pixels adjacent to contiguous pixels

        new_mask=fltarr(coherenceBoxLength,coherenceBoxLength)   ;identify new contiguous pixels
        new_mask(cnew)=1.
        new_mask=new_mask*(1.-coh_mask)

        new=where(new_mask gt 0.,new_count)
       
        count=0
        FOR ic=0,new_count-1 DO BEGIN
            icx=new(ic) mod coherenceBoxLength
            icy=new(ic)/coherenceBoxLength

            ii=ix-coherenceBoxLength/2+icx
            jj=iy-coherenceBoxLength/2+icy

            ;Exceptions here prevent out of bounds values and also 
            ;performing correlation with pixels with no signal
            IF (ii LT 0) OR (ii GT nx-1) OR (jj LT 0) OR (jj GT ny-1) OR (rms_spec[ii,jj] EQ 0) or mask[ii,jj] eq 0 $ 
                 THEN coh(icx,icy)=0. ELSE BEGIN
            	       IF coh(icx,icy) EQ 0. THEN BEGIN
              	         
                         spec2c=reform(conjspec(ii,jj,*))
                         cspec=spec1*spec2c
                         cspec=smooth(cspec,coherenceSmoothing)
                         g2=reform(lar_g[ii,jj,*])
                         
                         ;Don't do calculation with DC component - it is zero and leads to maths errors                                         
                         coh(icx,icy)=mask[ii,jj]*total( filter[1:-1]*( abs(cspec[1:-1])/sqrt(g1[1:-1]*g2[1:-1]) ) )
                         ;see McIntosh et al. 2008. Sol. Phys.
                      ENDIF

                      IF coh(icx,icy) GT coherenceLimit THEN BEGIN
                  	     coh_mask(icx,icy)=1.
                  	     count=count+1
                      ENDIF
            ENDELSE
         ENDFOR
         
         good=where(coh_mask eq 1.,ncoh) 

   ENDWHILE

   IF keyword_set(debug) THEN BEGIN
    wset,1
    tv,bytscl(rebin(coh*coh_mask,f*coherenceBoxLength,f*coherenceBoxLength,/sample),0.,1.)
   ENDIF   

   ;IF (keyword_set(save_coh) and (ncoh gt minNumberCoherentPixels)) THEN BEGIN
   IF (ncoh gt minNumberCoherentPixels) THEN BEGIN
     island = coh*coh_mask
     ellipsepts = fit_ellipse(where(island ne 0.),axes=axes,xsize=coherenceBoxLength,ysize=coherenceBoxLength)
     coh_measure[ix,iy,*] = axes
     if ((size(ellipsepts))[0] gt 0) and keyword_set(debug) then plots, f*ellipsepts, color=254, /device
   ENDIF

   ;================================================================ 
   ;  find angle of coherence island if enough good points
   ;================================================================

   IF ncoh GT minNumberCoherentPixels THEN BEGIN

 
       weight=coh(good)^2
       theta=min_dist_fit(xbox(good),ybox(good),error=error) ; analytical solution to the minimum distance problem
       ;theta=mle_angle(xbox(good),ybox(good)) 
       ;error=0.

       IF keyword_set(altangle) THEN BEGIN

           xe=xbox(good) & ye=ybox(good)
           sxy=total((xe-mean(xe))*(ye-mean(ye)))

           IF round(sxy) NE 0. THEN BEGIN
    
                 sixlin,xbox(good),ybox(good),a,siga,b,sigb
                 theta=[atan(b[2])/(!dtor)]
                 error=[atan(sigb[2])/(!dtor)]
           ENDIF
       ENDIF
        
       wave_angle(ix,iy)=theta
       angle_error(ix,iy)=error


       
       IF keyword_set(debug) THEN BEGIN
          xfit=findgen(2*coherenceBoxHalfLength+1)-coherenceBoxHalfLength
          coef=[0.,tan(theta*!dtor)]
          yfit=poly(xfit,coef)  ;C0 + C1*x + C2*x^2

          wset,6
          plot,xbox(good),ybox(good),psym=6,xr=[-coherenceBoxHalfLength,coherenceBoxHalfLength],yr=[-coherenceBoxHalfLength,coherenceBoxHalfLength]
          oplot,xfit,yfit
          wset,4
          tv,bytscl(wave_angle,-90,90)
          wset,5
          tv,bytscl(angle_error,0,10)
          wait,0.005
          print,theta,error
       ENDIF

  ENDIF
  
ENDIF


END
