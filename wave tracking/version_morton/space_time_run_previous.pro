pro space_time_run, date, velocity, angle, index,config, name_addon, $
	                   debug=debug, wr_speeds=wr_speeds
	                   

;PURPOSE:  Procedure to compute space-time diagram of velocity time series over a whole image.
;
;
;OUTLINE:  The track is computed from the angle map. Velocities are then interpolated onto a time-distance
;          map. The mean value of the velocity for the map is subtracted. The temporal mean for each position
;          along the track is also subtracted before FFT. FFT of map taken and high spatial and temporal
;          frequencies are removed. Then separates prograde and retrograde components (i.e., where w=ck & w=-ck).
;          Inverse fourier transform is taken for velocities which are then sent to compute_speed to compute the
;          phase speed. 
;
;
;INPUTS: data- date of observations
;        velocity - Aligned + trimmed Doppler velocity cube
;        angle - array of angles calculated form coherence islands (see wave_tracking.pro + McIntosh et al. 2008)
;        config - general useful variables
;        name_addon - additional quantifier to be added to the name of the .sav filter 
;
;OPTIONAL INPUTS: /debug - Use to debug code + process - opens multiple windows + prints to screen
;                 /wr_speed - Saves pro/retro - grade velocity maps
;                 npt_init - length of track to use for cross-correlation, maximum useful length appears to 
;                            be around 21. Signal only correlates well over sort distances. Specify in comp_config.pro file
;                
;
;CALLS: compute_speed.pro 
;
;
;HISTORY: Created S Tomczyk
;        - Fixed minor bug for dealing with tracks that are shorter than specified track length.
;          Turned continue's to break's - if track hits boundary it stops R Morton 11/2014 
;        - Added if statements to debug only sections of code - hopefully reduce run time.
;          Edited out power calculation to reduce run time  
;          Attempted to remove bad series from calculation by looking at RMS values
;          - tiny values of RMS suggest incomplete series.
;          Changed start pixel selection criteria based on angle value to mask value 
;          Added break statement if track hits bad mask value
;          Applied Median smoothing of wave angle map to reduce noise  R Morton 03/2016
;        - v0.3 (RJM 05/2019) Introducing versioning to keep track of code development. Given past major changes, 
;          the first numbered version is as version 0.3.
;
;TO DO: Is there a smarter way to interpolate velocity onto track - at the moment there is a for loop over time.
;       * Fourier filtering causes artefacts - 'curled up/down' ends due to edge effects - increases gradient!
;       									 -  max propagation speed of ~1347 km/s (highest frequency * track length). 

ans=' '

;Define constants & variables from index
nx=index.NAXIS1          ; x dimension of velocity data cube
ny=index.NAXIS2          ; y dimension
ntime=index.NFRAMES         ; time dimension

normalCadence=index.normalCadence
xscale=index.xscale

lowerMaskRadius=index.lowerMaskRadius+config['maxOcculterOffset']
upperMaskRadius=index.upperMaskRadius
mask=index.mask


;number of time points to use (make even)
if ((ntime mod 2) eq 0) then nt = ntime else nt = ntime-1
velocity=temporary(velocity[0:nx-1,0:ny-1,0:nt-1])

;;; *** choose length of track here ***
npt_init=config['maxTrackLength'] ;number of steps to map along track (make odd)
dx=1.0      ;step size along track

;===========================
npt= npt_init
ang_track=fltarr(npt)   ;angle along track
vmap=fltarr(nt,npt)     ;array to hold space-time diagram

f=findgen(nt)/(float(nt)*normalCadence)       ;compute temporal frequency scale
for i=0,nt/2 do f(nt-i-1)=f(i+1)
freq=rebin(f,nt,npt)                         ;along track,

sf=findgen(npt)/(float(npt)*xscale)          ;compute spatial frequency scale ;length is in Mm (res=4.46*0.725 Mm)
for i=0,npt/2 do sf(npt-i-1)=sf(i+1)         ;make retrograde frequencies negative
sp_freq=rebin(transpose(sf),nt,npt)          ;along track

phase_speed=freq/sp_freq                     ;compute phase speed = (w/k)

;  identify acceptable phase speed regions
;good_phase=where(freq gt 0.0025 and freq lt 0.010,complement=bad_phase)

w=0.001      ; create filter to remove low frequency signal;
             ; Note: This is not a gaussian filter
             ; like in the other filtering case.
             ; Here it just removes the low frequency
filter=1.-exp(-freq^2/w^2)
;filter2=1.-exp(-(nt/2.-findgen(nt))^2/10.^2)
;filter2=1.rebin(filter2,nt,npt)

angle=median(angle,3) ;Median smoothing of wave angle map to reduce noise

;##############################################################
; Define x and y limits for map
; entire map
;##############################################################

;Already done in wave_tracking_v2
x=rebin(findgen(ny)-(nx/2.)+0.5,nx,ny)
y=transpose(x)
rad=sqrt(x^2+y^2)

;mask[where(rad lt lowerMaskRadius)]=0.
;mask[where(rad gt upperMaskRadius)]=0.

;Define good mask point based on coherence measure from wave angle calculation




;THIS IS UNDER TEST
;Aiming to reject poor data values instead of processing them
;Plotting RMS velocity for non zero pixels reveals bimodal
;distribution
;Those in first distribution have tiny RMS values, so we should rule them out as being poor series
;The following implements a very basic attempt to separate into two distinct values
;May not work for all data sets - a better method is required. 
;
;summ=sqrt(total(velocity^2,3))/nt
;plothist,summ[where(summ gt 0)]<0.1,xhist,yhist,/noplot,bin=1e-3
;r=min(yhist,loc)
;cutoff=xhist[loc]
;mask[where(summ lt cutoff)]=0.

;good=where(mask eq 1, complement=bad)
;angle[bad]=0. ; Edited out by RJM, use mask values rather than angles
               ;Some zero angles may be real!!    



;new limits should help speed up routine (marginally) by only looping over valid elements (RJM)
xstart=310-upperMaskRadius ;0 
xend=310+upperMaskRadius ;nx-1
ystart=310-upperMaskRadius;0
yend=310+upperMaskRadius;ny-1

;IF NOT keyword_set(debug) THEN BEGIN
 pixcounter = 0.
 old_perc_proc = 0
 tot_pix = (xend-xstart)*1L*(yend-ystart) ;n_elements(where(mask ne 0))
;ENDIF


IF keyword_set(debug) THEN BEGIN
	 window,0,xs=nx,ys=ny,retain=2,xpos=0,ypos=700
	 window,1,xs=nx,ys=ny,retain=2,xpos=nx+10,ypos=700
	 ;window,2,xs=nt*2,ys=npt*2,title='Space-Time Diagram'
	 window,3,xs=nt*2,ys=npt*2,xpos=0,ypos=2*npt+45,title='Prograde Filtered'
	 window,4,xs=nt*2,ys=npt*2,xpos=0,ypos=4*npt+90,title='Retrograde Filtered'
	 window,8,xs=900,ys=600

	 device,decomposed=0

	 wset,0    ;display images
	 img1=velocity[*,*,0]
	 blue_red
	 tvlct, brr,brg,brb, /get
	 brr[255] = 0 & brg[255] = 0 & brb[255] = 0
	 tvlct, brr,brg,brb
	 img1 = bytscl(img1,-8,8,top=254)
	 tv, img1
	 wset,1
	 loadct, 13, /silent
	 tvlct, r,g,b, /get
	 r[255] = 255 & g[255] = 255 & b[255]= 255
	 tvlct, r,g,b
	 img2=bytscl(angle,-90,90,top=254)
	 tv,img2
ENDIF


;##############################################################
;Define arrays to store prograde and retrograde quantities 
;##############################################################

pro_power=fltarr(nx,ny)
ret_power=fltarr(nx,ny)
sign_power=fltarr(nx,ny)
pro_speed=fltarr(nx,ny)
ret_speed=fltarr(nx,ny)
phase_speed=fltarr(nx,ny)
pro_speed_err=fltarr(nx,ny)
ret_speed_err=fltarr(nx,ny)
phase_speed_err=fltarr(nx,ny)
pro_no_tr=fltarr(nx,ny)
ret_no_tr=fltarr(nx,ny)


IF wr_speeds EQ 1 THEN BEGIN
	 p_vel=fltarr(nx,ny,nt,npt)
	 r_vel=fltarr(nx,ny,nt,npt)
ENDIF




;##############################################################
;Main section of the routine
;##############################################################

message, /cont,'Starting processing...'
old_perc_proc=0
FOR iy=ystart,yend DO BEGIN
  FOR ix=xstart,xend DO BEGIN
      
      IF keyword_set(debug) THEN BEGIN
         perc_proc = round((float(pixcounter)/float(tot_pix))*100.)
         pixcounter = pixcounter+1.
         if perc_proc gt old_perc_proc then begin
            if strmid(getenv('TERM'),0,5) eq 'xterm' then begin
               progress_comp, perc_proc, 100, msg='Computing phase speeds'
            endif else begin
               message, /cont, strcompress(string(perc_proc,format='(I3)'),/rem)+'% processed.'
            endelse
            old_perc_proc = perc_proc
         ENDIF
      ENDIF

      IF mask[ix,iy] NE 0. THEN BEGIN ;changed from angle
         
         ;Here the condition of the each hemisphere is taken which will be used later.
         IF ix le nx/2-1 THEN hemisphere=1 ELSE hemisphere=-1

         IF keyword_set(debug) THEN BEGIN
            img1 = velocity[0:nx-1,0:ny-1,0]   ;reinitialize velocity image
            img1 = bytscl(img1,-8,8,top=254)
            img2 = bytscl(angle,-90,90,top=254)
         ENDIF

         npt=npt_init         ;every times we have to do this
         xtrack=fltarr(npt)   ;x and y position along track
         ytrack=fltarr(npt)   ;............................
         ang_track=fltarr(npt)   ;angle along track

         imid=fix(npt/2)

         ;##############################################################
	  	 ; use angle at each point to define track, starting at cursor position as central point
	  	 ;
	  	 ;Place initial position +angle at centre of arrays xtrack + ytrack +ang_track
	  	 ;##############################################################
	  
	  	 xtrack(imid)=float(ix) & ytrack(imid)=float(iy)
	  	 ang_track(imid)=angle(ix,iy)

	  	 ntrack=0        ;this index will count the number of pixels on the track

		 ;this condition gives the lower limit of the track
	  	 IF sqrt( (float(ix)-nx/2.+0.5)^2 + (float(iy)-ny/2.+0.5)^2 )-imid LT lowerMaskRadius THEN $
	             short_track=1 else short_track=-1 
  
	  	 ;##############################################################
	  	 ;  first, move out from cursor position
	  	 ;##############################################################
	  
	  	 FOR i=imid+1,npt-1 DO BEGIN

	    	     ;Calculates which x + y direction to follow
	    	     ;Uses angle calculated from cross-correlation in wave_tracking
	    	     ;Use of hemisphere ensure signal always steps away from occulting disk
	 	     ;This means value should always greater than 0 and less than nx-1 (ny-1)
	 	     xtrack(i)=0. > xtrack(i-1)-hemisphere*dx*cos( ang_track(i-1)*!dtor) < float(nx-1)
	 	     ytrack(i)=0. > ytrack(i-1)-hemisphere*dx*sin( ang_track(i-1)*!dtor) < float(ny-1)

		     intx=0 > fix(xtrack(i)+0.5) < (nx-1)
		     inty=0 > fix(ytrack(i)+0.5) < (ny-1)
	    

	             ;this means if you are near the lower boundary
	             pixcheck=sqrt((intx-nx/2.+0.5)^2 + (inty-ny/2.+0.5)^2)
		     IF pixcheck lt lowerMaskRadius THEN BEGIN
	       	        ;this will not allow to track further down in the inward direction!
		        xtrack(i)=0.
	 	        ytrack(i)=0.

		        BREAK
	                
		     ENDIF
	   
		     IF pixcheck GT upperMaskRadius THEN BEGIN
	 	        xtrack(i)=0.
		        ytrack(i)=0.

		        BREAK
		     ENDIF

	             ;IF mask[intx,inty] LT 1 THEN BEGIN ;don't use bad pixels
	 	     ;   xtrack(i)=0.
		     ;   ytrack(i)=0.
		     ;   BREAK
		     ;ENDIF

		     ntrack=ntrack+1

	         IF keyword_set(debug) THEN BEGIN
			 	img1(intx,inty)=255 ;put track into displayed images
			 	img2(intx,inty)=255
		     ENDIF
		     
	   	     ang_track(i)=angle(intx,inty)  ; nearest neighbor interp for angle
	                                   ; (works better than bilinear at boundaries)

	   	     ; if angle is zero (not defined), use previous value
		     ;IF ang_track(i) eq 0. THEN ang_track(i)=ang_track(i-1)

		     ; if big difference in adjacent angles, resolve 180 degree ambiguity
		     ; by choosing angle closest to ang_track(i-1)
		     IF abs(ang_track(i)-ang_track(i-1)) gt 90. THEN BEGIN
		       IF ang_track(i)-ang_track(i-1) GT 0. THEN ang_track(i)=ang_track(i)-180. ELSE ang_track(i)=ang_track(i)+180.
		     ENDIF

         ENDFOR

		 ;##############################################################
		 ;  next, move in from cursor position
		 ;##############################################################
		 FOR i=imid-1,0,-1 DO BEGIN
		    factor=(ang_track(i+1)+180.)*!dtor
		    xtrack(i)=0. > xtrack(i+1)-hemisphere*dx*cos(factor ) < float(nx-1)
		    ytrack(i)=0. > ytrack(i+1)-hemisphere*dx*sin( factor) < float(ny-1)

		    intx=0 > fix(xtrack(i)+0.5) < (nx-1)
		    inty=0 > fix(ytrack(i)+0.5) < (ny-1)

		    ;this means if you are near the lower boundary
	            pixcheck=sqrt((intx-nx/2.+0.5)^2 + (inty-ny/2.+0.5)^2)
		    IF pixcheck lt lowerMaskRadius THEN BEGIN
		       ;this will not allow to track further down in the inward direction!
		      xtrack(i)=0.
		      ytrack(i)=0.
			  BREAK ; stops tracks when hits lower boundary
	        ENDIF
	   
		    IF pixcheck GT upperMaskRadius THEN BEGIN
		        xtrack(i)=0.
		        ytrack(i)=0.
			    BREAK ; stops tracks when hits upper boundary
    		ENDIF

	           ;  IF mask[intx,inty] LT 1 THEN BEGIN ;don't use bad pixels
	 	   ;     xtrack(i)=0.
		   ;     ytrack(i)=0.
		   ;     BREAK
		   ;  ENDIF

		    ntrack=ntrack+1

	        IF keyword_set(debug) THEN BEGIN
		       img1(intx,inty)=255 ;put track into displayed images
		       img2(intx,inty)=255 
	        ENDIF

		    ang_track(i)=angle(intx,inty)  ; nearest neighbor interp for angle
	                                   ; (works better than bilinear at boundaries)

		    ; if angle is zero (not defined), use previous value
		    ;IF ang_track(i) EQ 0. THEN ang_track(i)=ang_track(i+1)

		    ; if big difference in adjacent angles, resolve 180 degree ambiguity
	            ; by choosing angle closest to ang_track(i-1)
		    IF abs(ang_track(i)-ang_track(i+1)) GT 90. THEN BEGIN
		      IF ang_track(i)-ang_track(i+1) GT 0. THEN ang_track(i)=ang_track(i)-180. ELSE ang_track(i)=ang_track(i)+180.
		    ENDIF

		 ENDFOR


         IF ntrack LT 4 THEN CONTINUE ;this tells that we are taking at least 5 points in the track
         
         IF ntrack NE ntrack/2*2+1 THEN ntrack=ntrack/2*2+1  ;Gives the total number of pixels in the track.
                                                             ;For normal case it is npt, but when the track
                                                             ;is shorter, it is < npt
	                				       ;This also makes npt odd
      	 IF ntrack LT npt THEN BEGIN    ;......If the track is shorter, then it will do the followings
		       in=where(xtrack ne 0)
		       newxtrack=xtrack(in) ;......it creates new track of shorter length
		       xtrack=newxtrack
		       newytrack=ytrack(in)
		       ytrack=newytrack
		       npt=n_elements(in)
		       ang_track=ang_track(in)
	 	 ENDIF ELSE BEGIN
	       	   npt=npt_init                    ; We need the following lines for the new track
	 	 ENDELSE

	 	 vmap=fltarr(nt,npt)           ;Again we need to make it (array to hold space-time diagram)

		 ; interpolate velocity onto track
		 FOR i=0,nt-1 DO vmap[i,0:npt-1]=interpolate(velocity[0:nx-1,0:ny-1,i],xtrack,ytrack,cubic=-0.5)
	         

	  
		 IF keyword_set(debug) THEN BEGIN
		    wset,0                      ;re-display velocity image with track
		    blue_red
		    tvlct, brr,brg,brb, /get
		    brr[255] = 0 & brg[255] = 0 & brb[255] = 0
		    tvlct, brr,brg,brb
		    tv,img1
		    wset,1                      ;re-display angle image with track
		    loadct, 13, /silent
		    tvlct, r,g,b, /get
		    r[255] = 255 & g[255] = 255 & b[255]= 255
		    tvlct, r,g,b
		    tv,img2
		 ENDIF
	  
		 ;##############################################################
		 ; Filter velocity map & calculate pro/retro - grade velocity & power
		 ;##############################################################
		  vmap=vmap-mean(vmap)
	  	  vmap=vmap-transpose(rebin(mean(vmap,dim=1),npt,nt),[1,0]) ;remove temporal mean
	 

	  	 trans=fft(vmap,-1)       ;compute fourier transform
		 mag=abs(trans)         ;compute magnitude and phase
		 phase=atan(imaginary(trans),float(trans))
		 back=rebin( rebin(mag(nt/2-(nt/2-2):nt/2+(nt/2-2),*),1,npt) ,nt,npt)
	   	 mag=mag-back       ;remove high temporal frequency noise
	 
		 back=rebin(mag(*,npt/2),nt,npt)
	     mag=mag-back       ;remove high spatial frequency noise
	 	 mag=mag > 0.   ;insure that magnitude is positive
	 
	   	 pw_trans=complex(mag*cos(phase),mag*sin(phase))   ;recompute transform
	     pw_trans=pw_trans*filter      ;remove low temporal frequencies;  filter is defined at top

		 pro_trans=trans             ;select prograde waves (assume nt even, npt odd)
	  	 pro_trans(1:nt/2-1,0:npt/2)=0.
		 pro_trans(nt/2+1:nt-1,npt/2+1:npt-1)=0.
		 pro_vel=float(fft(pro_trans,1))

		 pro_pw_trans=pw_trans             ; separate calculation for the wave power
		 pro_pw_trans(1:nt/2-1,0:npt/2)=0.
		 pro_pw_trans(nt/2+1:nt-1,npt/2+1:npt-1)=0.
		 ;pro_pw_trans(bad_phase)=0.

		 ret_trans=trans             ;select retrograde waves
	  	 ret_trans(1:nt/2-1,npt/2+1:npt-1)=0.
	  	 ret_trans(nt/2+1:nt-1,0:npt/2)=0.
	  	 ret_vel=float(fft(ret_trans,1))

	  	 ret_pw_trans=pw_trans             ; separate calculation for the wave power
	  	 ret_pw_trans(1:nt/2-1,npt/2+1:npt-1)=0.
	  	 ret_pw_trans(nt/2+1:nt-1,0:npt/2)=0.
	  	 ;ret_pw_trans(bad_phase)=0.

		 ; display space time diagram and prograde/retrograde wave components
		 IF keyword_set(debug) THEN BEGIN
		    loadct, 70, /silent
		    wset,3
		    tvscl,bytscl(rebin(pro_vel,nt*2,npt*2,/sample),-2,2)
		    wset,4
		    tvscl,bytscl(rebin(ret_vel,nt*2,npt*2,/sample),-2,2)
		    loadct, 39, /silent
		    wset,8
		    !p.multi=[0,2,2]
		 ENDIF

		 IF wr_speeds eq 1 THEN BEGIN
		    tmp_siz = (size(pro_vel))[2]
		    p_vel[ix,iy,*,0:tmp_siz-1] = pro_vel
		    r_vel[ix,iy,*,0:tmp_siz-1] = ret_vel
		 ENDIF


		 ;##############################################################
	 	 ; compute phase speeds of prograde, retrograde and combination
		 ;##############################################################
	          
	          
		  ;r=compute_speed_new(pro_vel,xscale,normalCadence,no_tr=no_tr)
		  r=compute_speed(pro_vel,xscale,normalCadence)
		  pro_speed(ix,iy)=r(0)
		  pro_speed_err(ix,iy)=r(1)
		  ;pro_no_tr[ix,iy]=no_tr
	    
		  r=compute_speed(ret_vel,xscale,normalCadence,/ret)
		  ;r=compute_speed_new(ret_vel,xscale,normalCadence,/ret,no_tr=no_tr)
		  ret_speed(ix,iy)=r(0)
		  ret_speed_err(ix,iy)=r(1)
	      ;ret_no_tr[ix,iy]=no_tr

		  phase_speed(ix,iy)=(pro_speed(ix,iy)/pro_speed_err(ix,iy)^2 + abs(ret_speed(ix,iy))/$
	                      ret_speed_err(ix,iy)^2)/(1.0/pro_speed_err(ix,iy)^2 + 1.0/ret_speed_err(ix,iy)^2)
		  phase_speed_err(ix,iy)=(1.0/pro_speed_err(ix,iy)^2 + 1.0/ret_speed_err(ix,iy)^2)^(-0.5)

		 ;##############################################################
		 ; compute prograde and retrograde power
		 ;##############################################################
		  pro_pow=abs(pro_pw_trans)^2
		  ;pro_power(ix,iy)=mean(pro_pow(good_phase))
		  pro_power(ix,iy)=total(pro_pow)
		  ret_pow=abs(ret_pw_trans)^2
		  ;ret_power(ix,iy)=mean(ret_pow(good_phase))
		  ret_power(ix,iy)=total(ret_pow)

		 ;##############################################################
		 ; determine whether track is outward or inward, which determines sign of prograde and retrograde
		 ;##############################################################

		  sign_power(ix,iy)=-(rad(xtrack(0),ytrack(0))-rad(xtrack(npt-1),ytrack(npt-1)))/ $
	                      abs(rad(xtrack(0),ytrack(0))-rad(xtrack(npt-1),ytrack(npt-1)))

		  if keyword_set(debug) then !p.multi=0
	          
	  ENDIF
      counter,ix-xstart+(xend-xstart)*1L*(iy-ystart),tot_pix,/percent
   ENDFOR
ENDFOR

;##############################################################
;Save outputs
;##############################################################

nwlst = strcompress(string(config['numberWaveLengths']),/rem)
save, file=config['outpath']+'speeds_'+name_addon+date+'_'+nwlst+'_'+config['freqst']+'.sav',pro_speed,ret_speed, phase_speed,$
      phase_speed_err,pro_speed_err,ret_speed_err,phase_speed_err,pro_power,ret_power,sign_power ;,pro_no_tr,ret_no_tr

if wr_speeds eq 1 then begin
 save,file=config['outpath']+'speeds_testing_'+name_addon+date+'_'+nwlst+'_'+config['freqst']+'.sav',p_vel, r_vel
endif

end
