FUNCTION debugSetPlots, nx,ny,nt,maxTrackLength,velocityImage
	; set displays for space time diagram and prograde/retrograde wave components
	window,0,xs=nx,ys=ny,retain=2,xpos=0,ypos=700
	window,1,xs=nx,ys=ny,retain=2,xpos=nx+10,ypos=700
	;window,2,xs=nt*2,ys=npt*2,title='Space-Time Diagram'
	window,3,xs=nt*2,ys=maxTrackLength*2,xpos=0,ypos=2*maxTrackLength+45,title='Prograde Filtered'
	window,4,xs=nt*2,ys=maxTrackLength*2,xpos=0,ypos=4*maxTrackLength+90,title='Retrograde Filtered'
	window,8,xs=900,ys=600

	device,decomposed=0

	wset,0    & img1=velocityImage
	blue_red & tvlct, brr,brg,brb, /get
	brr[255] = 0 & brg[255] = 0 & brb[255] = 0
	tvlct, brr,brg,brb
	img1 = bytscl(img1,-8,8,top=254)
	tv, img1
	wset,1
	loadct, 13, /silent & tvlct, r,g,b, /get
	r[255] = 255 & g[255] = 255 & b[255]= 255
	tvlct, r,g,b
	img2=bytscl(angle,-90,90,top=254)
	tv,img2

	return, {r:r,g:g,b:b,brr:brr,brg:brg,brb:brb}
END

FUNCTION debugPlot,img1,angle,tempTrack,nx,ny,nt,debugPlotColors,pro_vel,ret_vel,currentTrackLength
	; display space time diagram and prograde/retrograde wave components

	img1 = bytscl(img1,-8,8,top=254)
	img2 = bytscl(angle,-90,90,top=254)

	img1(tempTrack.xCoordinates,tempTrack.yCoordinates)=255 ;put track into displayed images
	img2(tempTrack.xCoordinates,tempTrack.yCoordinates)=255

	wset,0                      ;re-display velocity image with track
	tvlct, debugPlotColors.brr,debugPlotColors.brg,debugPlotColors.brb
	tv,img1
	wset,1                      ;re-display angle image with track
	tvlct, debugPlotColors.r,debugPlotColors.g,debugPlotColors.b
	tv,img2
	loadct, 70, /silent
	wset,3 & tvscl,bytscl(rebin(pro_vel,nt*2,currentTrackLength*2,/sample),-2,2)
	wset,4 & tvscl,bytscl(rebin(ret_vel,nt*2,currentTrackLength*2,/sample),-2,2)
	loadct, 39, /silent
	wset,8
	!p.multi=[0,2,2]
	return,0
END

FUNCTION findTrack,tempTrack,hemisphere,nx,ny,i,factor,noPixelsTrack,angle,stepSize,lowerMaskRadius,upperMaskRadius
;Calculates which x and y direction to follow

	returnValue=0
	;Use of hemisphere ensure signal always steps away from occulting disk
	IF factor EQ 0 THEN sign=-1 ELSE sign=1
	anglePrevious=tempTrack['angles',i+sign]
	tempTrack['xCoordinates',i]=0. > tempTrack['xCoordinates',i+sign]-hemisphere*stepSize*cos( (anglePrevious+factor)*!dtor) < float(nx-1)
	tempTrack['yCoordinates',i]=0. > tempTrack['yCoordinates',i+sign]-hemisphere*stepSize*sin( (anglePrevious+factor)*!dtor) < float(ny-1)

	nextX=0 > fix(tempTrack['xCoordinates',i]+0.5) < (nx-1)
	nextY=0 > fix(tempTrack['yCoordinates',i]+0.5) < (ny-1)
	    
	pixcheck=sqrt((nextX-nx/2.+0.5)^2 + (nextY-ny/2.+0.5)^2)
	IF (pixcheck lt lowerMaskRadius) or (pixcheck GT upperMaskRadius) THEN BEGIN
	   	tempTrack['xCoordinates',i]=0.
	    tempTrack['yCoordinates',i]=0.
		returnValue=1
	ENDIF

	;  IF mask[nextX,nextY] LT 1 THEN BEGIN ;don't use bad pixels
	;     xtrack(i)=0.
	;     ytrack(i)=0.
	;     BREAK
	;  ENDIF

    IF returnValue EQ 0 THEN BEGIN
		noPixelsTrack=noPixelsTrack+1
		     
		tempTrack['angles',i]=angle(nextX,nextY)  

		; if big difference in adjacent angles, resolve 180 degree ambiguity
		; by choosing angle closest to ang_track(i-1)
		angleDifference=tempTrack['angles',i]-tempTrack['angles',i+sign]			     
		IF abs(angleDifference) gt 90. THEN $
			IF angleDifference GT 0. THEN tempTrack['angles',i]-=180. $
			ELSE tempTrack['angles',i]+=180.
	ENDIF

	return, returnValue 

END


FUNCTION filterVelocity, velocityTrack,currentTrackLength,nt

	velocityTrack=velocityTrack-mean(velocityTrack)
	velocityTrack=velocityTrack-transpose(rebin(mean(velocityTrack,dim=1),currentTrackLength,nt),[1,0]) ;remove temporal mean
	 
	velocityFFT=fft(velocityTrack)

	pro_trans=velocityFFT            ;select prograde waves (assume nt even, currentTrackLength odd)
	pro_trans[1:nt/2-1, 0:currentTrackLength/2]=0.
	pro_trans[nt/2+1:-1, currentTrackLength/2+1:-1]=0.
	progradeVelocity=real_part(fft(pro_trans,/inverse))

	ret_trans=velocityFFT            ;select retrograde waves
	ret_trans[1:nt/2-1,currentTrackLength/2+1:-1]=0.
	ret_trans[nt/2+1:-1,0:currentTrackLength/2]=0.
	retrogradeVelocity=real_part(fft(ret_trans,/inverse))

	return,{retrograde:retrogradeVelocity,prograde:progradeVelocity}

END


FUNCTION calcPhaseSpeed,speeds

	relativeSpeedPrograde=speeds.prograde/speeds.progradeError
	relativeSpeedRetrograde=abs(speeds.retrograde)/speeds.retrogradeError
	relativeError=1.0/speeds.progradeError^2 + 1.0/speeds.retrogradeError^2

	speeds.phaseSpeed=(relativeSpeedPrograde^2 + relativeSpeedRetrograde^2)/relativeError
	speeds.phaseSpeedError=relativeError

	return,0
END


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
;                 maxTrackLength - length of track to use for cross-correlation, maximum useful length appears to 
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


;Define constants & variables from index
nx=index.NAXIS1          ; x dimension of velocity data cube
ny=index.NAXIS2          ; y dimension
ntime=index.NFRAMES         ; time dimension



normalCadence=index.normalCadence
spatialSampling=index.xscale

lowerMaskRadius=index.lowerMaskRadius+config['maxOcculterOffset']
upperMaskRadius=index.upperMaskRadius
mask=index.mask

maxTrackLength=config['maxTrackLength'] ;number of steps to map along track (make odd)
stepSize=1.0      ;step size along track

numberLagValues=config['numberLagValues']

;number of time points to use (make even)
if ((ntime mod 2) eq 0) then nt = ntime else nt = ntime-1
velocity=temporary(velocity[0:nx-1,0:ny-1,0:nt-1])

;===========================
;Empty track template
track=DICTIONARY("length",maxTrackLength,"angles",fltarr(maxTrackLength),$
       "xCoordinates",fltarr(maxTrackLength),"yCoordinates",fltarr(maxTrackLength) )
midPoint=fix(maxTrackLength/2)

; x = FINDGEN((nt - 1)/2) + 1
; is_nt_even = (nt MOD 2) EQ 0
; if (is_nt_even) then $
;   freq = [0.0, X, nt/2, -nt/2 + X]/(nt*normalCadence) $
; else $
;   freq = [0.0, X, -(nt/2 + 1) + X]/(nt*normalCadence)

; track.freq=rebin(freq,nt,npt)                         

; x = FINDGEN((maxTrackLength - 1)/2) + 1
; is_maxTrackLength_even = (maxTrackLength MOD 2) EQ 0
; if (is_maxTrackLength_even) then $
;   freq = [0.0, X, maxTrackLength/2, -maxTrackLength/2 + X]/(maxTrackLength*spatialSampling) $
; else $
;   freq = [0.0, X, -(maxTrackLength/2 + 1) + X]/(maxTrackLength*spatialSampling)

; ;sp_freq
; track.spatialFreq=rebin(transpose(sf),nt,npt)        

angle=median(angle,3) ;Median smoothing of wave angle map to reduce noise


; Define x and y limits for map
locationLimits=where(mask eq 1)
xstart=min(locationLimits mod nx)  
xend=max(locationLimits mod nx)
ystart=min(locationLimits/nx)
yend=max(locationLimits/nx)
pixels2process = n_elements(locationLimits )

;Already done in wave_tracking_v2
x=rebin(findgen(ny)-(nx/2.)+0.5,nx,ny)
y=transpose(x)
rad=sqrt(x^2+y^2)

;Define arrays to store prograde and retrograde quantities 
speeds={prograde:fltarr(nx,ny),retrograde:fltarr(nx,ny),phaseSpeed:fltarr(nx,ny),progradeError:fltarr(nx,ny),$
		retrogradeError:fltarr(nx,ny),phaseSpeedError:fltarr(nx,ny),sign:fltarr(nx,ny)}


IF wr_speeds EQ 1 THEN BEGIN
	 progradeCubeV=fltarr(nx,ny,nt)
	 retrogradeCubeV=fltarr(nx,ny,nt)
ENDIF


;##############################################################
;Main section of the routine
;##############################################################

message, /cont,'Starting processing...'
count=0

IF keyword_set(debug) THEN debugPlotColors=debugSetPlots(nx,ny,nt,maxTrackLength,velocity[*,*,0])

FOR yPosition=ystart,yend DO BEGIN
  FOR xPosition=xstart,xend DO BEGIN
      
      

      IF mask[xPosition,yPosition] EQ 1 THEN BEGIN 
      	 counter,count,pixels2process,/percent 
         
	  	 tempTrack=track[*] ;reinitialise temporary dictionary

	  	 tempTrack['xCoordinates',midPoint]=float(xPosition) & tempTrack['yCoordinates',midPoint]=float(yPosition)
	  	 tempTrack['angles',midPoint]=angle[xPosition,yPosition]

	  	 noPixelsTrack=0
  
	  	 ;  first, move out from cursor position 
	  	 IF xPosition le nx/2-1 THEN hemisphere=1 ELSE hemisphere=-1
	 	 factor=0
	  	 FOR i=midPoint+1,tempTrack.length-1 DO BEGIN
	  	 	res=findTrack(tempTrack,hemisphere,nx,ny,i,factor,noPixelsTrack,angle,stepSize,lowerMaskRadius,upperMaskRadius)
	  	 	IF res eq 1 THEN BREAK
	  	 ENDFOR

		 ;  next, move in from cursor position
		 factor=180
		 FOR i=midPoint-1,0,-1 DO BEGIN
		 	res=findTrack(tempTrack,hemisphere,nx,ny,i,factor,noPixelsTrack,angle,stepSize,lowerMaskRadius,upperMaskRadius)
		 	IF res eq 1 THEN BREAK
	  	 ENDFOR

         IF noPixelsTrack LT 4 THEN CONTINUE ;this tells that we are taking at least 5 points in the track
         
      	 IF noPixelsTrack LT maxTrackLength THEN BEGIN ;Trim track arrays if necessary 
		 	in=where(tempTrack.xCoordinates ne 0)
		    xTrack=tempTrack['xCoordinates',in]
		    yTrack=tempTrack['yCoordinates',in]
		    currentTrackLength=n_elements(in)
	 	 ENDIF ELSE BEGIN
	 	 	currentTrackLength=maxTrackLength
	 	 	xTrack=tempTrack['xCoordinates']
		    yTrack=tempTrack['yCoordinates']
	 	 ENDELSE                    

	 	 velocityTrack=fltarr(nt,currentTrackLength)         
		 FOR i=0,nt-1 DO velocityTrack[i,0:-1]=interpolate(velocity[0:-1,0:-1,i],xTrack,yTrack,cubic=-0.5)
	     velocityFiltered=filterVelocity(velocityTrack,currentTrackLength,nt)

		 IF keyword_set(debug) THEN resDebug=debugPlot(velocity[*,*,0],angle,tempTrack,nx,ny,nt,debugPlotColors,pro_vel,ret_vel,currentTrackLength)

		 ;##############################################################
	 	 ; compute phase speeds of prograde, retrograde and combination
	          
		  ;r=compute_speed_new(pro_vel,xscale,normalCadence,no_tr=no_tr)
		  r=compute_speed(velocityFiltered.prograde, spatialSampling, normalCadence);,numberLagValues)
		  speeds.prograde(xPosition,yPosition)=r(0)
		  speeds.progradeError(xPosition,yPosition)=r(1)
	    
		  r=compute_speed(velocityFiltered.retrograde ,spatialSampling, normalCadence,/ret);, numberLagValues,/ret)
		  ;r=compute_speed_new(ret_vel,xscale,normalCadence,/ret,no_tr=no_tr)
		  speeds.retrograde(xPosition,yPosition)=r(0)
		  speeds.retrogradeError(xPosition,yPosition)=r(1)


		 ; determine whether track is outward or inward, which determines sign of prograde and retrograde
		  speeds.sign(xPosition,yPosition)=-(rad(xTrack(0),yTrack(0))-rad(xTrack(currentTrackLength-1),yTrack(currentTrackLength-1)))/ $
	                      abs(rad(xTrack(0),yTrack(0))-rad(xTrack(currentTrackLength-1),yTrack(currentTrackLength-1)))
	     
	     IF wr_speeds eq 1 THEN BEGIN
	     	val=where(xTrack eq xPosition and yTrack eq yPosition)
		    p_vel[xPosition,yPosition,0:-1] = velocityFiltered.prograde[val]
		    r_vel[xPosition,yPosition,0:-1] = velocityFiltered.retrograde[val]
		 ENDIF
		 
		 IF keyword_set(debug) THEN !p.multi=0
	     count=count+1  
	     ;IF (1.*count)/pixels2process GT 0.2 THEN stop
	  ENDIF    
   ENDFOR
ENDFOR

result=calcPhaseSpeed(speeds)

;##############################################################
;Save outputs
;##############################################################

nwlst = strcompress(string(config['numberWaveLengths']),/rem)
save, file=config['outpath']+'speeds_'+name_addon+date+'_'+nwlst+'_'+config['freqst']+'.sav',pro_speed,ret_speed, phase_speed,$
      phase_speed_err,pro_speed_err,ret_speed_err,phase_speed_err,pro_power,ret_power,sign_power ;,pro_no_tr,ret_no_tr

if wr_speeds eq 1 then begin
 save,file=config['outpath']+'speeds_testing_'+name_addon+date+'_'+nwlst+'_'+config['freqst']+'.sav',p_vel, r_vel
endif

END

