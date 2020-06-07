pro space_time_interact_tracks, date, frequency=frequency, nwl=nwl, tracklen=tracklen, on_velo=on_velo, on_enh_int=on_enh_int

; PURPOSE:  Procedure to compute space-time diagram of velocity time series. The central point of
;  the track is selected interactively and the track is computed from the angle map. This
;  routine calls compute_speed to compute the phase speed and display_kimage to display
;  the k-omega diagram.
;
;
;INPUTS: date - date of observations
;        frequency - select which wave angle map to use, e.g. frequency=3.5 will look for map calculated with data
;                    filtered at 3.5 mHz from wave_tracking.pro
;
;OPTIONAL INPUTS: nwl - number of wavelength positions, default is 3
;                 tracklen - maximum length of space-time track, default is 51. The routine will always make the
;                            track length odd, even if an even number is input
;                 /on_velo - plots velocity map instead of intensity map
;                 /on_enh_int - plots enhanced intensity map instead of intensity map - only if on HAO server
;
;EXAMPLE CALL:    space_time_interact_tracks,'20120327',freq=3.5,tracklen=31
;
;                 
;
;RESTRICTIONS: Can only use certain commands if on HAO servers
;              Need to define the directory for COMP data on personal machine in the code
;
;CALLS: compute_speed.pro
;
;HISTORY: Author: S. Tomczyk
;  	  Modified by: C. Bethge
;         R Morton 11/2014 - Fixed minor bug for dealing with tracks that are shorter than specified track length.
;                            Commented code.
;         R Morton 2016 - updated code to work with output from wave_tracking_v2                 
;



;#####################################
;Restore files from wave_tracking.pro
;#####################################

IF n_elements(nwl) EQ 0 THEN nwl = 3
IF n_elements(frequency) EQ 0 THEN frequency = -1
nwlst = strcompress(string(nwl),/rem)

outpath = 'analysis/comp/wave_tracking_output/'+date+'/'   ;<- same as in wave_tracking.pro


ans=' '
torad=!pi/180.

;=======================================
; *** Things that need to be defined ***
;=======================================
case frequency of 
  1.5: freqst = '1.5mHz'
  3.5: freqst = '3.5mHz'
  5.5: freqst = '5.5mHz'
  else: begin
         freqst = 'manual_freq_setting'  
        end
endcase

restore,  outpath+'cube_ivw_'+date+'.sav'
restore,  outpath+'wave_angle_'+date+'_'+nwlst+'_'+freqst+'.sav'

;#####################################
;USER DEFINED BOUNDARIES - i.e., where occulting disk is etc.
;#####################################
case date of
'20120327': begin
             lower_r = 229.
             upper_r = 276.
            end
 else:      begin
             lower_r = 230.
             upper_r = 292.
            end
endcase



;#####################################
;if you are not at HAO and want to use the enhanced intensity, you
;need to provide a valid path for img1 here (also if you *are* at
;HAO and want a different file than he first 3pt file of the day)
;#####################################

IF KEYWORD_SET(on_enh_int) THEN BEGIN
  datpath  = '/hao/acos/'+strmid(date,8,4,/reverse)+'/'+strmid(date,3,2,/reverse)+'/'+strmid(date,1,2,/reverse)+'/'
  flist = file_list(datpath,'*comp.1074.dynamics.3.fts.gz')
  img1 = readfits(flist[0],ext=2)

ENDIF

velocity = cube_v
angle=wave_angle

s=size(velocity)
nx=s(1)          ;x dimension of velocity data cube
ny=s(2)          ;y dimension


;#####################################
;number of time points to use (make even)
;#####################################
IF (s[3] mod 2) EQ 1 THEN nt=s[3]-1 ELSE nt=s[3]    


;#####################################
;Define track length
;#####################################

IF n_elements(tracklen) EQ 0 THEN BEGIN
   npti = 51  ;number of steps to map along track (make odd)
ENDIF ELSE BEGIN
  IF (tracklen MOD 2) EQ 0 THEN tracklen = tracklen+1
  npti = tracklen
ENDELSE

npt=npti   
dx=1.0     ;step size along track (pixels)


;==============================================
; *** End of things that need to be defined ***
;==============================================


;#####################################
;Zero pixels that lie outside boundaries 
;#####################################

xx=rebin(findgen(ny)-(nx/2.)+0.5,nx,ny)
yy=transpose(xx)
rr=sqrt(xx^2+yy^2)
good=where(rr ge lower_r and rr le upper_r and index.mask gt 0, comp=bad)

FOR ii=0,nt-1 DO BEGIN
 temp = reform(cube_v[*,*,ii])
 temp[bad] = 0
 cube_v[*,*,ii] = temp
ENDFOR
angle[bad] = 0


;#####################################
;Set arrays to store variables
;#####################################
ang_track=fltarr(npt)   ;angle along track
vmap=fltarr(nt,npt)     ;array to hold space-time diagram

f=findgen(nt)/(float(nt)*index.norm_cadence)   ;compute temporal frequency scale
FOR i=0,nt/2 DO f(nt-i-1)=f(i+1)
freq=rebin(f,nt,npt)            


;#####################################
;Create filter to remove low frequency signal
;#####################################
w=0.001     
filter=1.-exp(-freq^2/w^2)

;freq=findgen(nt)/(float(nt)*30.0)     ;compute frequency scale
;for i=0,nt/2 do freq(nt-i-1)=freq(i+1)
;good_phase=where(freq gt 0.0025 and freq lt 0.010,complement=bad_phase)    
; taking cutoff in both side when power is calculated!


;#####################################
;Load up windows for interactive display
;#####################################

; display
loadct, 13, /silent
window,0,xs=nx,ys=ny,retain=2,xpos=0,ypos=400
window,1,xs=nx,ys=ny,retain=2,xpos=nx,ypos=400
window,2,xs=nt*2,ys=npt*2,title='Space-Time Diagram',xpos=0, ypos=2*(npt*2)+66
window,3,xs=nt*2,ys=npt*2,title='Prograde Filtered',xpos=0, ypos=(npt*2)+44
window,4,xs=nt*2,ys=npt*2,title='Retrograde Filtered',xpos=0, ypos=22
window,5,xs=900,ys=600,xpos=2*nt,ypos=0
;window,6,xpos=2*nx,ypos=0
device, decomposed=0

wset,0    ;display images
IF KEYWORD_SET(on_velo) THEN BEGIN
   img1=velocity[*,*,0]
   blue_red
   tvlct, brr,brg,brb, /get
   brr[255] = 0 & brg[255] = 0 & brb[255] = 0 
   tvlct, brr,brg,brb
   img1 = bytscl(img1,-8,8,top=254)
ENDIF ELSE BEGIN
   IF KEYWORD_SET(on_enh_int) THEN BEGIN
      img1 = bytscl(img1,top=254)
   ENDIF ELSE BEGIN
      img1=cube_i[*,*,0]
      img1 = bytscl(sqrt(img1),min=1,max=5,top=254)
   ENDELSE
   loadct, 0
   tvlct, brr,brg,brb, /get
   brr[255] = 255 & brg[255] = 0 & brb[255] = 0 
   tvlct, brr,brg,brb
ENDELSE
tv, img1
wset,1
loadct, 4, /silent
tvlct, war,wag,wab, /get
war[255] = 255 & wag[255] = 255 & wab[255]= 255 
tvlct, war,wag,wab
img2=bytscl(angle,-90,90,top=254)
tv,img2

track_counter = 0
track_counter_old = 0
wave_info = ''
wave_info_old = ''

continue_choosing_tracks:
!mouse.button = 0
quit_trigger  = 0
WHILE quit_trigger EQ 0 DO BEGIN
  print,'******************************************************************'
  print,'* move cursor to region to trace and press left button           *'
  print,'* click right to undo last selection                             *' 
  print,'* move cursor to upper left corner and press left button to quit *'
  print,'******************************************************************'
  wset,1
  cursor,ix,iy,/down,/dev     ;get starting position from cursor

  IF !mouse.button EQ 4 THEN BEGIN
    track_counter = track_counter_old
    wave_info = wave_info[0:track_counter-1]
    wset,0
    tvlct, brr,brg,brb
    tv, img1
    FOR ii=0,track_counter-1 DO plots, wave_info[ii].xtrack, wave_info[ii].ytrack, psym=3, color=255, /device
    wset,1
    tvlct, war,wag,wab
    
    FOR ii=0,track_counter-1 DO plots, wave_info[ii].xtrack, wave_info[ii].ytrack, psym=3, color=255, /device
    print, 'Number of selected tracks: '+string(track_counter)
    goto, sayonara
  ENDIF ELSE BEGIN
    track_counter_old = track_counter
    IF (ix le 9 and iy ge 610) THEN BEGIN
        quit_trigger = 1
        goto, sayonara
    ENDIF

    track_counter = track_counter+1   
    print, 'Number of selected tracks: '+string(track_counter)
 
    npt=npti       ;we need to define these again
    imid=fix(npt/2)
    xtrack=fltarr(npt)   ;x and y position along track
    ytrack=fltarr(npt)
    xtrack(imid)=float(ix) & ytrack(imid)=float(iy)
    ang_track=fltarr(npt) 
    ang_track(imid)=angle(ix,iy)

    ;Here the condition of the each hemisphere is taken which will be used later.
    IF ix LE nx/2-1 THEN hemisphere=1 ELSE hemisphere=-1

    ; use angle at each point to define track, starting at cursor position as central point

    ntrack=0                      ;this index will count the number of pixels on the track
    IF sqrt( (float(ix)-nx/2.+0.5)^2 + (float(iy)-ny/2.+0.5)^2 )-imid LT lower_r THEN $
     short_track=1 ELSE short_track=-1 ; this condition gives the lower limit of the track 
    ; this value (225) we have to take carefully!
    ; first, move out from cursor position

    FOR i=imid+1,npt-1 DO BEGIN

        xtrack(i)=0. > xtrack(i-1)-hemisphere*dx*cos( ang_track(i-1)*torad) < float(nx-1)  
        ;This means value should always greater
        ytrack(i)=0. > ytrack(i-1)-hemisphere*dx*sin( ang_track(i-1)*torad) < float(ny-1)  
        ;than 0 and less than nx-1

        intx=0 > fix(xtrack(i)+0.5) < (nx-1)
        inty=0 > fix(ytrack(i)+0.5) < (ny-1)

        IF short_track EQ 1 THEN BEGIN                                         
           ; this means if u are near the lower boundary
           IF sqrt( (intx-nx/2.+0.5)^2 + (inty-ny/2.+0.5)^2 ) LT lower_r THEN BREAK  
           ; this means it will not allow to track further down in the inward direction!
        ENDIF                                                                 
        ; actually this is not required, but sometimes the (in the loop) 
        ; move outward also can go inside
        
        ntrack=ntrack+1

        ang_track(i)=angle(intx,inty)  
        ; nearest neighbor interp for angle (works better than bilinear at boundaries)

        ; if angle is zero (not defined), use previous value
        IF ang_track(i) EQ 0. THEN ang_track(i)=ang_track(i-1)

       ; if big difference in adjacent angles, resolve 180 degree ambiguity by choosing angle closest to ang_track(i-1)
       IF abs(ang_track(i)-ang_track(i-1)) gt 90. THEN BEGIN
             if ang_track(i)-ang_track(i-1) gt 0. THEN ang_track(i)=ang_track(i)-180. ELSE ang_track(i)=ang_track(i)+180.
       ENDIF  

    ENDFOR

    ; next, move in from cursor position

    FOR i=imid-1,0,-1 DO BEGIN

        xtrack(i)=0. > xtrack(i+1)-hemisphere*dx*cos( (ang_track(i+1)+180.)*torad) < float(nx-1)
        ytrack(i)=0. > ytrack(i+1)-hemisphere*dx*sin( (ang_track(i+1)+180.)*torad) < float(ny-1)

        intx=0 > fix(xtrack(i)+0.5) < (nx-1)
        inty=0 > fix(ytrack(i)+0.5) < (ny-1)

        ;print,sqrt( (float(intx)-nx/2.+0.5)^2 + (float(inty)-ny/2.+0.5)^2 )
        IF short_track EQ 1 THEN BEGIN       
        ; this means if you are near the lower boundary
          IF sqrt( (intx-nx/2.+0.5)^2 + (inty-ny/2.+0.5)^2 ) LT lower_r THEN BREAK  
          ;this means it will not allow to track further down in the inward direction!
        ENDIF
        ;+++++++++++++++++++++
        ;UNDER TEST - RM
        ;Some tracks reach boundary and then do a u-turn - trying to stop this
        ;doesn't work well for cureved loops
        ;r=sqrt( (float(intx)-nx/2.+0.5)^2 + (float(inty)-ny/2.+0.5)^2 )
        ;r1=sqrt( (float(xtrack(i+1))-nx/2.+0.5)^2 + (float(ytrack(i+1))-ny/2.+0.5)^2 )
        ntrack=ntrack+1
        
        ;IF r gt r1 THEN BREAK
        

        ang_track(i)=angle(intx,inty)  
 
        ; nearest neighbor interp for angle (works better than bilinear at boundaries)
        ; if angle is zero (not defined), use previous value
        IF ang_track(i) eq 0. THEN ang_track(i)=ang_track(i+1)

        ; if big difference in adjacent angles, resolve 180 degree ambiguity by choosing angle closest to ang_track(i-1)
        IF abs(ang_track(i)-ang_track(i+1)) GT 90. THEN BEGIN
          IF ang_track(i)-ang_track(i+1) GT 0. THEN ang_track(i)=ang_track(i)-180. ELSE ang_track(i)=ang_track(i)+180.
        ENDIF
        ;IF ang_track(i) GT 70. THEN BREAK
    ENDFOR

    IF ntrack LT 4 THEN continue    ;this tells that we are taking atleast 5 points in the track

    ; ntrack=ntrack+1   
    ; this gives the total number of pixels
    ; in the track. for normal case it is =npt,
    ; but when the track is shorter it is < npt
    IF ntrack NE ntrack/2*2+1 THEN ntrack=ntrack/2*2+1
    ;print,ntrack,npt
     
    IF ntrack LT npt THEN BEGIN    ;......If the track is shorter, then it will do the followings
       in=where(xtrack ne 0)
       newxtrack=xtrack(in) ;......it creates new track of shorter length
       xtrack=newxtrack
       newytrack=ytrack(in)
       ytrack=newytrack
       npt=n_elements(in)
       ang_track=ang_track(in)
    ENDIF ELSE BEGIN
       npt=npti                    ; We need the following lines for the new track
    ENDELSE
    ;print,npt,ntrack,n_elements(xtrack)
    vmap=fltarr(nt,npt)          ; We need to make it again (array to hold space-time diagram)
    

    ; interpolate velocity onto track

    FOR i=0,nt-1 DO BEGIN
        vel=velocity(*,*,i)
        v=interpolate(vel,xtrack,ytrack,cubic=-0.5)
        vmap(i,*)=v
    ENDFOR
    vmap=vmap-mean(vmap)
   
    wset,0                       ;redisplay images with track
    tvlct, brr,brg,brb
    tv,img1
    IF track_counter GT 1 THEN BEGIN 
       FOR ii=0,track_counter-2 DO plots, wave_info[ii].xtrack, wave_info[ii].ytrack, psym=3, color=255, /device
    ENDIF
    plots, xtrack, ytrack, psym=3, color=255, /device
    wset,1
    tvlct, war,wag,wab
    tv,img2
    IF track_counter GT 1 THEN BEGIN
        FOR ii=0,track_counter-2 DO plots, wave_info[ii].xtrack, wave_info[ii].ytrack, psym=3, color=255, /device
    ENDIF  
    plots, xtrack, ytrack, psym=3, color=255, /device

   

    ;#####################################
    ; filter velocity map
    ;#####################################

    ;remove temporal mean
    FOR i=0,npt-1 DO vmap(*,i)=vmap(*,i)-mean(vmap(*,i)) 

    trans=fft(vmap,-1)            ;compute fourier transform

    mag=abs(trans)         ;compute magnitude and phase
    phase=atan(imaginary(trans),float(trans))
 
    back=rebin( rebin(mag(nt/2-(nt/2-2):nt/2+(nt/2-2),*),1,npt) ,nt,npt)
    mag=mag-back       ;remove high temporal frequency noise
 
    back=rebin(mag(*,npt/2),nt,npt)
    mag=mag-back       ;remove high spatial frequency noise
 
    mag=mag > 0.   ;insure that magnitude is positive
 
    pw_trans=complex(mag*cos(phase),mag*sin(phase))   ;recompute transform
 
    pw_trans=pw_trans*filter       ;remove low temporal frequencies;  filter is defined at top

    ;#####################################
    ;select prograde waves (assume nt even, npt odd)
    ;#####################################
    pro_trans=trans               
    pro_trans(1:nt/2-1,0:npt/2)=0.
    pro_trans(nt/2+1:nt-1,npt/2+1:npt-1)=0.
    pro_vel=float(fft(pro_trans,1))

    pro_pw_trans=pw_trans             ; separate calculation for the wave power
    pro_pw_trans(1:nt/2-1,0:npt/2)=0.
    pro_pw_trans(nt/2+1:nt-1,npt/2+1:npt-1)=0.

    ;#####################################
    ;select retrograde waves
    ;#####################################
    ret_trans=trans               
    ret_trans(1:nt/2-1,npt/2+1:npt-1)=0.
    ret_trans(nt/2+1:nt-1,0:npt/2)=0.
    ret_vel=float(fft(ret_trans,1))

    ret_pw_trans=pw_trans             ; separate calculation for the wave power
    ret_pw_trans(1:nt/2-1,npt/2+1:npt-1)=0.
    ret_pw_trans(nt/2+1:nt-1,0:npt/2)=0.

    vel=float(fft(trans,1))


    ; display space time diagram and prograde/retrograde wave components
    loadct, 13, /silent
    wset,2
    tvscl,bytscl(rebin(vel,2*nt,2*npt,/sample),-2,2)
    wset,3
    tvscl,bytscl(rebin(pro_vel,nt*2,npt*2,/sample),-2,2)
    wset,4
    tvscl,bytscl(rebin(ret_vel,nt*2,npt*2,/sample),-2,2)
    loadct, 39, /silent
    wset,5
    !p.multi=[0,2,2]

    ;#####################################
    ; compute phase speeds of prograde, retrograde and combination
    ;#####################################

    pro_speed=compute_speed(pro_vel,index.xscale,index.norm_cadence,debug=1)
    ret_speed=compute_speed(ret_vel,index.xscale,index.norm_cadence,/ret,debug=1)

    phase_speed = pro_speed-pro_speed
    phase_speed[0]=(pro_speed[0]/pro_speed[1]^2 + abs(ret_speed[0])/$
                      ret_speed[1]^2)/(1.0/pro_speed[1]^2 + 1.0/ret_speed[1]^2)
    phase_speed[1]=(1.0/pro_speed[1]^2 + 1.0/ret_speed[1]^2)^(-0.5)

    ;#####################################
    ; compute prograde and retrograde power
    ;#####################################
    pro_pow=abs(pro_pw_trans)^2
    ;pro_power=mean(pro_pow(good_phase))
    pro_power=mean(pro_pow)
    ret_pow=abs(ret_pw_trans)^2
    ;ret_power=mean(ret_pow(good_phase))
    ret_power=mean(ret_pow)

    !p.multi=0
  
    ;#####################################
    ; compute k-omega diagram (not used in current output)
    ;#####################################
    
    power=abs(trans)^2
    power=power(*,0:npt/2-1)
    power=shift(power,nt/2,0)

    power=rebin(power,3*nt,6*fix(npt/2),/sample) > 0.
    power=alog10(power)

    ;  phase_sp=(abs(pro_speed[0])+abs(ret_speed[0]))/2.
    loadct, 0, /silent
    ;  display_kimage,power,-6.,-4.,tit='k-!4x!3 Diagram',units='log!d10!n(Power) (km!u2!n s!u-2!n)',$
    ;   table=0,output='s',windnum=6,wherebar='right',scale=1.,xt='Frequency (Hz)',$
    ;   yt='Spatial Frequency (Mm!u-1!n)',phase_speed=phase_speed[0]

    print,n_elements(xtrack)
    for ij=0,track_counter-2 do print,n_elements(wave_info[ij].xtrack)

    str = {power:power, pro_power:pro_power, ret_power:ret_power, vel:vel, pro_vel:pro_vel, ret_vel:ret_vel, $
         pro_speed:pro_speed, ret_speed:ret_speed, phase_speed:phase_speed, xtrack:xtrack, ytrack:ytrack, $
         ang_track:ang_track, nt:nt, npt:npt}
    IF track_counter EQ 1 THEN wave_info = str
    IF track_counter GT 1 THEN wave_info = merge_struct(wave_info, str)
  ENDELSE

sayonara:
ENDWHILE

answ  = ''
answs = ''
ansfn = ''
new_answ:  
read, answ, prompt='Do you really want to quit? (y/n): '  
IF answ NE 'y' AND answ ne 'n' THEN goto, new_answ
IF answ EQ 'n' THEN goto, continue_choosing_tracks

IF track_counter GE 1 THEN BEGIN
   new_answs:
   read, answs, prompt='Save results? (y/n): '  
   IF answs NE 'y' AND answs NE 'n' THEN goto, new_answs
   IF answs eq 'y' THEN BEGIN
      new_ansfn:
      read, ansfn, prompt='Please type in filename (e.g. mytracks.sav): '  
      IF strmid(ansfn,3,4,/reverse) ne '.sav' THEN BEGIN
         print, 'filename has to end with .sav, please type again'
         goto, new_ansfn
      ENDIF
      velocity   = img1
      wave_angle = img2
      save, wave_angle, velocity, wave_info, filename=outpath+ansfn
    ENDIF 
ENDIF

FOR ii=0,5 DO wdelete, ii

END

