;pro space_time_run_2,delt=delt

;## LATEST version based on Richard Morton's codes.

;  Procedure to compute space-time diagram of velocity time series over a whole image.
;  The track is computed from the angle map. This routine calls compute_speed to
;  compute the phase speed.

;  slightly changed by Hui Tian on June 13, 2012
set_plot,'x'
;compile COMPUTE_SPEED
; need to set r_inner and r_upper
r_inner=229.
r_upper=299.

debug='no'
cd,'/Data/comp_2019/'
;Device, RETAIN=2 ;If X windows protocol error, execute this command in IDL before any graphics windows are opened
;!path=!path+expand_path('+/home/htian/Scratch')
ans=' '
torad=!pi/180.
; date='20161014';'20111230';
date='20161014'
; restore,'./'+date+'/velocity_interp.sav'
restore,'./'+date+'/velocity_filtered.sav'
restore,'./'+date+'/wave_angle_new.sav'
angle=wave_angle

s=size(velo)
print,s
nx=s(1)          ;x dimension of velocity data cube
ny=s(2)          ;y dimension
ntime=s(3)     ;time dimension

if ((ntime mod 2) eq 0) then nt = ntime else nt = ntime-1    ;number of time points to use (make even)
max_npt=31
npt=max_npt    ;number of steps to map along track (make odd)

dx=1.0     ;step size along track

; xscale=3.22*dx   ;horizontal scale in Mm/sample along track
xscale=3.15*dx

xtrack=fltarr(npt)   ;x and y position along track
ytrack=fltarr(npt)
ang_track=fltarr(npt)   ;angle along track

; vmap=fltarr(nt,npt)     ;array to hold space-time diagram

delt=30.

f=findgen(nt)/(float(nt)*30.)               ;compute temporal frequency scale
for i=0,nt/2 do f(nt-i-1)=f(i+1)
freq=rebin(f,nt,npt)

sf=findgen(npt)/(float(npt)*3.22)          ;compute spatial frequency scale
for i=0,npt/2 do sf(npt-i-1)=sf(i+1)         ;make retrograde frequencies negative
sp_freq=rebin(transpose(sf),nt,npt)

phase_speed=freq/sp_freq            ;compute phase speed

;  identify acceptable frequency/phase speed regions

;good_phase=where(freq gt 0.0025 and freq lt 0.010,complement=bad_phase)
;good_phase=where(phase_speed gt 0.1 and $
; freq gt 0.0025 and freq lt 0.010,complement=bad_phase)


w=0.001      ;create filter to remove low frequency signal
filter=1.-exp(-freq^2/w^2)

angle=median(angle, 3) ;median smoothing of wave angle to reduce noise (from R. Morton's code)


pro_power=fltarr(nx,ny)     ;create arrays for prograde and retrograde power
ret_power=fltarr(nx,ny)
sign_power=fltarr(nx,ny)


mask=intarr(620,620)    ;mask out pixels to invert
x=rebin(findgen(620)-310.5,620,620)
; y=transpose(x)
y=transpose(rebin(findgen(620)-310.5,620,620))
;x=x(0:249,30:349)
;y=y(0:249,30:349)
r=sqrt(x^2+y^2)
rad=r

;good=where(r gt 230. and r lt 320. and x+3.*y gt -870,complement=bad)
good=where(r gt r_inner and r lt r_upper,complement=bad)
mask[good]=1
nmask=intarr(620,620,nt)
for tt=0,nt-1 do begin
  nmask[*,*,tt]=mask
endfor
angle(bad)=0.

bad_velo=where(finite(velo,/nan))
velo[bad_velo]=0
velo=velo*nmask


; display first velocity image for interactive input

window,0,xs=nx,ys=ny,retain=2,xpos=0,ypos=700
window,1,xs=nx,ys=ny,retain=2,xpos=nx+10,ypos=700
window,2,xs=nt*2,ys=npt*2,title='Space-Time Diagram'
window,3,xs=nt*2,ys=npt*2,xpos=0,ypos=2*npt+45,title='Prograde Filtered'
window,4,xs=nt*2,ys=npt*2,xpos=0,ypos=4*npt+90,title='Retrograde Filtered'
window,8,xs=600,ys=800
device,decomposed=0

wset,0    ;display images
img1=velo(*,*,0)
tvscl,img1<9>(-9)
wset,1
img2=bytscl(angle,-90,90)
tvscl,img2

pro_speed=fltarr(nx,ny)
ret_speed=fltarr(nx,ny)
speed=fltarr(nx,ny)
pro_speed_err=fltarr(nx,ny)
ret_speed_err=fltarr(nx,ny)
speed_err=fltarr(nx,ny)

;  x and y limits for map

; xstart=45  ;south part of cavity region
; xend=55
; ystart=330
; yend=350

; xstart=310-r_upper ;0 
; xend=310+r_upper ;nx-1
; ystart=310-r_upper;0
; yend=310+r_upper;ny-1

xstart=0 
xend=nx-1
ystart=0
yend=ny-1

 pixcounter = 0.
 old_perc_proc = 0
 tot_pix = (xend-xstart)*1L*(yend-ystart)


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


for ix=xstart,xend do for iy=ystart,yend do if angle[ix,iy] ne 0. then begin

  img1=velo(*,*,0)   ;reinitialize velocity image

  npt=31              ;redefine npt and x/ytrack due to possible change of npt when clost to the boundary.(note: it is not always good to shorten the track, npt=31 is a good option for most times)
  xtrack=fltarr(npt)   ;x and y position along track
  ytrack=fltarr(npt)   ;............................

  imid=fix(npt/2)

  xtrack(imid)=float(ix) & ytrack(imid)=float(iy)
  ang_track(imid)=angle(ix,iy)


;  use angle at each point to define track, starting at cursor position as central point
                                                                                     
if ix le nx/2-1 then leftright=1 else leftright=-1

ntrack=0
if sqrt((float(ix)-nx/2.-0.5)^2.+(float(iy)-ny/2.-0.5)^2.)-imid lt r_inner then short_track=1 else short_track=-1  ;condition for position close to the boundary

;  first, move out from cursor position

  for i=imid+1,npt-1 do begin
    xtrack(i)=0. > xtrack(i-1)-leftright*dx*cos( ang_track(i-1)*torad) < float(nx-1)
    ytrack(i)=0. > ytrack(i-1)-leftright*dx*sin( ang_track(i-1)*torad) < float(ny-1)

    intx=0 > fix(xtrack(i)+0.5) < (nx-1)
    inty=0 > fix(ytrack(i)+0.5) < (ny-1)

    pixcheck=sqrt((intx-nx/2.+0.5)^2 + (inty-ny/2.+0.5)^2)
		     IF pixcheck lt r_inner THEN BEGIN
	       	        ;this will not allow to track further down in the inward direction!
		        xtrack(i)=0.
	 	        ytrack(i)=0.

		        BREAK
	                
		     ENDIF
	   
		     IF pixcheck GT r_upper THEN BEGIN
	 	        xtrack(i)=0.
		        ytrack(i)=0.

		        BREAK
		     ENDIF
    ntrack=ntrack+1
    if debug eq 'yes' then begin
    img1(intx,inty)=min(img1) ;put track into displayed images
    endif
    ang_track(i)=angle(intx,inty)  ;  nearest neighbor interp for angle (works better than bilinear at boundaries)
;  if angle is zero (not defined), use previous value

    if ang_track(i) eq 0. then ang_track(i)=ang_track(i-1)

;  if big difference in adjacent angles, resolve 180 degree ambiguity by choosing angle closest to ang_track(i-1)

    if abs(ang_track(i)-ang_track(i-1)) gt 90. then begin
      if ang_track(i)-ang_track(i-1) gt 0. then ang_track(i)=ang_track(i)-180. else ang_track(i)=ang_track(i)+180.
    endif

  endfor

;  next, move in from cursor position

  for i=imid-1,0,-1 do begin
    ; factor=(ang_track(i+1)+180.)*torad
    xtrack(i)=0. > xtrack(i+1)-leftright*dx*cos( (ang_track(i+1)+180.)*torad) < float(nx-1)
    ytrack(i)=0. > ytrack(i+1)-leftright*dx*sin( (ang_track(i+1)+180.)*torad) < float(ny-1)

    intx=0 > fix(xtrack(i)+0.5) < (nx-1)
    inty=0 > fix(ytrack(i)+0.5) < (ny-1)

	  pixcheck=sqrt((intx-nx/2.+0.5)^2 + (inty-ny/2.+0.5)^2)
		    IF pixcheck lt r_inner THEN BEGIN
		       ;this will not allow to track further down in the inward direction!
		      xtrack(i)=0.
		      ytrack(i)=0.
			  BREAK ; stops tracks when hits lower boundary
	        ENDIF
	   
		    IF pixcheck GT r_upper THEN BEGIN
		        xtrack(i)=0.
		        ytrack(i)=0.
			    BREAK ; stops tracks when hits upper boundary
    		ENDIF
    print,'r=',sqrt((intx-nx/2.-0.5)^2.+(inty-ny/2.-0.5)^2.)
    ntrack=ntrack+1

    if debug eq 'yes' then begin
    img1(intx,inty)=min(img1) ;put track into displayed images
    endif


    ang_track(i)=angle(intx,inty)  ;  nearest neighbor interp for angle (works better than bilinear at boundaries)

;  if angle is zero (not defined), use previous value

    if ang_track(i) eq 0. then ang_track(i)=ang_track(i+1)

;  if big difference in adjacent angles, resolve 180 degree ambiguity by choosing angle closest to ang_track(i-1)

    if abs(ang_track(i)-ang_track(i+1)) gt 90. then begin
      if ang_track(i)-ang_track(i+1) gt 0. then ang_track(i)=ang_track(i)-180. else ang_track(i)=ang_track(i)+180.
    endif

  endfor

  if ntrack lt 4 then continue
  if ntrack mod 2 eq 0 then ntrack=ntrack+1   ;make it odd

  if ntrack lt max_npt then begin
  in=where(xtrack ne 0)
  newxtrack=xtrack(in)
  xtrack=newxtrack
  newytrack=ytrack(in)
  ytrack=newytrack
  npt=n_elements(in)
  ;ang_track=ang_track(in)
  endif else begin
  npt=max_npt
  endelse
  vmap=fltarr(nt,npt)
;  interpolate velo onto track

  for i=0,nt-1 do begin
    vel=velo(*,*,i)
    v=interpolate(vel,xtrack,ytrack,cubic=-0.5)
    vmap(i,*)=v
  endfor

  if debug eq 'yes' then begin
    wset,0       ;redisplay velocity image with track
  tvscl,img1
  endif

vmap=vmap-mean(vmap)
vmap=vmap-transpose(rebin(mean(vmap,dim=1),npt,nt),[1,0])




;  filter velocity map

  ;for i=0,npt-1 do vmap(*,i)=vmap(*,i)-mean(vmap(*,i)) ;remove temporal mean

  trans=fft(vmap,-1)       ;compute fourier transform

  mag=abs(trans)         ;compute magnitude and phase
  phase=atan(imaginary(trans),float(trans))

  ;back=rebin( rebin(mag(nt/2-50:nt/2+50,*),1,npt) ,nt,npt)
  ;back=rebin( rebin(mag(nt/2-nt/4:nt/2+nt/4,*),1,npt) ,nt,npt)
  back=rebin( rebin(mag(nt/2-(nt/2-2):nt/2+(nt/2-2),*),1,npt) ,nt,npt)
  mag=mag-back       ;remove high temporal frequency noise

  back=rebin(mag(*,npt/2),nt,npt)
  mag=mag-back       ;remove high spatial frequency noise

  mag=mag > 0.   ;insure that magnitude is positive

  pw_trans=complex(mag*cos(phase),mag*sin(phase))   ;recompute transform


  pw_trans=pw_trans*filter       ;remove low temporal frequencies


  pro_trans=trans             ;select prograde waves (assume nt even, npt odd)
  pro_trans(1:nt/2-1,0:npt/2)=0.
  pro_trans(nt/2+1:nt-1,npt/2+1:npt-1)=0.
  ;pro_trans(bad_phase)=0.
  pro_vel=real_part(fft(pro_trans,1))

	pro_pw_trans=pw_trans             ; separate calculation for the wave power
	pro_pw_trans(1:nt/2-1,0:npt/2)=0.
	pro_pw_trans(nt/2+1:nt-1,npt/2+1:npt-1)=0.

  ret_trans=trans             ;select retrograde waves
  ret_trans(1:nt/2-1,npt/2+1:npt-1)=0.
  ret_trans(nt/2+1:nt-1,0:npt/2)=0.
  ;ret_trans(bad_phase)=0.
  ret_vel=real_part(fft(ret_trans,1))

  ret_pw_trans=pw_trans             ; separate calculation for the wave power
	ret_pw_trans(1:nt/2-1,npt/2+1:npt-1)=0.
	ret_pw_trans(nt/2+1:nt-1,0:npt/2)=0.

  ; trans(bad_phase)=0.
  ; vel=float(fft(trans,1))


;  display space time diagram and prograde/retrograde wave components
if debug eq 'yes' then begin
  ;   wset,2
  ; tvscl,bytscl(rebin(vel,2*nt,2*npt,/sample),-2,2)
  wset,3
  tvscl,bytscl(rebin(pro_vel,nt*2,npt*2,/sample),-2,2)
  wset,4
  tvscl,bytscl(rebin(ret_vel,nt*2,npt*2,/sample),-2,2) 
   wset,8
  !p.multi=[0,3,2,0,1]
endif




;  compute phase speeds of prograde, retrograde and combination



  r=compute_speed_1(pro_vel,xscale,delt)
  pro_speed(ix,iy)=r(0)
  pro_speed_err(ix,iy)=r(1)

  r=compute_speed_1(ret_vel,xscale,delt,/ret)
  ret_speed(ix,iy)=r(0)
  ret_speed_err(ix,iy)=r(1)

  phase_speed(ix,iy)=(pro_speed(ix,iy)/pro_speed_err(ix,iy)^2 + abs(ret_speed(ix,iy))/$
	                      ret_speed_err(ix,iy)^2)/(1.0/pro_speed_err(ix,iy)^2 + 1.0/ret_speed_err(ix,iy)^2)
	phase_speed_err(ix,iy)=(1.0/pro_speed_err(ix,iy)^2 + 1.0/ret_speed_err(ix,iy)^2)^(-0.5)

;  print,pro_speed(ix,iy),pro_speed_err(ix,iy),ret_speed(ix,iy),ret_speed_err(ix,iy),$
;   speed(ix,iy),speed_err(ix,iy)

;  compute prograde and retrograde power

  pro_pow=abs(pro_pw_trans)^2
  ;pro_power(ix,iy)=mean(pro_pow(good_phase))
  pro_power(ix,iy)=total(pro_pow)
  ret_pow=abs(ret_pw_trans)^2
  ;ret_power(ix,iy)=mean(ret_pow(good_phase))
  ret_power(ix,iy)=total(ret_pow)

;  determine whether track is outward or inward, which determines sign of prograde and retrograde

  sign_power(ix,iy)=-(rad(xtrack(0),ytrack(0))-rad(xtrack(npt-1),ytrack(npt-1)))/ $
   abs(rad(xtrack(0),ytrack(0))-rad(xtrack(npt-1),ytrack(npt-1)))

  !p.multi=0
endif


; save,file='./'+date+'/1074.l2/proret.sav',pro_speed,ret_speed,speed,pro_speed_err,ret_speed_err,speed_err,pro_power,ret_power,sign_power
save,file='./'+date+'/proret_new_test.sav',pro_speed,ret_speed,speed,pro_speed_err,ret_speed_err,speed_err,pro_power,ret_power,sign_power,trans,mag, phase, pw_trans

end
