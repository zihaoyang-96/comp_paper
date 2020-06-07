pro find_angle_v2,spec,dt,xyrange,wave_angle,angle_error

;## LATEST version based on Richard Morton's codes.

;  Procedure to determine the wave propagation angle from the cross correlation of the
;  velocity time series of a test point with adjacent time series. The data cube of
;  velocity spectra is input and the computation is performed in the Fourier domain.
;  Only contiguous points having coherence greater than a threshold (typically 0.5)
;  are used. The coherence is computed over a range defined by a frequency filter.
;
;  modified from find_angle.pro (Steve) by Hui Tian on June 12, 2012
;  input:
;  spec: fft spectrum of the velocity cube
;  dt: time spacing
;  xyrange: only do the calculation for the region of [xmin,xmax,ymin,ymax], in the unit of pixel
;  output: wave_angle,angle_error
;
rmin=229.
rmax=299.

ans=' '
torad=!pi/180.
debug='no'; 'yes';  ;debug mode, yes or no

s=size(spec)
if debug eq 'yes' then print,s
nx=s(1)
ny=s(2)
nspec=s(3)

;  create filter
;  gaussian filter peaked at freq0 with width of fwidth
freq=findgen(nspec)/(float(nspec*2)*dt) ;note that nspec=ntime/2, see prewavetracking.pro
nfreq=n_elements(freq)
freq0=.0035
fwidth=.0015
filter=exp( -(freq-freq0)^2/fwidth^2 )
filter(0)=0.    ;set dc to zero
filter(where (freq lt .001))=0.     ;set low frequencies to zero
filter=filter/total(filter)

if debug eq 'yes' then begin
  window
  plot,freq,filter
  read,'enter return',ans
endif

dx=20
nbox=2*dx+1   ;box size in pixels

nsmooth= 15;21 ;51  ;smoothing width for cross spectra (nominally 11)
; nsmooth=3
limit= 0.5;0.6; 0.9;    ;set coherence limit for acceptance

mask=intarr(620,620)    ;mask out pixels to invert
x=rebin(findgen(620)-310.5,620,620)
y=transpose(rebin(findgen(620)-310.5,620,620))
r=sqrt(x^2+y^2)
good=where(r gt rmin and r lt rmax)
mask(good)=1

xstart=xyrange[0]
xend=xyrange[1]
ystart=xyrange[2]
yend=xyrange[3]

; open windows
if debug eq 'yes' then begin
  f=8
window,0,xs=nx,ys=ny,retain=2,xpos=400,ypos=0
window,1,xs=nbox*f,ys=nbox*f,xpos=0,ypos=0
window,2,xs=nbox*f,ys=nbox*f,xpos=0,ypos=f*nbox+50
window,3,xs=nbox*f,ys=nbox*f,xpos=0,ypos=2*f*nbox+100
window,4,xs=nx,ys=ny,retain=2,xpos=400+nx+20,ypos=0
window,5,xs=nx,ys=ny,retain=2,xpos=400,ypos=ny+50
window,6,xs=512,ys=512,xpos=400,ypos=800
device,decomposed=0
endif

if debug eq 'no' then begin
  pixcounter = 0L
  old_perc_proc=0
  tot_pix=long(n_elements(where (mask eq 1)))
endif

conjspec=conj(spec)

lar_g=fltarr(nx,ny,nspec)
FOR iy=0,ny-1 DO FOR ix=0,nx-1 DO $
    IF (mask[ix,iy] EQ 1) THEN lar_g[ix,iy,*]=smooth(real_part(reform(spec[ix,iy,*]*conjspec[ix,iy,*])),nsmooth)

xbox=rebin(findgen(nbox)-dx,nbox,nbox)  ;coordinates of points in box
ybox=transpose(xbox)

wave_angle=fltarr(nx,ny)
angle_error=fltarr(nx,ny)
coh_measure = fltarr(nx,ny,2)


;calculates rms value of spectrum
;Used to provide secondary condition on for loops
;and exclude pixels with low signal
rms_spec=total(abs(spec),3)

for ix=xstart,xend do for iy=ystart,yend do if (mask(ix,iy) eq 1) and (rms_spec[ix,iy] gt 1e-7) then begin

  if debug eq 'no' then begin
    perc_proc=round((float(pixcounter)/float(tot_pix))*100.)
    pixcounter=pixcounter+1
    IF perc_proc GT old_perc_proc THEN BEGIN
            IF strmid(getenv('TERM'),0,5) EQ 'xterm' THEN BEGIN
     		progress_comp, perc_proc, 100, msg='Computing wave propagation angles'
    	    ENDIF ELSE message, /cont, strcompress(string(perc_proc,format='(I3)'),/rem)+'% processed.'
    	    old_perc_proc = perc_proc
     ENDIF
  endif
  
;  compute coherence for contiguous pixels with coherence greater than limit

  coh=fltarr(nbox,nbox)           ;initialize coherence array to zero
  coh_mask=fltarr(nbox,nbox)      ;initialize mask of good points in coherence

  spec1=reform(spec[ix,iy,*])
  

  ;g1=real_part(spec1*conjspec)
  ;lar_g=smooth(real_part(spec*conjspec),[1,1,nsmooth])
  ;g1=smooth(g1,[1,1,nsmooth])
  g1=reform(lar_g[ix,iy,*])

  coh(nbox/2,nbox/2)=1.
  coh_mask(nbox/2,nbox/2)=1.

  count=1
  ncoh=1
  while count gt 0 and ncoh lt 100 do begin
    cont=where(coh_mask eq 1.)
    cnew=[cont-1,cont+1,cont+nbox,cont-nbox]  ;test pixels adjacent to contiguous pixels

    new_mask=fltarr(nbox,nbox)   ;indentify new contiguous pixels
    new_mask(cnew)=1.
    new_mask=new_mask*(1.-coh_mask)

    new=where(new_mask gt 0.,new_count)

    if debug eq 'yes' then begin
      wset,2
      tv,bytscl(rebin(new_mask,f*nbox,f*nbox,/sample),0.,1.)
    endif

    count=0
    for ic=0,new_count-1 do begin
      icx=new(ic) mod nbox
      icy=fix(new(ic)/nbox)

      ii=ix-nbox/2+icx
      jj=iy-nbox/2+icy

      if (ii lt 0) or (ii gt nx-1) or (jj lt 0) or (jj gt ny-1) then coh(icx,icy)=0. else begin
        if coh(icx,icy) eq 0. then begin
          ; spec2=spec(ii,jj,*)
          spec2c=reform(conjspec(ii,jj,*))
          cspec=spec1*spec2c ;cross-spectrum: FFT of reference pixel * conjugation of FFt of other pixel
          cspec=smooth(cspec,nsmooth)
          g2=lar_g[ii,jj,*]
          ; g2=smooth(g2,nsmooth)

          coh(icx,icy)=mask(ii,jj)*total( filter[1:-1]*( abs(cspec[1:-1])/sqrt(g1[1:-1]*g2[1:-1]) ) ,/nan) ;see McInotosh et al. 2008. Sol. Phys.
        endif

        if coh(icx,icy) gt limit then begin
          coh_mask(icx,icy)=1.
          count=count+1
        endif

      endelse

    endfor

    if debug eq 'yes' then begin
      wset,3
      tv,bytscl(rebin(coh_mask,f*nbox,f*nbox,/sample),0.,1.)
      read,'enter return',ans
    endif

   good=where(coh_mask eq 1.,ncoh)

   endwhile

   ; wset,1
   ; tv,bytscl(rebin(coh*coh_mask,f*nbox,f*nbox,/sample),0.,1.)
   ;IF (keyword_set(save_coh) and (ncoh gt npix)) THEN BEGIN
  ;  IF (ncoh gt 10) THEN BEGIN
  ;    island = coh*coh_mask
  ;    ellipsepts = fit_ellipse(where(island ne 0.),axes=axes,xsize=nbox,ysize=nbox)
  ;    coh_measure[ix,iy,*] = axes
  ;    if ((size(ellipsepts))[0] gt 0) and debug eq 'yes' then plots, f*ellipsepts, color=254, /device
  ;  ENDIF

;  find angle of coherence island if enough good points

  if ncoh gt 10 then begin
    if debug eq 'yes' then begin
      wset,6
    plot,xbox(good),ybox(good),psym=6,xr=[-dx,dx],yr=[-dx,dx]

    endif
    
    
    weight=coh(good)^2
    theta=min_dist_fit(xbox(good),ybox(good),error=error) ;analytical solution to the minimum distance problem
    print,theta,error

    wave_angle(ix,iy)=theta
    angle_error(ix,iy)=error
    if debug eq 'yes' then begin
    xfit=findgen(2*dx+1)-dx
    coef=[0.,tan(theta*torad)]
    yfit=poly(xfit,coef)  ;C0 + C1*x + C2*x^2
    oplot,xfit,yfit

    wset,4
    tv,bytscl(wave_angle,-90,90)

    wset,5
    tv,bytscl(angle_error,0,10)
    endif
  endif

  wait,0.005
end

print,'done'

end
