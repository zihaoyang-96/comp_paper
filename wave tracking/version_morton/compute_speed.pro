
;PURPOSE: Subroutine to compute phase speed by cross correlating rows of the velocity
;         space-time diagram (which is input). The time series used for cross correlation 
;         is formed by collapsing all the rows of the diagram after shifting by an 
;         estimate of the phase speed. This estimate is obtained using only the central 
;         rows of the diagram.
;  
;
;
;
;  INPUTS: vel - filtered time-distance velocity map. Either pro or retrograde.
;          delx - pixel size in Mm
;          delt - cadence of observations
;
;  KEYWORDS:   debug - obvious...
;
;              ret   - set this if the computed speed is a retrograde
;                      speed (for a better initial guess for the fit).
;         
;              force_zero - force the fit to go through (0,0), i.e.
;                           through the reference pixel. By default, it
;                           it is turned off in the wave tracking code,
;                           because it gives higher errors even though
;                           the fit might still be good (i.e. same
;                           slope, but with an offset in x/y). 
;
;CALLS: parabola.pro, mpfit.pro, mpfitfunc.pro, lin_func_fit.pro
;
;
;
;HISTORY: - Created by S Tomczyk
;         - Corrected in Feb 2013 by C. Bethge. More accurate phase speeds and
;           errors now.
;         - Added code to deal with even length tracks.
;           Added extra constraint for CC - max CC value has to be
;           above probability of that of random signals.
;           Corrected bugs in phase speed error calculation
;           Added new method for rejection of outliers
;           Added correction to t-factor for phase speed error to account for varying track lengths - R Morton 11/2014
;           Fixed bugs in errors for phase speed fitting -  R Morton 05/2015
;         - v0.3 (RJM 05/2019) Introducing versioning to keep track of code development. Given past major changes, 
;           the first numbered version is as version 0.3.


FUNCTION compute_speed, vel, delx, delt, debug=debug, ret=ret, force_zero=force_zero


s=size(vel)
nt=s(1)
npt=s(2)
IF npt mod 2 EQ 0 THEN npt=npt-1 ;deals with even length tracks
icent=npt/2


nlag=npt  ;number of points in cross correlation (make odd)
lag=indgen(nlag)-fix(nlag/2.)
n = 19 < npt  ;number of points along track for initial guess (make odd)
cent=fltarr(n)



;##############################################################
;Cross-correlate adjacent rows in time-distance map
;##############################################################

FOR i=-n/2,n/2 DO BEGIN
  ccor=c_correlate(vel(*,icent),vel(*,icent+i),lag)
  m=max(ccor,imax)
  imax=1 > imax < (nlag-2)
  cent(i+n/2)=parabola(lag(imax-1:imax+1),1.-ccor(imax-1:imax+1))
;  cent(i+n/2)=total(lag(imax-1:imax+1)*ccor(imax-1:imax+1))/total(ccor(imax-1:imax+1))
ENDFOR


;TESTING USING FFT
;cent_test=fltarr(n,2)
;FOR i=-n/2,n/2 DO BEGIN
;  TT=secondary_compute_speed(vel(*,icent),vel(*,icent+i),f)
;  cent_test(i+n/2,*)=TT
;ENDFOR


;##############################################################
;Calculate derivative of the lags giving 1/c_ph & take median value
;##############################################################
c=[0.,median(deriv(cent))]   ;use robust estimate for mean shift;  collapse all time series along track
col=fltarr(nt)


;##############################################################
;Assuming phase speed is constant along track, calculate the shift
;between each row and collapse to produce 'super' row, i.e.
;create an improved noise time-series with which to correlate
;##############################################################
FOR i=0,npt-1 DO BEGIN
  shift=c(1)*(i-npt/2)+c(0)
  col=col+interpolate(vel(*,i),findgen(nt)+shift)
ENDFOR
col=col/float(npt)


;##############################################################
;Use collapsed time series to perform cross correlation with all time series along track
;##############################################################
cent=fltarr(npt)
dcent=fltarr(npt)+.1

IF keyword_set(debug) THEN BEGIN
 plot,lag,/nodata,xr=[-nlag/2.,nlag/2.],yr=[-npt/4.-1,npt/4.+1],ysty=1,xsty=1,xtit='Lag (Time Steps)',$
 	ytit='Correlation',chars=1.5
ENDIF


;##############################################################
;First iteration: 
;For an initial guess to determine a useful
;range for searching the position of max correlation
;##############################################################
FOR i=0,npt-1 DO BEGIN
  ccor=c_correlate(col,vel(*,i),lag)
  IF keyword_set(debug) THEN BEGIN
   oplot,lag,ccor-npt/4.+i*0.5, col=50
  ENDIF
  m=max(ccor,imax)
  imax=1 > imax < (nlag-2)
  cent(i)=parabola(lag(imax-1:imax+1),1.-ccor(imax-1:imax+1))
  
ENDFOR

;cent_init2=cent
;FOR i=0,npt-1 DO BEGIN
;  TT=secondary_compute_speed(col_test,vel(*,i))
;   cent_test(i,*)=TT
;ENDFOR


;##############################################################
;Previous round of cc provides new cc positions wrt to the central
;pixel.
;These central values are used to provide an estimate (init_guess) for the shift
;which is assumed constant.
;The initial guess is compared to the cc calculated shifts and the
;residuals are found (pdistance)
;##############################################################

shifts=deriv(cent)
init_guess=lag*median(shifts)

pdistance = fltarr(npt)
FOR ii=0,npt-1 DO BEGIN
 pdistance[ii] = pnt_line([lag[ii],cent[ii]], [lag[0],init_guess[0]], [lag[npt-1],init_guess[npt-1]])
ENDFOR


;##############################################################
;Residuals used to select good pixels for fitting
;##############################################################
good_p = where(pdistance le 2.*median(pdistance))
err = 1.


;##############################################################
;Fit line with constraints on fit
;The y intercept is constrained to go through 0 (.fixed)
;The gradient is confined to be between 0.001 and 20 (.limited & .limits)
;##############################################################
parinfo = replicate({value:0.D, fixed:0, limited:[0,0], limits:[0.D,0.D]}, 2)
parinfo[1].limited[0] = 1 
parinfo[1].limited[1] = 1 
parinfo[0].fixed = 1

if keyword_set(ret) then begin
 parinfo[1].limits[0]  = -20D
 parinfo[1].limits[1]  = -0.001D
 parinfo[*].value      = [0D,-0.5D] 
endif else begin
 parinfo[1].limits[0]  = 0.001D
 parinfo[1].limits[1]  = 20D
 parinfo[*].value      = [0D,0.5D] 
endelse

param = mpfitfun('lin_func_fit', lag[good_p], cent[good_p], err, parinfo=parinfo, /quiet,xtol=1e-4)
i_yfit = param[0] + lag*param[1]

;##############################################################
;Fit results used to define bounds for which to fit cross-correlation
;curve with parabola to find maximum
;Debug option plots the bounds (triangles), the cc max (diamonds)
;and the maximum from the parabola fit (triangle)
;##############################################################
cor_val=fltarr(npt)
tstep_tolerance = 4
for i=0,npt-1 do begin
  ccor=c_correlate(col,vel(*,i),lag)
  lo_l = i_yfit[i]-tstep_tolerance+npt/2. ; lower limit
  up_l = i_yfit[i]+tstep_tolerance+npt/2. ; upper limit
  if (lo_l gt 1 and up_l lt nlag-2) then begin 
    m=max(ccor[lo_l:up_l],imax)
    imax=imax+i_yfit[i]-tstep_tolerance+npt/2.
    if keyword_set(debug) then plots, imax-nlag/2., m -npt/4. + i*0.5, psym = 4, color=100
    if keyword_set(debug) then plots, lo_l-nlag/2., ccor[lo_l] -npt/4. + i*0.5, psym = 5, color=250
    if keyword_set(debug) then plots, up_l-nlag/2., ccor[up_l] -npt/4. + i*0.5, psym = 5, color=250

  endif else begin
   m=max(ccor,imax)
   imax=1 > imax < (nlag-2)
   if keyword_set(debug) then plots, imax-nlag/2., m -npt/4. + i*0.5, psym = 4, color=185
  endelse
  cent(i)=parabola(lag(imax-1:imax+1),1.-ccor(imax-1:imax+1))
  cor_val(i)=m
if keyword_set(debug) then plots, cent[i], m -npt/4. + i*0.5,psym = 5,color=200
endfor



data = cent*delt

xfit=(findgen(npt)-fix(npt/2))*delx
if keyword_set(debug) then begin
 plot,xfit,data,psym=1,xtit='Position Along Track (Mm)',ytit='Lag (s)',chars=1.5,$
     xr=[-npt/2,npt/2]*delx*1.5,xsty=1,yr=[-npt,npt]*4,ysty=1
;********

endif

;##############################################################
;Do the fitting again to find new constrained slope
;##############################################################
shifts=deriv(cent)*delt/delx
init_guess=xfit*median(shifts)

;first iteration
pdistance = fltarr(npt)
for ii=0,npt-1 do begin
 pdistance[ii] = pnt_line([xfit[ii],data[ii]], [xfit[0],init_guess[0]], [xfit[npt-1],init_guess[npt-1]])
endfor

good_p = where((pdistance le 2.*median(pdistance)) and (cor_val gt 0.2));2./sqrt(nt)))
err = 1.


parinfo = replicate({value:0.D, fixed:0, limited:[0,0], limits:[0.D,0.D]}, 2)
parinfo[1].limited[0] = 1 
parinfo[1].limited[1] = 1 
if keyword_set(force_zero) then begin
  parinfo[0].fixed = 1
endif

if keyword_set(ret) then begin
 parinfo[1].limits[0]  = -20D
 parinfo[1].limits[1]  = -0.001D
 parinfo[*].value      = [0D,-0.5D] 
endif else begin
 parinfo[1].limits[0]  = 0.001D
 parinfo[1].limits[1]  = 20D
 parinfo[*].value      = [0D,0.5D] 
endelse

param = mpfitfun('lin_func_fit', xfit[good_p], data[good_p], err, parinfo=parinfo, /quiet,xtol=1e-4)
i_yfit = param[0] + xfit*param[1]

;##############################################################
;second iteration
;
;RJM - attempt to make rejections statistically significant 
;##############################################################

pdistance = fltarr(npt)
;for ii=0,npt-1 do begin
; pdistance[ii] = pnt_line([xfit[ii],data[ii]], [xfit[0],i_yfit[0]], [xfit[npt-1],i_yfit[npt-1]])
;endfor

;pnt_line only gives perpendicular distance from line - what we really want is residuals
pdistance=i_yfit-data

;More statistically significant tests for data points
;First assess the residuals:
;- remove values that don't meet correlation test significance (see below) 
;- calculate standard deviation for residuals
;- remove data points who residuals suggest that data points are 'outliers' (or erroneous)
;
;The selection of 'outlies' requires a cutoff value that corresponds to a probability of occurrence. 
;The probability is a user defined limit.
;Assuming we want a 99.7% cutoff (3sigma for a normal distribution), i.e., points removed with a 0.3 percent chance of occurrence.
;We assume data points are sampled from a normal distribution, hence statistics follow t-distribution
;The cutoff value is then given by standard deviation*t_cvf(0.3/2.,v), where t_cvf calculates the
;t values for a two tailed test (hence the divide by 2) for v degrees of freedom. v is n-1 where
;n is the number of samples
;
;
;Correlation value has to be greater than that of random series. Cross correlation
;of random series should have sample correlations of mean 0 and variance 1/N
;N is sample number
;Using Normal CDF - erf^-1( p) *sqrt(2)/sqrt(N) you obtain a cutoff value for cross correlation
;A cutoff of 2.5sigma is p=0.987 giving significance level of ~2.58/sqrt(N)
;A cutoff of 3sigma is p=0.997 giving significance level of ~2.99997/sqrt(N)
;A cutoff of 4sigma is p=0.9999 giving significance level of ~3.99996/sqrt(N)
;A cutoff of 5 sigma is p=0.999 999 giving significance level of ~4.9999/sqrt(N)
;

;good_p = where((pdistance le 10.*median(pdistance)) and (cor_val gt 0.2),comp=bad_p)
good_pi=where(cor_val gt 0.2)
ngp=n_elements(good_pi)
IF ngp GT 2 THEN BEGIN
    sigma=sqrt( (moment(pdistance[good_pi]))[1] )
    good_p = where( (abs(pdistance) le t_cvf(0.03/2.,ngp-1.)*sigma ) and (cor_val gt 2.9997/sqrt(nt)),comp=bad_p)
ENDIF ELSE good_p=0

num_good_p = n_elements(good_p)

;Can only fit a straight line with 3 points!
IF num_good_p GT 2 THEN BEGIN

	;Count number of finite values in pdistance
	fin_val = where(finite(pdistance) eq 1,num_fin_val)
	if num_fin_val eq n_elements(pdistance) then begin
 		err = abs(pdistance[good_p]) ;added absolute RJM 29/05/2015
	endif else begin
 		err=fltarr(n_elements(pdistance[good_p]))
 		err[*] = 1.
	endelse

	;In some cases the error is set to zero and mpfit is not impressed
	;added to avoid this
       ;tweaked so that error takes minimum value of non-zero errors rather than max as before  RJM 29/05/2015
	in=where(err eq 0,comp=val)
	err[in]=min(err[val])

	param = mpfitfun('lin_func_fit', xfit[good_p], data[good_p], err, parinfo=parinfo, /quiet,xtol=1e-4)
	new_yfit = param[0] + xfit*param[1]


	;##############################################################
	;DERVIVATION OF ERROR ON FIT
	;From Bevington & Robinson (2003) - Chapter 6 107-111
	;Estimate for variance of data and error in gradient
	;Assume Chi^2_v is 1 so no information on quality of fit 
	;
	;Bevington & Robinson miss correction factor for sampled data
	;See Acton - Analysis of straight line data (1966)
	;
	; ~68% confidence level is given by
	;t_factor given by t_cvf((1.-0.683)/2.,num_good_p-1)
	;Different t_factors for different confidence levels
	; speed_error =((sqrt(total((data[good_p]-new_yfit[good_p])^2.)/(num_good_p-2)))/
	;(sqrt(s_xx)))*t_factor 
	;##############################################################


	if (num_good_p gt 0.3333*num_fin_val) and (num_good_p gt 2) then begin
		 phase_speed = 1./param[1] 
		 x_m = total(xfit[good_p])^2./float(num_good_p)
		 s_xx = (total(xfit[good_p]^2.)-x_m)
 		  t_factor =t_cvf((1.-0.683)/2.,num_good_p-1) ;1.02 - corrected RJM

 		;Correction to phase speed error - not the same as error in gradient -RJM
 		grad_error = ( (sqrt( total( (data[good_p]-new_yfit[good_p])^2. )/ (float(num_good_p)-2.) ) )/(sqrt(s_xx)))*t_factor
 		speed_error=grad_error/param[1]^2
	endif else begin
 
 		phase_speed = !Values.F_NAN
 		speed_error = !Values.F_NAN
	endelse



	if keyword_set(debug) then begin
	 if n_elements(bad_p) eq 1 then begin 
 		 if bad_p[0] ne -1 then oplot,xfit[bad_p],data[bad_p], color=230, psym=2
	       
 	  endif else oplot,xfit[bad_p],data[bad_p], color=230, psym=2

         oplot,xfit,new_yfit, color = 150,linestyle=2
         oplot,xfit,init_guess, color = 185
	endif
ENDIF ELSE BEGIN

   phase_speed = !Values.F_NAN
   speed_error = !Values.F_NAN
ENDELSE

return,[phase_speed,speed_error]


end
