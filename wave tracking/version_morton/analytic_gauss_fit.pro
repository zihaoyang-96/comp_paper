;+
; :Description:

;   procedure to compute analytic fit to a gaussian sampled at three points

;  input:
;   i1, i2, i3:  the intensity at three points in the line profile increasing monotonically in wavelength
;   d_lambda:  the wavelegth spacing of the samples

;  output:
;   doppler_shift: the shift of the fit gaussian from the center wavelength in the same units as d_lambda
;   width:  the linewidth in the same units as d_lambda
;   i_cent: the central intensity of the gaussian in the same units as i1, i2, i3
;
; :Params:
;    profile
;    d_lambda
;    doppler_shift
;    width
;    i_cent
;
;
;
; :Author: Christian Bethge
;-
pro analytic_gauss_fit,profile,d_lambda,doppler_shift,width,i_cent

  if (profile[0] lt 0 or profile[1] lt 0 or profile[2] lt 0) then begin
    
    width = 0D
    doppler_shift = 0D
    i_cent = 0D
    goto, skip_calc
  endif
  
  i1=profile[0]
  i2=profile[1]
  i3=profile[2]
  
  a=alog(i3/i2)
  b=alog(i1/i2)
  
  if (-2D*d_lambda^2D/(a+b)) lt 0 then begin
    
    width = 0D
    doppler_shift = 0D
    i_cent = 0D
    goto, skip_calc
  endif
  
  width=sqrt( -2D*d_lambda^2D/(a+b) )
  doppler_shift=width^2D/(4D*d_lambda)*(a-b)
  i_cent=i2*exp(doppler_shift^2D/width^2D)
  
  skip_calc:
end
