function min_dist_fit,x,y,error=error, slope = slope, status = status

;  function fits a line to a set of points under the constraint that
;  1) the line goes through the point x=0, y=0, and
;  2) the total perpindicular distance from the set of points to the line is minimized.

;  This function returns the angle (degrees) of the line. The error on the angle (degrees) is
;  returned in the keyword error.

torad=!pi/180.
slope = 0.
;  compute sums

sxsq=total(x^2)
sysq=total(y^2)
sxy=total(x*y)

if sxy ne 0. then begin

; compute two possible choices for slope

  arg1=(-sxsq+sysq)/(2.*sxy)
  arg2=sqrt( (sxsq-sysq)^2 +4.*sxy^2 )/(2.*sxy)
  m1=arg1-arg2
  m2=arg1+arg2

;  find sum of distance squared for two slopes and choose slope which minimizes sum of distance squared

  distsq1=(m1^2*sxsq  -2.*m1*sxy +sysq)/(1.+m1^2)
  distsq2=(m2^2*sxsq  -2.*m2*sxy +sysq)/(1.+m2^2)

  if distsq1 le distsq2 then m=m1 else m=m2
  ang=atan(m)
  ss = [m1, m2]
  slope = min(ss) + (max(ss)-min(ss))/2.

endif else if sxsq ge sysq then ang=0. else ang=!pi/2.

;  find error

if sxy ne 0. or sxsq ne sysq then begin

  dang=2.*torad     ; 2 degree perturbation in radians
  ang_plus=ang+dang
  ang_minus=ang-dang

  m=tan(ang)
  mplus=tan(ang_plus)
  mminus=tan(ang_minus)

  d0=(m^2*sxsq  -2.*m*sxy +sysq)/(1.+m^2)
  dplus=(mplus^2*sxsq  -2.*mplus*sxy +sysq)/(1.+mplus^2)
  dminus=(mminus^2*sxsq  -2.*mminus*sxy +sysq)/(1.+mminus^2)

  error=1./(sqrt( (dplus+dminus-2.*d0)/dang^2 ))

endif else error=0.

ang=ang/torad   ;convert to degrees
error=error/torad

status = 0
if finite(ang) and finite(error) and finite(slope) then status = 1

return,ang
end