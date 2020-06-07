;PURPOSE: Defines a 1D apodization window
;
;INPUTS: tn - length of time-series
;
;OPTIONAL INPUTS: apod - sets strength of apodisation window, default=0.2 
;
;OUTPUTS: returns apodisation window function
;         cpg- Coherent Power Gain - correction factor needed for amplitude/power of time-series after apodisation
;
;HISTORY: Created by R J Morton - borrowed apod function from R Rutten's plotkowpower
;


FUNCTION do_apod,tn,cpg,apod=apod
 
      apodt = fltarr(tn)+1
      IF n_elements(apod) EQ 0 THEN apod=0.2
      apodrimt = tn*apod
      apodt[0] = (sin(!pi/2.*findgen(apodrimt)/apodrimt))^2
      apodt = apodt*shift(rotate(apodt,2),1)
     
      ;Coherent Power Gain - correction factor needed for power after apodisation 
      CPG=total(apodt)/tn 
RETURN,apodt
END