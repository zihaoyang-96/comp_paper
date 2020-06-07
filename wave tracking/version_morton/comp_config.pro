FUNCTION comp_config

; I/O
root_dir='/Data/'
numberWaveLengths = 3            ;INTEGER, 3 or 5. (3pt or 5 pt data).
filenameExtra=''                 ;Set this to a string that is added to the name of the sav-files and plots, to
                                ;distinguish between different runs.

; Cross-correlation parameters
crossCorr={threshold:0.4,maxNumberIterations:200}


;Interploation of missing frame
;Number of coefficients to use for Maximum entropy method
maxentCoefficientNumber=10

; Sun parameters
radiusSunMm = 695.7 ;solar radius in Mm

; Wave angle calculations
coherence={BoxHalfLength:10,smoothing:15,limit:0.5,minNumberPixels:10}
filter={width:0.0015,centralFrequency:0.0035}


; Basic mask parameters 
upperMaskRadius = 280 ; pixel radius just greater than 1.3 Rsun
lowerMaskRadius = 220 ; pixel radius just greater than 1.02 Rsun

;Space time run parameters
maxTrackLength=25  ;choose path length for phase speed measurement (default 25). Make odd.
numberLagValues=29 ;number of lag values to calculate cross-correlation at. Make odd.

return,DICTIONARY('numberWaveLengths', numberWaveLengths, 'radiusSunMm', radiusSunMm, 'lowerMaskRadius', lowerMaskRadius,$ 
                  'upperMaskRadius', upperMaskRadius, 'crossCorr', crossCorr,'coherence',coherence, $
                   'filter', filter, 'maxentCoefficientNumber',maxentCoefficientNumber, $
                   'root_dir',root_dir,'filenameExtra',filenameExtra,'maxTrackLength',maxTrackLength,'numberLagValues',numberLagValues)


END
