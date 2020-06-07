;
;PURPOSE: Loads suitable CoMP files and fills in gaps. Default gap filling by linear interpolation.
;         Called by wave_tracking.pro
;
;OPTIONAL INPUTS: /debug - for debugging
;                 frameLimitInterpolation  - sets number of missing data files allowed to interpolate over.
;                        Value is set in wave_tracking- default is 2, i.e. interpolate over one missing frame
;                 max_ent - use maximum entropy method for gap filling as opposed to linear interpolation
;                        
;
;OUTPUTS: cube_i - intensity cube
;         cube_v - velocity cube
;         cube_w - doppler width cube
;
;
;CALLS: file_list.pro,headfits.pro,sxpar.pro,readfits.pro
;
;HISTORY: - Written by R J Morton 05/2015 using parts of wave_tracking.pro by S. Tomczyk
;         - Made data interpolation routine generic - RJM 06/2015
;         - Created index & temp_index to keep together useful quantities from wave tracking RJM 04/2016
;         - Fixed bug for files with no missing frames - RJM/C. Bethge 06/2016
;         - Edited as part of Python port. Fixed some bugs in file loading - RJM 03/2019
;         - v0.3 (RJM 05/2019) Introducing versioning to keep track of code development. Given past major changes, 
;           the first numbered version is as version 0.3.

PRO comp_load_files,cube_i,cube_v,cube_w,index,config,frameLimitInterpolation=frameLimitInterpolation,max_ent=max_ent,no_hard_mask=no_hard_mask

comm_string='Processed with v0.3 of CoMP WTC on '+systime()+'. '

nwlst = strcompress(string(config['numberWaveLengths']),/rem)
fileList   = file_search('~/'+config['inpath']+'*.comp.1074.dynamics.'+nwlst+'.fts.gz',/expand_tilde)
if config['startFile'] ne 0 then begin
 fileList = fileList[config['startFile']:-1]
endif
numberFiles = n_elements(fileList)
cadence = intarr(numberFiles-1)
hdr     = headfits(fileList[0],ext=1)
nx      = sxpar(hdr,'NAXIS1')
ny      = sxpar(hdr,'NAXIS2')

index=fitshead2struct(headfits(fileList[0]) )


;================================================================
;compute cadence
;================================================================
FOR ii=0,numberFiles-2 do begin
   hdr       = headfits(fileList[ii])
   dateObservations  = sxpar(hdr,'DATE-OBS')
   timeObservations  = sxpar(hdr,'TIME-OBS')
   IF n_elements(time) EQ 0 THEN time=timeObservations ELSE time=[temporary(time),timeObservations]
   new_time  = anytim2tai(dateObservations+' '+timeObservations) 
   IF ii GT 0 THEN cadence[ii-1] = fix(new_time-old_time)
   old_time = new_time
ENDFOR


normalCadence = median(cadence)

;================================================================
; Determine which files to read in based on cadence
; look for the start of the normal cadence in case the sequence 
;================================================================
sequenceStart = 0
;sequenceEnd   = numberFiles-1

in=where(cadence EQ normalCadence)
sequenceStart=in[0]

IF n_elements(frameLimitInterpolation) EQ 0 THEN frameLimitInterpolation=2

;Determine final image in series either last frame or frame before large gap
;Locations where the jump is larger than the frameLimitInterpolationit x normal cadence
largeGapLocation=where( cadence[sequenceStart:numberFiles-2] GT (frameLimitInterpolation)*normalCadence )
IF largeGapLocation[0] EQ -1 THEN sequenceEnd=numberFiles-1 $
                             ELSE sequenceEnd=largeGapLocation[0]+sequenceStart


;Always plot info to output file
maxCadence=max(cadence)<1000
cad_plot=plot(cadence,xtitle='Frame No.',ytitle='Cadence (s)',/buffer,yrange=[0,maxCadence],title='Cadence of entire data set')
cad_plot_end=plot([sequenceEnd,sequenceEnd],[0,maxCadence],'r',line=2,/overplot,name='Series end')
leg=legend(target=[cad_plot_end],/normal,position=[0.8,0.8])
cad_plot.save,config['outpath']+'data_gaps.png'


;Cuts down file list to those values being loaded in
fileList = fileList[sequenceStart:sequenceEnd]
cadence=cadence[sequenceStart:sequenceEnd-1]

;find locations where cadence is greater than the normal cadence
locationOfGaps=where(cadence GT 1.1*normalCadence)
gapStart=locationOfGaps

;Check to see whether any missing frames and if any gaps are too close to end of sequence
IF where(locationOfGaps[0] lt sequenceEnd) EQ -1 THEN numberGaps=0 $
ELSE IF gapStart[0] GE 0 THEN BEGIN 
        endpoint=where(locationOfGaps GT sequenceEnd-frameLimitInterpolation,complement=keepGaps)
        IF endpoint[0] gt 0 THEN BEGIN
            sequenceEnd=locationOfGaps[endpoint[0]]
            numberGaps=n_elements(keepGaps)
            locationOfGaps=temporary(locationOfGaps[keepGaps])
            gapStart=locationOfGaps
        ENDIF
     ENDIF ELSE numberGaps=0 ; number of gaps

print,'======================'
print,'Files to be read in: Start of sequence',sequenceStart,' || End of sequence',sequenceEnd

;================================================================
;Works out which positions need filling
;including the extension of the array for previously 
;interpolated frames
;================================================================
IF numberGaps GT 0 THEN BEGIN
      gapSize=cadence[locationOfGaps]/normalCadence
      cumulativeGapSize=total( [0, gapSize]-1 ,/cum)+1
      gapStart=fix(locationOfGaps+1+cumulativeGapSize)
      numberMissingFrames=cumulativeGapSize[-1]
      print,'Location of missing frames',gapStart
      print,'Size of gap', fix(gapSize)-1
      print,'Total missing frames',numberMissingFrames
      print,'======================'

      gapEnd=gapStart+fix(gapSize)-2

      add=0
      frames_to_skip=intarr(numberMissingFrames)
      FOR i=0,numberGaps-1 DO BEGIN 
          num_skip=indgen(gapSize[i]-1)
          frames_to_skip[i+add:i+add+num_skip[-1]]=gapStart[i]+num_skip
          add+=num_skip[-1]
      ENDFOR
ENDIF ELSE numberMissingFrames=0

;================================================================
; fill cubes with data and interpolate if necessary
; Should now be able to interpolate over any number of frame gaps…
;…IF DESIRED - 'Authenticity' decreases with increasing interpolation gaps
;================================================================


nt = sequenceEnd-sequenceStart+1+numberMissingFrames
; manually exclude bad data
IF config['numberFilesToProcess'] NE 0 THEN BEGIN
  IF config['numberFilesToProcess'] GT nt THEN config['numberFilesToProcess']=nt $
  ELSE nt=config['numberFilesToProcess']
ENDIF



cube_i = fltarr(nx,ny,nt)
cube_v = fltarr(nx,ny,nt)
cube_w = fltarr(nx,ny,nt)


xscale = config['radiusSunMm']*index.CDELT1/index.RSUN ;Mm
time_test = strarr(nt)+'ZERO' ; (in)sanity check

;First read in files and put in cube in their place, leaving gaps for missing frames
k=0
add=0
print,"Reading in data"
FOR ii=0,nt-1 DO BEGIN
   
    IF ii NE (frames_to_skip[k]) THEN BEGIN
        cube_i[*,*,ii] = readfits(fileList[ii-add],ext=1,/silent)
        cube_v[*,*,ii] = readfits(fileList[ii-add],ext=3,/silent)
        cube_w[*,*,ii] = readfits(fileList[ii-add],ext=4,/silent)
        time_test[ii]  = sxpar(headfits(fileList[ii-add]),'TIME-OBS')

    ENDIF ELSE add=add+1
    
    If numberMissingFrames GT 0 THEN IF k LT numberMissingFrames-1 THEN IF ii EQ frames_to_skip[k] THEN k=k+1
     counter,ii+1,nt,/percent 
ENDFOR


; Create mask - mask out pixels not to process
mask=intarr(nx,ny)
x=rebin(findgen(ny)-(nx/2.)+0.5,nx,ny)
y=transpose(x)
rPolar=sqrt(x^2+y^2)

av=median(cube_i,dim=3)
good=where(rPolar ge config['lowerMaskRadius']+config['maxOcculterOffset'] and rPolar le config['upperMaskRadius'] and median(av,5) gt 0 and av gt 0)
mask[good]=1


IF NOT keyword_set(no_hard_mask) THEN BEGIN
  print,'Using hard mask'
  FOR i=0,nx-1 DO FOR j=0,ny-1 DO BEGIN
      IF mask[i,j] eq 1 THEN BEGIN
         zeros=where(cube_i[i,j,*] eq 0)

         IF n_elements(zeros) gt numberMissingFrames THEN mask[i,j]=0
      ENDIF
  ENDFOR
ENDIF

cube_i*=rebin(mask,nx,ny,nt)
cube_v*=rebin(mask,nx,ny,nt)
cube_w*=rebin(mask,nx,ny,nt)


;Fill in missing frames
IF numberMissingFrames GT 0 THEN BEGIN
  IF NOT keyword_set(max_ent) THEN BEGIN
          ii=0
          nx_arr=findgen(nx)
          ny_arr=findgen(ny)
          IF numberGaps GT 0 THEN BEGIN
             print,"Interpolating data"
              WHILE ii LT numberGaps DO BEGIN
                   gapStartIndex=gapStart[ii]
                   gapEndIndex=gapSize[ii]
                   z=findgen(gapEndIndex)/gapEndIndex
                   z=z[1:gapEndIndex-1]
                   ;print,gapStartIndex,gapEndIndex,ii
                   dumi=[[[cube_i[0:-1,0:-1,gapStartIndex-1]]],[[cube_i[0:-1,0:-1,gapStartIndex+gapEndIndex-1]]]]
                   dumv=[[[cube_v[0:-1,0:-1,gapStartIndex-1]]],[[cube_v[0:-1,0:-1,gapStartIndex+gapEndIndex-1]]]]
                   dumw=[[[cube_w[0:-1,0:-1,gapStartIndex-1]]],[[cube_w[0:-1,0:-1,gapStartIndex+gapEndIndex-1]]]]
                   cube_i[0:-1,0:-1,gapStartIndex:gapStartIndex+gapEndIndex-2] = interpolate(dumi,nx_arr,ny_arr,z,/grid)
                   cube_v[0:-1,0:-1,gapStartIndex:gapStartIndex+gapEndIndex-2] = interpolate(dumv,nx_arr,ny_arr,z,/grid)
                   cube_w[0:-1,0:-1,gapStartIndex:gapStartIndex+gapEndIndex-2] = interpolate(dumw,nx_arr,ny_arr,z,/grid)
                   ;ii+=gapEndIndex-1
                   ii+=1
              ENDWHILE
            
          ENDIF
          comm_string=comm_string+' Gaps filled with linear interpolation. '
  ENDIF ELSE BEGIN

          print,'Interpolating data via Maximum Entropy'
          numberGaps=n_elements(gapStart)
          type_arr=fltarr(nx,ny,numberGaps,3)
          coef_num=config['maxentCoefficientNumber']
          FOR i=0,nx-1 DO FOR j=0,ny-1 DO BEGIN
              
              IF mask[i,j] EQ 1 THEN BEGIN
                 ;cube_i --------------
                 temp=reform(cube_i[i,j,*])
                 fill_gaps,temp,gapStart,gapEnd,numberGaps,coef_num=coef_num,types=types
                 cube_i[i,j,*]=temp

                 ;Saves type of 'filling' used for each gap
                 IF (where(types eq 0))[0] NE -1 THEN types[where(types eq 0)]=6
                 type_arr[i,j,0:-1,0]=types

                 ;cube_v --------------
                 temp=reform(cube_v[i,j,*])
                 fill_gaps,temp,gapStart,gapEnd,numberGaps,coef_num=coef_num,types=types
                 cube_v[i,j,*]=temp
                 IF (where(types eq 0))[0] NE -1 THEN types[where(types eq 0)]=6
                 type_arr[i,j,0:-1,1]=types

                 ;cube_w --------------
                 temp=reform(cube_w[i,j,*])
                 fill_gaps,temp,gapStart,gapEnd,numberGaps,coef_num=coef_num,types=types
                 cube_w[i,j,*]=temp
                 ;Saves type of 'filling' used for each gap
                 IF (where(types eq 0))[0] NE -1 THEN types[where(types eq 0)]=6
                 type_arr[i,j,0:-1,2]=types

              ENDIF
              counter,j+i*ny,nx*ny,/percent 
          ENDFOR
          comm_string=comm_string+' Gaps filled with Maximum Entropy method. '

  ENDELSE
ENDIF




index=create_struct(temporary(index), 'NAXIS1', nx, 'NAXIS2', ny, 'NFRAMES', nt, $ 
                    'TIME_D$OBS_LIST',time_test, 'XSCALE', xscale, 'normalCadence',  normalCadence, $
                    'numberWaveLengths', config['numberWaveLengths'], 'START_FILE', config['startFile'], 'lowerMaskRadius', config['lowerMaskRadius'],$
                    'upperMaskRadius', config['upperMaskRadius'], 'maxOcculterOffset', config['maxOcculterOffset'], $
                    'coordinatesCCBox', config['coordinatesCCBox'], 'rmsOfCrossCorr', fltarr(2),'missingFrames', gapStart,$
                    'gapSize', gapSize-1, 'mask', mask, 'COMMENTS', comm_string)



END