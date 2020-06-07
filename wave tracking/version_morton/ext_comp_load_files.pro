;
;PURPOSE: Loads suitable CoMP files and interpolates. Called by wave_tracking.pro
;
;OPTIONAL INPUTS: /debug - for debugging
;                 lim  - sets number of missing data files allowed to interpolate over.
;                        Value is set in wave_tracking- default is 2
;                        code cannot deal with more than this yet (probably sensible anyway!) 
;
;OUTPUTS: cube_i - intensity cube
;         cube_v - velocity cube
;         cube_w - doppler width cube
;
;
;CALLS: file_list.pro,headfits.pro,sxpar.pro,readfits.pro
;
;HISTORY: Written by R J Morton 05/2015 using parts of wave_tracking.pro by S. Tomczyk
;         Made data interpolation routine generic - RJM 06/2015
;         Created index & temp_index to keep together useful quantities from wave tracking RJM 04/2016
;         Fixed bug for files with no missing frames - RJM/C. Bethge 06/2016

PRO ext_comp_load_files,inpath=inpath,cube_i,cube_v,cube_w,index,debug=debug,lim=lim,cadence=cadence

rsun_mm = 695.7 ;solar radius in Mm

flist   = find_files('*.fts.gz',inpath)
;if temp_index.start_file ne 0 then begin
; flist = flist[temp_index.start_file:*]
;endif
fnumber = n_elements(flist)
cadence = intarr(fnumber-1)
hdr     = headfits(flist[0],ext=1)
nx      = sxpar(hdr,'NAXIS1')
ny      = sxpar(hdr,'NAXIS2')

index=fitshead2struct(headfits(flist[0]) )


;================================================================
;compute cadence
;================================================================
FOR ii=0,fnumber-2 do begin
 hdr       = headfits(flist[ii])
 date_obs  = sxpar(hdr,'DATE-OBS')
 time_obs  = sxpar(hdr,'TIME-OBS')
 IF n_elements(time) EQ 0 THEN time=time_obs ELSE time=[temporary(time),time_obs]
 new_time  = anytim2tai(date_obs+' '+time_obs) 
 IF ii GT 0 THEN cadence[ii-1] = fix(new_time-old_time)
 old_time = new_time
  
ENDFOR


IF keyword_set(debug) THEN plot,cadence,xtitle='Frame No.',ytitle='Cadence (s)'
norm_cadence = median(cadence)

;================================================================
; look for the start of the normal cadence in case the sequence 
; starts with files with an abnormal cadence
; RJM updated section - MAY2015
;================================================================
seq_start = 0
seq_end   = fnumber-1


;================================================================
; cut out contiguous file list and prepare
; the interpolated data cube
; RJM updated section - MAY2015
;================================================================
print,'Start sequence',seq_start,' End sequence',seq_end

nt = seq_end-seq_start+1

cube_i = fltarr(nx,ny,nt)
cube_v = fltarr(nx,ny,nt)
cube_w = fltarr(nx,ny,nt)


xscale = rsun_mm*index.CDELT1/index.RSUN ;Mm
time_test = strarr(nt)+'ZERO' ; (in)sanity check

;================================================================
; fill cubes with data and interpolate if necessary
; RJM updated section - MAY2015
; Should now be able to interpolate over any number of frame gaps…
;…IF DESIRED - 'Authenticity' decreases with increasing interpolation gaps
;================================================================

;First read in files and put in cube in their place, leaving gaps for missing frames
k=0
add=0
FOR ii=0,nt-1 DO BEGIN
   
    
        cube_i[*,*,ii] = readfits(flist[ii-add],ext=1,/silent)
        cube_v[*,*,ii] = readfits(flist[ii-add],ext=3,/silent)
        cube_w[*,*,ii] = readfits(flist[ii-add],ext=4,/silent)
        time_test[ii]  = sxpar(headfits(flist[ii-add]),'TIME-OBS')

ENDFOR



index=create_struct(temporary(index), 'NAXIS1', nx, 'NAXIS2', ny, 'NFRAMES', nt, $ 
                    'TIME_D$OBS_LIST',time_test[0:nt-1], 'XSCALE', xscale, 'NORM_CADENCE',  norm_cadence)

save,cube_i,cube_v,cube_w,index,filename=inpath+'cubes.sav'

END