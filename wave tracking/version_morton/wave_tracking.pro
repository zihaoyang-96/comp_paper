;  PURPOSE: Procedure to determine the wave propagation angle from the cross
;  correlation of the velocity time series of a test point with
;  adjacent time series. The computation is performed in the Fourier
;  domain. Only contiguous points having coherence greater than a
;  threshold (typically 0.5) are used. SSWIDL is needed.
;
;  The default is to use CoMP L2 data, in which case you just need to
;  provide a date (and the correct inpath/outpath); the code (ideally)
;  does the rest. You can also provide your own data (see keyword
;  "own_data"). The code starts looking for data early in the day,
;  which is usually the best due to the better seeing conditions. The
;  code then looks for the longest contiguous time sequence, and
;  interpolates over data gaps of one missing file, which is quite
;  often the case. The data is filtered, the default is a central
;  filter position of 3.5 mHz and filter width of 1.5 mHz. Change the
;  variables freq0 and fwidth if you want to change the filtering. 
;  Also, check if the default values of lower_r and upper_r (lower and
;  upper radius of the field-of-view for masking) are suitable for
;  your data.
;
; NOTE: This version will only work with IDL v8.0+
;  
;  INPUT:      date            - STRING in the format yyyymmdd. 
;                                Will work on HAO's machines or if you
;                                provide CoMP L2 data in your own folder
;                                (change inpath in this case). 
;  
;  KEYWORDS:   
;              Keyword related to reading/writing in data
;              ------------------------------------------
;              init_load       - use this keyword for initial data load and gap filling.
;                                This keyword is required for first load of data. If not provided,
;                                default is to look for existing, gap-filled data (saved as cube_ivw_date).
;
;              frameLimitInterp - sets maximum allowed number of missing frames (gaps) in a row, default is 2.
;                                 Default gap filling method is linear interpolation.
;             
;              max_ent         - Uses maximum entropy method for gap filling
;
;              no_hard_mask    - turns off the hard masking of the data. Hard mask removes any pixel that has at least one frame
;                                with no signal
;
;              Keyword related to cross-correlation of data
;              ------------------------------------------
;
;              cross_corr      - Remove jitter from the CoMP data. Default operation is
;                                to use coordinate locations in find_man_date. If
;                                you do this, make sure to set the debug
;                                keyword for the first run, which will show
;                                you the box for which the cross-correlation 
;                                is computed. If the default box does not look 
;                                correct, set the values x1, x2, y1, y2
;                                manually in find_man_date.  
;
;              choose_corr_box - Select the cross-correlation box with
;                                the mouse.
;
;              dosob           - uses a sobel filtered image for alignment (typically better & faster results)
;
;              Keyword related to wave angle calculation
;              ------------------------------------------
;
;              compute_waveang - will compute the angles of wave propagation, i.e. orientation of magnetic field.
;                                Best to align data before calculating wave angle.
;
;              compute_speeds  - compute phase speeds after the
;                                computation of the wave propagation
;                                angle. 
;
;              plot_maps       - plot ps-files with all quantities after the
;                                computation.
;                              
             
;
;              save_coh        - Create a sav-file containing a measure of
;                                the coherence for each pixel. This is done
;                                by fitting an ellipse to the coherence
;                                island. The coherence for a pixel is good 
;                                when the ellipse is long and slim, i.e.
;                                the normalized ratio of the major and minor
;                                axis of the ellipse is a measure of the
;                                coherence. The sav-file will contain an
;                                array coh_measure, which is just the major
;                                (coh_measure[*,*,0]) and minor
;                                (coh_measure[*,*,1]) axis of the ellipse. 
;                                One possible way to compute a quantity for 
;                                the coherence quality would then be
;  coherence_quality = (coh_measure[*,*,0]-coh_measure[*,*,1])/(coh_measure[*,*,0]+coh_measure[*,*,1])
;
;             debug            - set this if you want to visually
;                                check what's going on. *MUCH MUCH*
;                                slower computation though, so don't
;                                use this as the default.
;                           
;             wr_speeds        - This is mainly intended for checking the
;                                computed speeds afterwards. It will write
;                                out *very* large arrays, so don't
;                                do this on a computer with insufficient
;                                memory or disk space, especially in
;                                combination with the save_coh keyword. Not
;                                really needed for everyday use, mostly for
;                                debugging. 
;
;              altangle - Use alternative estimate for wave angle based on 
;                         sixlin & choosing 
;              splt_cube - save cubes as separate files
;
;CALLS: find_man_date.pro, do_apod.pro, comp_load_files.pro, wave_angle_calc.pro, compute_speed.pro, dejitter.pro,
;        plot_wave_maps.pro
;
;
;
;
;
;HISTORY - written by S. Tomczyk & C. Bethge
;        - Modularised by R J Morton 2015
;        - Made data interpolation routine generic. Reduced computation time for wave_angle_calculation. 
;        - Added apodisation window for wave_angle_calculation to aid in minimising leakage. RJM 06/2015
;        - Added option to use Sobel filtered data for alignment - RJM 11/2015
;        - Created temp_index & index to keep together useful quantities from wave tracking - RJM 04/2016
;        - Redefinition of mask values using coherence measure from wave angles - RJM + CB 06/2016
;        - Changed COMMON blocks for dictionaries to prepare for Python port - RJM 03/2019
;        - v0.3 (RJM 05/2019) Introducing versioning to keep track of code development. Given past major changes, 
;          the first numbered version is as version 0.3.
;          This version also has a number of changes to basic functionality and feedback to user.  

function stop_and_save,cube_i,cube_v,cube_w,index,config,splt_cube=splt_cube
;A small function for saving
   print,'Stop to check data'
   print,"Type '.c' to continue or 'return,0' to quit"
   stop

   print,'Saving data'
   IF keyword_set(splt_cube) THEN BEGIN
      save, cube_i, index, filename=config['outpath']+'cube_i_'+config['date']+config['filenameExtra']+'.sav'
      save, cube_v, index, filename=config['outpath']+'cube_v_'+config['date']+config['filenameExtra']+'.sav'
      save, cube_w, index, filename=config['outpath']+'cube_w_'+config['date']+config['filenameExtra']+'.sav'
   ENDIF ELSE $
      save, cube_i, cube_v, cube_w, index, filename=config['outpath']+'cube_ivw_'+config['date']+config['filenameExtra']+'.sav'
  return,1
END


pro wave_tracking, date, init_load=init_load, frameLimitInterp=frameLimitInterp, max_ent=max_ent,$   ;related to reading/writing in data
                         splt_cube=splt_cube, no_hard_mask=no_hard_mask,$
                         cross_corr=cross_corr, $                          ;related to cross-correlation of data
                         choose_corr_box=choose_corr_box,dosob=dosob, $                  
                         compute_waveang=compute_waveang, altangle=altangle, $   ;related to wave angle calculation   
                         compute_speeds=compute_speeds, $ ;related to wave propagation speed calculation
                         plot_maps=plot_maps, save_coh = save_coh, $       ;additional things to save/plot
                         wr_speeds = wr_speeds,debug = debug              ;related to debugging
                         
                         

;Replace common blocks with dictionaries
config=comp_config()
;Enable easy comparison between IDL and Python versions
;wt_var=DICTIONARY('NWL', 0, 'START_FILE', 0, 'LOWER_R', 0, 'UPPER_R', 0, $
;                      'MAXOFFSET', 0.0, 'CC_COORD', fltarr(2,2), 'RMS_CC', fltarr(2),'DATE',date,'NAME_ADDON','') 
;Empty dictionary - to be filled and passed between modules
;wt_var_ang=DICTIONARY()


; paths where the CoMP data is located and where the .sav-files should be written to
; This will need editing by individual user
; config.inpath  = config['root_dir']+'CoMP/'+strmid(date,8,4,/reverse)+'/'+strmid(date,3,2,/reverse)+'/'+strmid(date,1,2,/reverse)+'/'
; config.outpath = config['root_dir']+'CoMP/wave_tracking_output/'+date+'/'
; config.date = date
config.inpath  = '/Data/CoMP/'+strmid(date,8,4,/reverse)+'/'+strmid(date,3,2,/reverse)+'/'+strmid(date,1,2,/reverse)+'/'
config.outpath = '/Data/CoMP/wave_tracking_output/'+date+'/'
config.date = date

;================================================================
;*** Here the details of the wave tracking need to be defined ***
;================================================================

; ; use 3 pt data by default
; if n_elements(nwl) eq 0 then config['numberWaveLengths'] = 3 $
; ELSE config['numberWaveLengths']=nwl


if keyword_set(wr_speeds) then wr_speeds = 1 else wr_speeds = 0


if file_test(config['outpath'],/DIR) eq 0 then spawn, 'mkdir '+config['outpath']
if keyword_set(plot_maps) then begin
 if file_test(config['outpath']+'plots/',/DIR) eq 0 then spawn, 'mkdir '+config['outpath']+'plots/'
endif


case (config['filter'].centralFrequency) of
 0.0015: freqst = '1.5mHz'
 0.0035: freqst = '3.5mHz'
 0.0055: freqst = '5.5mHz'
 else: begin
         freqst = 'manual_freq_setting'
       end
endcase
config['freqst']=freqst

;================================================================
;***      End of defining the wave tracking parameters        ***
;================================================================


;===============================================
;Read in files
;===============================================

find_man_date,date,config
;For initial load of fits data
IF keyword_set(init_load) THEN BEGIN
 
      ;Test file exists, if not return to user
      ft_res=file_test(config['inpath'])
      IF ft_res EQ 0 THEN BEGIN
        print,'Directory does not exist: '+config['inpath']
        print, 'Check date.'
        return
      END

      
      comp_load_files,cube_i,cube_v,cube_w,index,config,frameLimitInterp=frameLimitInterp,$
                      max_ent=max_ent,no_hard_mask=no_hard_mask


      message, /cont, 'Total # of frames: '+strcompress(string(index.NFRAMES,format='(F5.0)'),/rem)

      ;First save point
      res_ss=stop_and_save(cube_i,cube_v,cube_w,index,config,splt_cube=splt_cube)
      IF res_ss EQ 0 THEN return

ENDIF ELSE BEGIN
      ;For loading previous .sav files

       ;Test file exists, if not return to user
      ft_res=file_test(config['outpath']+'cube_ivw_'+date+config['filenameExtra']+'.sav')
      IF ft_res EQ 0 THEN BEGIN
        print,'File does not exist: cube_ivw_'+date+config['filenameExtra']+'.sav'
        print, 'Check date.'
        print, 'If date correct then run with /init_load keyword.' 
        return
      ENDIF
      print,'Restoring pre-loaded data: cube_ivw_'+date+config['filenameExtra']+'.sav'
      restore, config['outpath']+'cube_ivw_'+date+config['filenameExtra']+'.sav'
ENDELSE


;================================================================
; cross-correlate the data and align it
;================================================================
if keyword_set(cross_corr) then begin
  if keyword_set(debug) then ccdebug = 1

  if keyword_set(choose_corr_box) then begin
     IF keyword_set(dosob) THEN im=sobel(cube_i[*,*,0])<20 ELSE im=cube_i[*,*,0] >0 <25
     ONCE_MORE:
     x2 = -1
     y2 = -1
     window, /free, xs=620, ys=620
     win_nr = !d.window
     aia_lct, wave=193, /load
     tvscl, reform(im)
     message, /cont, 'Choose bottom left corner'
     cursor, x1, y1, /down, /device    
     message, /cont, 'Choose top right corner'

     IF !MOUSE.BUTTON EQ 1 THEN BEGIN 
        CURSOR, do_nothing_x, do_nothing_y, /UP, /WAIT
        !MOUSE.BUTTON = 0  
     ENDIF
     while !Mouse.Button NE 1 do begin
        cursor, tmp_x2, tmp_y2, /change, /device
        tvscl, reform(im)
        plots, [x1,tmp_x2], [y1,y1], /device, color=255 
        plots, [x1,tmp_x2], [tmp_y2,tmp_y2], /device, color=255  
        plots, [x1,x1], [y1,tmp_y2], /device, color=255 
        plots, [tmp_x2,tmp_x2], [y1,tmp_y2], /device, color=255 
     endwhile

     wdelete, win_nr
     if x1 EQ tmp_x2 OR y2 EQ tmp_y2 then begin
       message, /cont, 'Box must be at least 1x1 pixel' 
       goto, once_more
     endif
     if x1 GT tmp_x2 then begin
       x2 = x1
       x1 = tmp_x2
     endif else begin
       x2 = tmp_x2
     endelse
     if y1 GT tmp_y2 then begin
       y2 = y1
       y1 = tmp_y2
     endif else begin
       y2 = tmp_y2
     endelse
     message, /cont, 'Chosen box: x1='+strcompress(string(x1),/rem)+', x2='+strcompress(string(x2),/rem)+', y1='+$
              strcompress(string(y1),/rem)+', y2='+strcompress(string(y2),/rem)
     index.CC_COORD=[[x1,y1],[x2,y2]]
  endif

  print,'Beginning Alignment'
  dejitter,cube_i,cube_v,cube_w,index,config,ccdebug=ccdebug,dosob=dosob
  
  ;Update index with information
  index.COMMENTS=index.COMMENTS+'Data aligned with dejitter. '

   
  ;Second save point
  res_ss=stop_and_save(cube_i,cube_v,cube_w,index,config,splt_cube=splt_cube)
  IF res_ss EQ 0 THEN return
ENDIF





;================================================================
; Calculation of wave propagation angle via FFT of velocity data cube
;================================================================

IF keyword_set(compute_waveang) THEN BEGIN

      print,'Begin FFT of data cube


      nt=index.NFRAMES
      nx=index.NAXIS1
      ny=index.NAXIS2
      nspec = nt/2
      spec=complexarr(nx,ny,nspec)
      h=do_apod(nt,cpg) ; apodisation window 

      mn_val=cube_v-rebin(mean(cube_v,dim=3),nx,ny,nt)
      h_cub=transpose(rebin(h,nt,nx,ny),[1,2,0])
      spec=(fft(mn_val*h_cub,dim=3))[*,*,0:nspec-1]  

      ; FOR iy=0,ny-1 DO $
      ; FOR ix=0,nx-1 DO BEGIN

      ;   IF wt_var_ang.mask[ix,iy] EQ 1 THEN BEGIN ;RJM 29052015
      ; 	  d=reform(cube_v[ix,iy,*]-mean(cube_v[ix,iy,*])) 
      ;   	sp=fft(d*h)                ;multiply by apodisation window
      ;     spec[ix,iy,*]=sp[0:nspec-1]
      ;   ENDIF
      ; ENDFOR

      wave_angle_calc,spec,index,config,wave_angle=wave_angle,$
                          angle_error=angle_error,coh_measure=coh_measure,save_coh=save_coh,debug=debug


      ;update mask based on coherence results from wave angle calculation
      coh_qual = (coh_measure[*,*,0]-coh_measure[*,*,1])/(coh_measure[*,*,0]+coh_measure[*,*,1])
      coh_qual[where(finite(coh_qual) EQ 0)] = 0.
      index.mask[where(coh_qual LE 0)]=0.

      ;Makes for nicer plot - but not quite correct
      wave_ang_plot=wave_angle
      wave_ang_plot[where(wave_ang_plot EQ 0)]=-90
      ang_plot=image(wave_ang_plot,xtitle='Solar X',ytitle='Solar Y',/buffer,title='Wave angle map')
      ang_plot.save,config['outpath']+'pretty_wave_angle.png'

      index=create_struct(temporary(index),'COH_MEAS',coh_qual,'coh_dist',coh_measure)

      nwlst = strcompress(string(config['numberWaveLengths']),/rem)
      IF keyword_set(save_coh) THEN save, coh_measure, file=config['outpath']+'coherence_measure_'+config['filenameExtra']+date+'_'+nwlst+$
                                          '_'+freqst+'.sav'

      save,wave_angle,angle_error,index,file=config['outpath']+'wave_angle_'+config['filenameExtra']+date+'_'+nwlst+'_'+freqst+'.sav'

      ;Third save point
      ;Only difference is updated index
      res_ss=stop_and_save(cube_i,cube_v,cube_w,index,config,splt_cube=splt_cube)
      IF res_ss EQ 0 THEN return

ENDIF

IF keyword_set(compute_speeds) THEN BEGIN
     nwlst = strcompress(string(config['numberWaveLengths']),/rem)
     ft_res=file_test(config['outpath']+'wave_angle_'+config['filenameExtra']+date+'_'+nwlst+'_'+freqst+'.sav')
     IF ft_res EQ 0 THEN BEGIN
            print,'File does not exist: wave_angle_'+config['filenameExtra']+date+'_'+nwlst+'_'+freqst+'.sav'
            print, 'Check date or run wave angle calculation first'
            return
     ENDIF
     print,'Restoring wave angles: wave_angle_'+config['filenameExtra']+date+'_'+nwlst+'_'+freqst+'.sav'
     restore,config['outpath']+'wave_angle_'+config['filenameExtra']+date+'_'+nwlst+'_'+freqst+'.sav'

     print,'Calculating propagation velocities'
     space_time_run, date, cube_v, wave_angle, index, config, config['filenameExtra'], $
                            debug=debug, wr_speeds=wr_speeds
ENDIF


IF keyword_set(plot_maps) THEN BEGIN
   IF  keyword_set(save_coh) THEN BEGIN
       plot_wave_maps, date, nwl, freqst, outpath, config['filenameExtra'], /cohm
   ENDIF ELSE BEGIN 
       plot_wave_maps, date, nwl, freqst, outpath, config['filenameExtra']
   ENDELSE 
ENDIF

print, date+' - all done!'

END
