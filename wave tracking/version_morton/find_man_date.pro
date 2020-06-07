
;PURPOSE: File to store manual date information
;         Also contains default options
;
;
;Guide:   numberFilesToProcess - value given should be good frames + gaps
;
;
;          Modularised by R J Morton 2015
;          Created to keep together useful quantities from wave tracking RJM 04/2016
;          Edited as part of Python port - RJM 03/2019
;        - v0.3 (RJM 05/2019) Introducing versioning to keep track of code development. Given past major changes, 
;          the first numbered version is as version 0.3.


PRO find_man_date,date,config 

; *** Default values for a date that has not ***
; *** been checked before and is not defined ***
; *** below in the case statement            ***

; coordinates for cross-correlation box
config.coordinatesCCBox=[[45,230],[75,390]] ;[x1,y1],[x2,y2]


; Choose start file - used to skip bad files
config.startFile = 0

; Number of files to process. 0 means: let the 
; code choose the highest number possible. Useful 
; for manually cutting out bad data at the end.
config.numberFilesToProcess = 0
config.maxOcculterOffset = 0



CASE date OF
 ;example for a manual setting for a specific date

 '20111230': begin
                config.coordinatesCCBox=[[540,319],[575,371]]
                config.startFile = 0
                config.numberFilesToProcess = 179
                config.lowerMaskRadius = 228.
                config.upperMaskRadius = 290.
               end

 '20120327': begin
               config.coordinatesCCBox=[[55,220],[85,260]]
      	       config.startFile = 0
      	       config.numberFilesToProcess = 164
      	       config.lowerMaskRadius = 228.
      	       config.upperMaskRadius = 280.
              end

 '20120411': begin
              config.coordinatesCCBox=[[540,350],[580,390]]
      	      config.startFile = 0        ;361
      	      config.numberFilesToProcess = 360  ;532
      	      config.lowerMaskRadius = 228.
      	      config.upperMaskRadius = 290.
             end

 '20120410': begin
               config.coordinatesCCBox=[[535,390],[565,430]]
      	       config.startFile = 0
      	       config.numberFilesToProcess = 400
      	       config.lowerMaskRadius = 228.
      	       config.upperMaskRadius = 290.
      	      end

 '20120706': begin
               config.coordinatesCCBox=[[539,282],[563,310]]
               config.startFile = 42
               config.numberFilesToProcess = 132
               config.lowerMaskRadius = 228.
               config.upperMaskRadius = 290.
             end

 '20130213': BEGIN
              config.coordinatesCCBox=[[545,375],[580,410]]
              config.startFile = 0
              config.numberFilesToProcess = 178
              config.lowerMaskRadius = 228.
              config.upperMaskRadius = 288.
             END

 '20130502': begin
                config.coordinatesCCBox=[[540,220],[570,250]]
      	        config.startFile = 0
      	        config.numberFilesToProcess = 176
      	        config.lowerMaskRadius = 228.
      	        config.upperMaskRadius = 290.
            end

 '20130515': begin
                config.coordinatesCCBox=[[545,331],[573,372]]
                config.startFile = 0
                config.numberFilesToProcess = 173
                config.lowerMaskRadius = 228.
                config.upperMaskRadius = 290.
             end

 '20130708': begin
                config.coordinatesCCBox=[[57,241],[80,275]]
                config.startFile = 0
                config.numberFilesToProcess = 179
                config.lowerMaskRadius = 228.
                config.upperMaskRadius = 290.
             end                      

 '20130914': begin
                config.coordinatesCCBox=[[545,310],[580,345]]
                config.startFile = 0
      	        config.numberFilesToProcess = 152
      	        config.lowerMaskRadius = 228.
      	        config.upperMaskRadius = 290.
             end

 '20131223': begin
                config.coordinatesCCBox=[[61,178],[92,223]]
                config.startFile = 0
                config.numberFilesToProcess = 48
                config.lowerMaskRadius = 228.
                config.upperMaskRadius = 290.
             end

 '20140511': begin
                config.coordinatesCCBox=[[61,199],[89,231]]
                config.startFile = 0
                config.numberFilesToProcess = 178
                config.lowerMaskRadius = 228.
                config.upperMaskRadius = 290.
              end

 '20140510': begin
                config.coordinatesCCBox=[[538,345],[570,377]]
                config.startFile = 0
                config.numberFilesToProcess = 180
                config.lowerMaskRadius = 228.
                config.upperMaskRadius = 290.
             end

  '20150121': begin
                config.coordinatesCCBox=[[46,220],[77,259]]
                config.startFile = 0
                config.numberFilesToProcess = 169
                config.lowerMaskRadius = 228.
                config.upperMaskRadius = 290.
             end
  '20150122': begin
                config.coordinatesCCBox=[[46,220],[77,259]]
                config.startFile = 0
                config.numberFilesToProcess = 169
                config.lowerMaskRadius = 228.
                config.upperMaskRadius = 290.
             end                       

  '20150202': begin
                config.coordinatesCCBox=[[77,140],[112,176]]
                config.startFile = 0
                config.numberFilesToProcess = 169
                config.lowerMaskRadius = 228.
                config.upperMaskRadius = 290.
             end

  else:      begin
               message, /cont, 'Manual settings for this date were not found, using default values.'
             end
endcase



END