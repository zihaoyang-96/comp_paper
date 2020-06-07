;PURPOSE: Align data cubes
;
;INPUT: cube_i - intensity cube
;       cube_v - velocity cube
;       cube_w - line width cube
;       config - dictionary of variables shared by different modules
;
;OPTIONAL INPUTS: /dosob - uses sobel filtered i cube for alignment.
;                          More distinct features in gradient image for CC to 'lock-on to'
;
;OUTPUTS - the data cubes for i,v, and w are returned 'aligned'  
;
;
;CALLS: fg_rigidalign2.pro, tr_get_disp_2b.pro, shift_img.pro, add_tag.pro
;
;HISTORY: -Created by C Bethge 2012
;        - Re-wrote routine to use fgrigidalign2.pro & tided up - RJ Morton 06/2015
;        - Added option to align using Sobel Filtered image cube - RJ Morton 11/2015
;        - Edited as part of Python port - RJM 03/2019
;        - v0.3 (RJM 05/2019) Introducing versioning to keep track of code development. Given past major changes, 
;          the first numbered version is as version 0.3.




pro dejitter, cube_i, cube_v, cube_w,index, config,ccdebug=ccdebug, dosob=dosob



;define box for dejittering
IF NOT keyword_set(dosob) THEN cube_i_sm = $
       cube_i[index.coordinatesCCBox[0,0]:index.coordinatesCCBox[0,1],index.coordinatesCCBox[1,0]:index.coordinatesCCBox[1,1],*] $
ELSE BEGIN

  sob=cube_i
  FOR i=0,index.NFRAMES-1 DO sob[*,*,i]=sobel(cube_i[*,*,i])
  cube_i_sm=sob[index.coordinatesCCBox[0,0]:index.coordinatesCCBox[0,1],index.coordinatesCCBox[1,0]:index.coordinatesCCBox[1,1],*]

ENDELSE


result=0
ii=0
offset=dblarr(2,index.NFRAMES+1)


IF keyword_set(ccdebug) THEN BEGIN
    cube_i_sm_v = reform(cube_i[*,*,0])
    cube_i_sm_v[index.coordinatesCCBox[0,0],index.coordinatesCCBox[1,0]:index.coordinatesCCBox[1,1]] = max(cube_i_sm_v)
    cube_i_sm_v[index.coordinatesCCBox[0,1],index.coordinatesCCBox[1,0]:index.coordinatesCCBox[1,1]] = max(cube_i_sm_v)
    cube_i_sm_v[index.coordinatesCCBox[0,0]:index.coordinatesCCBox[0,1],index.coordinatesCCBox[1,0]] = max(cube_i_sm_v)
    cube_i_sm_v[index.coordinatesCCBox[0,0]:index.coordinatesCCBox[0,1],index.coordinatesCCBox[1,1]] = max(cube_i_sm_v)
  
    window, /free, xs=620, ys=620
    debwinnr_one = !d.window
    aia_lct, wave=193, /load
    tvscl, cube_i_sm_v > 0 < 25
    loadct, 0, /silent
    !p.multi=[0,1,2]
    window, /free
    debwinnr_two = !d.window
ENDIF

;RJM added JUL2014
bxsz=size(cube_i_sm)
dx=bxsz(1)
dy=bxsz(2)
nt=index.NFRAMES

repeat begin
 
 fg_rigidalign2,cube_i_sm,cube_i_out,dx=dx,dy=dy,dispvecs=temp,/quiet,x0=0,y0=0,/no_apply
 cube_i_sm=shift_img(temporary(cube_i_sm),temp)
 offset=offset+temp
 if (max(abs(0-temp[0,*])) lt (config['crossCorr'].threshold) and max(abs(0-temp[1,*])) lt (config['crossCorr'].threshold)) $
                              or ii eq (config['crossCorr'].maxNumberIterations)-1 then result=1
 ii=ii+1
 
 IF keyword_set(ccdebug) THEN BEGIN
    print, 'cross-correlation iterations (tr_get_disp_2b):'
    plot, temp[0,*], title='shift '+strcompress(ii)+' in x ' ; x shift
    oplot, findgen(2*nt), fltarr(2*nt)+(config['crossCorr'].threshold), linestyle=2
    oplot, findgen(2*nt), fltarr(2*nt)-(config['crossCorr'].threshold), linestyle=2
    plot, temp[1,*], title='shift '+strcompress(ii)+' in y ' ; y shift as a function of time
    oplot, findgen(2*nt), fltarr(2*nt)+(config['crossCorr'].threshold), linestyle=2
    oplot, findgen(2*nt), fltarr(2*nt)-(config['crossCorr'].threshold), linestyle=2
    wait, 0.1
    print,ii
 ENDIF

endrep until result

print,'======================================='
print,'Number of iterations performed ',ii
index.rmsOfCrossCorr=[sqrt(mean(offset[0,*]^2)),sqrt(mean(offset[1,*]^2))]
print,'Maximum CC value ', max(abs(offset))
print,'RMS CC value - x,y ',index.rmsOfCrossCorr[0],index.rmsOfCrossCorr[1]
print,'Maximum CC value on last iteration', max(abs(temp))
print,'RMS CC value on last iteration - x,y ',sqrt(mean(temp[0,*]^2)),sqrt(mean(temp[1,*]^2))

index=add_tag(index,offset,'imageShiftValues')

print,'Shifting cubes - may take a little while'
print,'======================================='
cube_i = shift_img(cube_i,offset)
cube_v = shift_img(cube_v,offset)
cube_w = shift_img(cube_w,offset)

if keyword_set(ccdebug) then begin
 !p.multi = 0
 wset, debwinnr_one
 aia_lct, wave=193, /load
 for ii=0,nt-1 do begin
  tvscl, cube_i[*,*,ii] > 0 < 25
  wait, 0.05
 endfor
 loadct, 0, /silent

 wdelete, debwinnr_one
 wdelete, debwinnr_two
endif

index.maxOcculterOffset = round(2.*max(offset))

end
