pro dejitter1, cube_i, cube_v, cube_w, nt=nt, ccdebug=ccdebug, box=box, ccthresh=ccthresh, $
              maxiter=maxiter, maxoffset=maxoffset

ss = size(cube_i)

; sob=cube_i
; for i=0,nt-1 do sob[*,*,i]=sobel(cube_i[*,*,i])

;define box for dejittering
;cube_i_sm = sob[box[0]:box[1],box[2]:box[3],*]
cube_i_sm = cube_i[box[0]:box[1],box[2]:box[3],*]
; cube_i_sm_v = reform(cube_i[*,*,0])
; cube_i_sm_v[box[0],box[2]:box[3]] = max(cube_i_sm_v)
; cube_i_sm_v[box[1],box[2]:box[3]] = max(cube_i_sm_v)
; cube_i_sm_v[box[0]:box[1],box[2]] = max(cube_i_sm_v)
; cube_i_sm_v[box[0]:box[1],box[3]] = max(cube_i_sm_v)

result=0
ii=0
offset=dblarr(2,ss[3]+1)

if keyword_set(ccdebug) then begin
    window, /free, xs=620, ys=620
    debwinnr_one = !d.window
    aia_lct, wave=193, /load
    tvscl, cube_i_sm_v > 0 < 25
    loadct, 0, /silent
    !p.multi=[0,1,2]
    window, /free
    debwinnr_two = !d.window
endif

bxsz=size(cube_i_sm)
dx=bxsz(1)
dy=bxsz(2)
nt=nt

repeat begin
 fg_rigidalign2,cube_i_sm,cube_i_out,dx=dx,dy=dy,dispvecs=temp,/quiet,x0=0,y0=0,/no_apply
 cube_i_sm=shift_img(temporary(cube_i_sm),temp)
 offset=offset+temp
 if (max(abs(0-temp[0,*])) lt ccthresh and max(abs(0-temp[1,*])) lt ccthresh) or ii eq maxiter-1 then result=1
 ii=ii+1
 
 IF keyword_set(ccdebug) THEN BEGIN
    print, 'cross-correlation iterations (tr_get_disp_2b):'
    plot, temp[0,*], title='shift '+strcompress(ii)+' in x ' ; x shift
    oplot, findgen(2*nt), fltarr(2*nt)+ccthresh, linestyle=2
    oplot, findgen(2*nt), fltarr(2*nt)-ccthresh, linestyle=2
    plot, temp[1,*], title='shift '+strcompress(ii)+' in y ' ; y shift as a function of time
    oplot, findgen(2*nt), fltarr(2*nt)+ccthresh, linestyle=2
    oplot, findgen(2*nt), fltarr(2*nt)-ccthresh, linestyle=2
    wait, 0.1
    print,ii
 ENDIF
endrep until result

print,'Number of iterations performed ',ii
print,'Maximum CC value ', max(abs(0-temp[0,*]))
help,offset
;stop
cube_i = shift_img(cube_i,offset)
cube_v = shift_img(cube_v,offset)
cube_w = shift_img(cube_w,offset)

if keyword_set(ccdebug) then begin
 !p.multi = 0
 wset, debwinnr_one
 aia_lct, wave=193, /load
 for ii=0,ss[3]-1 do begin
  tvscl, cube_i[*,*,ii] > 0 < 25
  wait, 0.05
 endfor
loadct, 0, /silent
endif

if keyword_set(ccdebug) then begin
 wdelete, debwinnr_one
 wdelete, debwinnr_two
endif
print,offset
maxoffset = round(2.*max(offset))
print,maxoffset
end
