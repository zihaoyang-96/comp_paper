pro plot_wave_maps, date, nwl, freqst, outpath, name_addon, cohm=cohm

case nwl of
 3: nwlst = '3'
 5: nwlst = '5'
endcase

if n_elements(outpath) eq 0 then begin
 outpath = '/Users/bethge/wave_tracking_output/'+date+'/'   ;<- same as in wave_tracking.pro
endif
if file_test(outpath+'plots/',/DIR) eq 0 then spawn, 'mkdir '+outpath+'plots/'
restore, outpath+'speeds_'+name_addon+date+'_'+nwlst+'_'+freqst+'.sav'
restore, outpath+'wave_angle_'+name_addon+date+'_'+nwlst+'_'+freqst+'.sav'
if keyword_set(cohm) then begin
 restore, outpath+'coherence_measure_'+name_addon+date+'_'+nwlst+'_'+freqst+'.sav'
endif

nx=(size(pro_speed))[1]
ny=(size(pro_speed))[2]
solar_x=findgen(nx)
solar_y=findgen(ny)

phase_speed     = fltarr(620,620)
phase_speed_err = fltarr(620,620)
;use the error weighted inward and outward phase speed to recalculate the phase speed 
for iy=0,ny-1 do begin
 for ix=0,nx-1 do begin 
  phase_speed[ix,iy]=(pro_speed[ix,iy]/pro_speed_err[ix,iy]^2 - ret_speed[ix,iy]/ret_speed_err[ix,iy]^2)/(1.0/pro_speed_err[ix,iy]^2 + 1.0/ret_speed_err[ix,iy]^2)
  phase_speed_err[ix,iy]=(1.0/pro_speed_err[ix,iy]^2 + 1.0/ret_speed_err[ix,iy]^2)^(-0.5)
 endfor
endfor

sub_bad=where(pro_speed le 0)
colbarpos = [0.24,0.45,0.24+0.525,0.45+0.03]

;wave_angle[where(mask eq 0)] = 0.
wave_angle = wave_angle mod 180.
bad_ang = where(wave_angle lt 0.)
if (size(bad_ang))[0] eq 1 then wave_angle[bad_ang]=wave_angle[bad_ang]+180.

set_plot, 'PS'

device,filename=outpath+'plots/wave_propagation_angle_'+name_addon+date+'_'+freqst+'.ps',xsize=13.,ysize=13.,bits_per_pixel=24, /color
loadct, 4, /silent
tvlct, r, g, b, /get
;r=reverse(r) & g=reverse(g) & b=reverse(b)
r[255]=255 & g[255]=255 & b[255]=255
r[254]=0 & g[254]=0 & b[254]=0
tvlct, r, g, b
wa = bytscl(wave_angle,0,180,top=253)
wa[where(wave_angle eq 0)] = 254
tv, wa
colorbar, position=colbarpos,chars=0.85, title='Wave propagation angle [degrees]', $
          range=[0,180], font=0, divisions=6, color=255, ncolors=253
xyouts, 0.83, 0.96, date, /norm, color=255
device,/close

device,filename=outpath+'plots/outward_phase_speed_'+name_addon+date+'_'+freqst+'.ps',xsize=13.,ysize=13.,bits_per_pixel=24, /color
pro_speed[sub_bad]=0
bpro_speed = bytscl(pro_speed*1000., min=200, max=1000, to=253)
bpro_speed[where(bpro_speed eq 0)] = 254
tv, bpro_speed
colorbar, position=colbarpos,chars=0.85, title='Outward phase speed [km/s]', $
          range=[200,1000], font=0, divisions=4, color=255, ncolors=253
xyouts, 0.83, 0.96, date, /norm, color=255
device,/close

device,filename=outpath+'plots/inward_phase_speed_'+name_addon+date+'_'+freqst+'.ps',xsize=13.,ysize=13.,bits_per_pixel=24, /color
ret_speed[sub_bad]=0
bret_speed = bytscl(-ret_speed*1000., min=200, max=1000, top=253)
tv, bret_speed
colorbar, position=colbarpos,chars=0.85, title='Inward phase speed [km/s]', $
          range=[200,1000], font=0, divisions=4, color=255, ncolors=253
xyouts, 0.83, 0.96, date, /norm, color=255
device,/close

device,filename=outpath+'plots/phase_speed_'+name_addon+date+'_'+freqst+'.ps',xsize=13.,ysize=13.,bits_per_pixel=24, /color
phase_speed[sub_bad]=0
b_phase_speed = bytscl(phase_speed*1000., min=200, max=1000, top=253)
b_phase_speed[where(b_phase_speed eq 0)] = 254
tv, b_phase_speed
colorbar, position=colbarpos,chars=0.85, title='Phase speed [km/s]', $
          range=[200,1000], font=0, divisions=4, color=255, ncolors=253
xyouts, 0.83, 0.96, date, /norm, color=255
device,/close

device,filename=outpath+'plots/outward_phase_speed_error_'+name_addon+date+'_'+freqst+'.ps',xsize=13.,ysize=13.,bits_per_pixel=24, /color
pro_speed_err[sub_bad]=0
bpro_speed_err = bytscl(pro_speed_err*1000., min=0, max=100, top=253)
bpro_speed_err[where(bpro_speed_err eq 0)] = 254
tv, bpro_speed_err
colorbar, position=colbarpos,chars=0.85, title='Outward phase speed error [km/s]', $
          range=[0,100], font=0, divisions=4, color=255, ncolors=253
xyouts, 0.83, 0.96, date, /norm, color=255
device,/close

device,filename=outpath+'plots/inward_phase_speed_error_'+name_addon+date+'_'+freqst+'.ps',xsize=13.,ysize=13.,bits_per_pixel=24, /color
ret_speed_err[sub_bad]=0
bret_speed_err = bytscl(ret_speed_err*1000., min=0, max=100, top=253)
bret_speed_err[where(bret_speed_err eq 0)] = 254
tv, bret_speed_err
colorbar, position=colbarpos,chars=0.85, title='Inward phase speed error [km/s]', $
          range=[0,100], font=0, divisions=4, color=255, ncolors=253
xyouts, 0.83, 0.96, date, /norm, color=255
device,/close

device,filename=outpath+'plots/phase_speed_error_'+name_addon+date+'_'+freqst+'.ps',xsize=13.,ysize=13.,bits_per_pixel=24, /color
phase_speed_err[sub_bad]=0
b_phase_speed_err = bytscl(phase_speed_err*1000., min=0, max=100, top=253)
b_phase_speed_err[where(b_phase_speed_err eq 0)] = 254
tv, b_phase_speed_err
colorbar, position=colbarpos,chars=0.85, title='Phase speed error [km/s]', $
          range=[0,100], font=0, divisions=4, color=255, ncolors=253
xyouts, 0.83, 0.96, date, /norm, color=255
device,/close

device,filename=outpath+'plots/outward_power_'+name_addon+date+'_'+freqst+'.ps',xsize=13.,ysize=13.,bits_per_pixel=24, /color
pro_power[sub_bad]=0
bpro_power = bytscl(alog10(pro_power), min=-6., max=-4.5, top=253)
bpro_power[where(bpro_power eq 0)] = 254
tv, bpro_power
colorbar, position=colbarpos,chars=0.85, title='Outward power [log!D10!N km!E2!Ns!E-2!N]', $
          range=[-6.,-4.5], font=0, divisions=4, color=255, ncolors=253, format='(f5.2)'
xyouts, 0.83, 0.96, date, /norm, color=255
device,/close

device,filename=outpath+'plots/inward_power_'+name_addon+date+'_'+freqst+'.ps',xsize=13.,ysize=13.,bits_per_pixel=24, /color
ret_power[sub_bad]=0
bret_power = bytscl(alog10(ret_power), min=-6., max=-4.5, top=253)
bret_power[where(bret_power eq 0)] = 254
tv, bret_power
colorbar, position=colbarpos,chars=0.85, title='Inward power [log!D10!N km!E2!Ns!E-2!N]', $
          range=[-6.,-4.5], font=0, divisions=4, color=255, ncolors=253, format='(f5.2)'
xyouts, 0.83, 0.96, date, /norm, color=255
device,/close

device,filename=outpath+'plots/normalized_power_ratio_'+name_addon+date+'_'+freqst+'.ps',xsize=13.,ysize=13.,bits_per_pixel=24, /color
pro_power[sub_bad]=0
ret_power[sub_bad]=0
pratio = (pro_power-ret_power)/(pro_power+ret_power)
pratio[sub_bad]=0
bpratio= bytscl(pratio, min=-0.5, max=0.5, top=253)
bpratio[where(pratio eq 0)] = 254
tv, bpratio
colorbar, position=colbarpos,chars=0.85, title='Normalized power ratio [(out-in)/(out+in)]', $
          range=[-0.5,0.5], font=0, divisions=4, color=255, ncolors=253, format='(f4.1)'
xyouts, 0.83, 0.96, date, /norm, color=255
device,/close

if keyword_set(cohm) then begin
 device,filename=outpath+'plots/coherence_measure_'+name_addon+date+'_'+freqst+'.ps',xsize=13.,ysize=13.,bits_per_pixel=24, /color
 coh_l = reform(coh_measure[*,*,0])
 coh_w = reform(coh_measure[*,*,1])
 coh_l[sub_bad]=0
 coh_w[sub_bad]=0
 coh_m = (coh_l-coh_w)/(coh_l+coh_w)
 coh_m[sub_bad]=0
 bcoh_m= bytscl(coh_m, min=0, max=1, top=253)
 bpratio[where(pratio eq 0)] = 254
 tv, bpratio
 colorbar, position=colbarpos,chars=0.85, title='Coh. measure [(length-width)/(length+width)]', $
           range=[0,1], font=0, divisions=2, color=255, ncolors=253, format='(f4.1)'
 xyouts, 0.83, 0.96, date, /norm, color=255
 device,/close
endif

set_plot, 'X'
;stop

end
