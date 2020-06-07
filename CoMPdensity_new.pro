;pro CoMPdensity_new
; NOTE: before running this program, you'll need the datacube containing the theoretical density ratio curves at different height.
;     To obtain this, you can use either the widget command 'dens_plotter', or the command 'density_ratios'. For 'dens_plotter', just type ...dens_plotter, 'fe_13'...,
;     for 'density_ratios', type ...density_ratios, 'fe_13', 10747., 10801., 5, 13, den, rat, desc, radtemp=5900, rphot=1.05... Remember that you need to select the 
;     two Fe XIII lines using either way. If using the widget, you should select 'include photoexcitation' and change the distance (in solar radii) accordingly; if using 
;     the command line, you should change rphot (representing the distance away from solar center, in soalr radii).

date='20161014'
cd,'/Data/comp_2019/'+date+'/dens_new'


mask=intarr(620,620)    ;mask out pixels to invert
x=rebin(findgen(620)-310,620,620)
y=transpose(rebin(findgen(620)-311,620,620))
r=sqrt(x^2+y^2)
good=where(r gt 231. and r lt 295.)
mask(good)=1

pmask=mask; this masks the occulter and the occulter post

f1074='./1074.FitI_avg*.sav' ;the average peak intensity of 1074.7 nm
f1079='./1079.FitI_avg*.sav' ;the average peak intensity of 1079.8 nm
restore,f1079
int1079=shift(fit[*,*,0],[1,2]) ;the data should be slightly shifted to coalign with the wave observation data

restore,f1074
int1074=shift(fit[*,*,0],[1,2])

rat=int1079/int1074 ;the observed intensity ratio


restore,'./rat_'+date+'.sav' ;the save file containing theoretical relationship from CHIANTI
; in the save file, 'h' is an array of heights in solar radii (from 1.05 to 1.40); 
;   'den' is an array of density (from 5.0 to 13.0), 'ratio_array' is a 2D array containing line ratios
;   corresponding to different density values at different heights.

nx=620  &  ny=620
rinner = 231. ;inner boundary
router = 295. ;outer boundary
rsun=220. ;solar radius, pix
density=dblarr(nx,ny)
for xx=0,nx-1 do begin
  for yy=0,ny-1 do begin
    rpos = sqrt((xx-float(nx/2)-1.)^2 + (yy-float(ny/2)-1)^2)
    subh=where(h ge rpos/rsun) ;find the corresponding height (in solar radii) for each pixel
    if (rpos ge rinner) and (rpos le router) then begin
      density[xx,yy]=interpol(den,float(reform(ratio_array[subh[0],*])),rat[xx,yy],/LSQUADRATIC,/NAN) ;find the corresponding density value from the density-ratio curve from CHIANTI
    endif else begin
      rat[xx,yy]=0 & density[xx,yy]=0 & int1074[xx,yy]=0 & int1079[xx,yy]=0
    endelse
  endfor
endfor

save,filename='./densitymap_new_'+date+'.sav',rat,density,int1074,int1079


new_int1074=int1074*pmask
new_int1079=int1079*pmask
new_rat=rat*pmask
new_dens=density*pmask
save,filename='./dens_map_new_'+date+'.sav',new_dens,density,pmask,new_int1079,new_int1074,new_rat
set_plot,'x'
end
