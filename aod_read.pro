;aod_read.pro
;EJL 1-29-17 
; this program reads aeronet data and puts it into a structure.

;function to fit trend and seasonal cycle; y = A+B*x+C*sin(D*x+E)
function myfit, x,p
  y=m*x+b+ar1*n[1:-1]
return, func
end


;Start of main program
pro aod_read

restore, 'aod_template.sav'
;spawn, 'ls AOT_Level2_All_Points/AOT/LEV20/ALL_POINTS/ > aot_files.txt'
;readcol, 'aot_files.txt', aot_files, format='(a)'
aot_files=file_search('DAILY/','*.lev20')
nf=n_elements(aot_files)

;!p.multi=[0,1,2]

;set_plot,'ps'
;plot ,findgen(30)+1990., findgen(70)/350.,/nodata, yrange=[-.01,.3], /ys

sigchop=0
noplot=0
;site_index=[31,75,103,119,150,156,157,345,380,393,410,431,473,506,544,548,561,576,598,614,616,713,752,761,779,781,807,822,846,883,889,937,942,944]
;site_index=[924];[11,186,354,545,625,692,708,757,833]
;restore,'aeronet_list.sav'
;site_index=list2

for si =0,n_elements(site_index)-1 do begin ;0,nf-1 do begin
  f=site_index[si]
  ;read data and plot
  path='AOT_Level2_All_Points/AOT/LEV20/ALL_POINTS/'
  data=read_ascii(path+aot_files[f], template=aod_temp, header=lev_head)
  header=strsplit(lev_head[4],',',/extract)

  loc=strsplit(lev_head[2],',=',/extract)
  site=loc[1]
  lon=loc[3]
  lat=loc[5]
  elev=loc[7]
  print, site

;  if (lon ge -160. and lon le -50) then begin
    dmy=strsplit(data.(0),':',/extract)
    dmy=dmy.toarray(type='INT')
    hms=strsplit(data.(1),':',/extract) ;GMT
    hms=hms.toarray(type='Float')
    aod500=float(data.(12))
    dataind=where(aod500 le 0.005, count, complement=gooddat)
    aod500[dataind]=!values.f_nan
    if (sigchop eq 1) then begin
      cutoff=mean(aod500,/nan)+3*stddev(aod500,/nan)
      toohigh=where(aod500 gt cutoff)
      aod500[toohigh]=!values.f_nan
    endif
    dates=double(dmy[*,2]+(dmy[*,1]-1.)/12.+dmy[*,0]/365.+(hms[*,0]-10.)/(365.*11.))
    diy=double((dmy[*,1]-1.)/12.+dmy[*,0]/365.+(hms[*,0]-10.)/(365.*11.))
    julday=float(data.(2))

   
    ; Get yearly data ====================================================
    yrs=uniq(dmy[*,2]) ; get month end indices
    nyr=n_elements(yrs)       ;number of months
    yrs=[-1,yrs]              ; add -1 for loop below
    nyrs=yrs[1:-1]-yrs[0:-2]
    aod500_yr=fltarr(nyr) 
    yr_dates=fltarr(nyr)
    for i=0,nyr-1 do begin
       yr_dates[i]=dmy[yrs[i+1],2]
       aod500_yr[i]=mean(aod500[yrs[i]+1:yrs[i+1]],/nan)
    endfor

    
    ;Get monthly climatology ======================================
    aod500_clim=fltarr(12)
    for i=0,11 do begin
      ind=where(dmy[*,1]-1 eq i)
      aod500_clim[i]=mean(aod500[ind],/nan)
    endfor
    
    ; Get montly data and deseasonalized data
    mons=uniq(dmy[*,1])         ; get month end indices
    if (mons[-1] ne n_elements(dmy[*,1])-1) then mons=[mons,n_elements(dmy[*,1])-1] ; if first and last month same, it doesn't show up, so add last number
    month=dmy[mons,1]
    month_yr=dmy[mons,2]
    nmon=n_elements(mons)       ;number of months
    ndat=mons[1:-1]-mons[0:-2]
    ndat=[mons[0]-0,ndat]
    mons=[-1,mons] ; add -1 for loop below
    aod500_mon=fltarr(nmon) 
    aod500_ds=fltarr(nmon)
    mon_dates=fltarr(nmon)
    aod500_clim2=fltarr(nmon)
    for i=0,nmon-1 do begin
       mon_dates[i]=dmy[mons[i+1],2]+(dmy[mons[i+1],1]-.5)/12.
       aod500_mon[i]=mean(aod500[mons[i]+1:mons[i+1]],/nan)
       aod500_ds[i]=mean(aod500[mons[i]+1:mons[i+1]],/nan) - aod500_clim[dmy[mons[i+1],1]-1]
       aod500_clim2[i]=aod500_clim[dmy[mons[i+1],1]-1]
    endfor


    ; Get daily mean data ======================================
    aod500_dayclim=fltarr(366)
    for i=0,356 do begin
      ind=where(fix(julday)-1 eq i)
      aod500_dayclim[i]=mean(aod500[ind],/nan)
    endfor

    ; Get monthly data and deseasonalized data
    jd=fix(julday)
    days=uniq(jd)         ; get month end indices
    if (days[-1] ne n_elements(jd)-1) then days=[days,n_elements(jd)-1] ; if first and last month same, it doesn't show up, so add last number
    nday=n_elements(days)       ;number of days
    ndat2=days[1:-1]-days[0:-2] ; number of reads per day
    ndat2=[days[0]-0, ndat2]
    diy=jd[days]
    day_yr=dmy[days,2]
    days=[-1,days] ; add -1 for loop below
    aod500_day=fltarr(nday) 
    aod500_dayds=fltarr(nday)
    day_dates=fltarr(nday)
    aod500_dayclim2=fltarr(nday)
    for i=0,nday-1 do begin
       day_dates[i]=dmy[days[i+1],2]+(jd[days[i+1]])/356. ;finding the decimal year of the data
       aod500_day[i]=mean(aod500[days[i]+1:days[i+1]],/nan) ; getting the daily average data
       aod500_dayds[i]=mean(aod500[days[i]+1:days[i+1]],/nan) - aod500_dayclim[jd[days[i+1]]-1] ;deseasonalizing the daily average with the appropriate day data
       aod500_dayclim2[i]=aod500_dayclim[jd[days[i+1]]-1] ; the climatology
    endfor


    ;Calculate trends from the data ==============================
    fin=where(finite(aod500) eq 1)
    trend=linfit(dates[fin]-dates[0], aod500[fin])
;    print, trend, 'all data trend'
    monfin=where(finite(aod500_mon) eq 1 )
    montrend=linfit(mon_dates[monfin]-mon_dates[0], aod500_mon[monfin], sigma=monsig)
;    print, montrend, monsig,'mon data trend'
    mondstrend=linfit(mon_dates[monfin]-mon_dates[0], aod500_ds[monfin], sigma=dsmonsig, yfit=monyfit)
;    print, mondstrend,dsmonsig, 'deseasonalized mon trend'
    dayfin=where(finite(aod500_day) eq 1)
    daytrend=linfit(day_dates[dayfin]-day_dates[0], aod500_day[dayfin], sigma=daysig)
;    print, daytrend,daysig, 'day data trend'
    daydstrend=linfit(day_dates[dayfin]-day_dates[0], aod500_dayds[dayfin], sigma=dsdaysig, yfit=dayyfit)
;    print, daydstrend, dsdaysig,'deseasonalized day trend'
    yrfin=where(finite(aod500_yr) eq 1)
    yrtrend=linfit(yr_dates[yrfin]-yr_dates[0], aod500_yr[yrfin]-mean(aod500_yr,/nan), sigma=yrsig, yfit=yryfit)
;    print, yrtrend, yrsig,'yearly trend'
;    print, ''


 if (n_elements(monfin) ge 36) then begin
    yr_acor=a_correlate(aod500_yr[yrfin],indgen(4)+1)
    mon_acor=a_correlate(aod500_mon[monfin],indgen(12)+1)
    mon_dsacor=a_correlate(aod500_ds[monfin],indgen(12)+1)
    day_acor=a_correlate(aod500_day[dayfin],indgen(365)+1)
    day_dsacor=a_correlate(aod500_dayds[dayfin],indgen(365)+1)


    mon_ds_unc=dsmonsig[1]/(n_elements(monfin)/12.)^(3./2.)*sqrt((1+mon_dsacor[0])/(1.-mon_dsacor[0]))
    day_ds_unc=dsdaysig[1]/(n_elements(dayfin)/365.)^(3./2.)*sqrt((1+day_dsacor[0])/(1.-day_dsacor[0]))     
    yr_unc=yrsig[1]/(n_elements(yrfin)/365.)^(3./2.)*sqrt((1+yr_acor[0])/(1.-yr_acor[0]))
 endif
;    print, monsig, mon_ds_unc, 'month trend uncertainty before and after autocorrelation correction'
;    print, daysig, day_ds_unc, 'day trend uncertainty before and after autocorrelation correction'
;    print, yrsig, 'year trend uncetainty'

if (noplot eq 0) then begin  
  !p.multi=[0,1,4]
  loadct, 39, ncolors=256
  set_plot,'ps'
  device,filename='plots/'+site+'_trend_detail.png',/encapsulated,xsize=6*2.54, ysize=8*2.54,/color
  plot, dates, aod500, psym=3, thick=3, charthick=3, xrange=[floor(min(dates)), floor(max(dates))+1],/xs, ytitle='AOD (500nm)'
  oplot, day_dates, aod500_day, thick=1, psym=1, symsize=.2, color=60
  oplot, mon_dates, aod500_mon, thick=3, color=120
  oplot, yr_dates, aod500_yr, thick=3, psym=7,symsize=.7, color=245
  legend2, ['Raw data','Daily Averages', 'Monthly averages','Yearly Averages'], psym=[3, 1, 0, 7], color=[0,60,120,245 ], box=0, thick=3, charthick=3
  
  plot, day_dates, aod500_dayds, psym=1, thick=1, charthick=3, xrange=[floor(min(dates)), floor(max(dates))+1],/xs, ytitle='deseasonalized daily AOD (500nm) with trendline'
  oplot, day_dates[dayfin], dayyfit, thick=3, color=60
  xyouts, .15, .685, 'Trend = '+string(daydstrend[1])+' +/-'+string(day_ds_unc), charthick=3,/normal
  
  plot, mon_dates, aod500_ds, thick=2, psym=1, charthick=3, xrange=[floor(min(dates)), floor(max(dates))+1],/xs, ytitle='deseasonalized monthly AOD (500nm) with trendline'
  oplot, mon_dates[monfin], monyfit, thick=3 , color=120
  xyouts, .15, .45, 'Trend = '+string(mondstrend[1])+' +/-'+string(mon_ds_unc), charthick=3,/normal
  
  plot, yr_dates, aod500_yr-mean(aod500_yr,/nan), thick=2, psym=1, charthick=3, xrange=[floor(min(dates)), floor(max(dates))+1],/xs, ytitle='yearly AOD (500nm) with trendline'
  oplot, yr_dates[yrfin], yryfit, thick=3, color=245
  xyouts, .15, .185, 'Trend = '+string(yrtrend[1])+' +/-'+string(yr_unc), charthick=3,/normal
  device,/close
  

  !p.multi=[0,1,2]
  device,filename='plots/'+site+'_autocorrelation.png',/encapsulated,xsize=6*2.54, ysize=6*2.54,/color 
  plot, mon_acor, thick=3,charthick=3, xtitle='Month', ytitle='Autocorrelation'
  oplot, mon_dsacor, thick=3, linestyle=2
  legend2, ['Raw data', 'Deseasonalized data'], linestyle=[0,2], thick=3, charthick=3, box=0
  
  plot, day_acor, thick=3,charthick=3, xtitle='Day', ytitle='Autocorrelation'
  oplot, day_dsacor, thick=3, linestyle=2  
  device,/close

  mname=['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']
  !p.multi=0
  loadct, 39, ncolors=14
  device,filename='plots/'+site+'_allmonths_trend.eps',/encapsulated,xsize=6*2.54, ysize=6*2.54,/color
  ind=where(month eq 4)
  plot, month_yr[ind],aod500_mon[ind], thick=3,charthick=3, xtitle='Year', ytitle='Monthly AOD (500nm)',/nodata, xrange=[min(month_yr), max(month_yr)], yrange=[0,max(aod500_mon)]
  for i=1,12 do begin
    ind=where(month eq i)
    if (ind[0] ne -1) then oplot, month_yr[ind],aod500_mon[ind], thick=3, color=i
  endfor
  legend2, [mname],linestyle=fltarr(12), color=[indgen(12)+1], thick=3, charthick=3, box=0
  device,/close

endif ;noplot

  ;Write the daily and monthly mean data to a file. 
  openw, lun, 'plots/newdat/'+site+'_daily_AOD500nm.txt',/get_lun
  printf, lun, 'Day   ','Year    ','AOD 500nm ','N_pts   ' 
  for i=0,n_elements(aod500_day)-1 do begin
    printf,lun, diy[i], day_yr[i], aod500_day[i],ndat2[i], format='(I8, I8, f10.3,I8)'
  endfor
  close,/all

  openw, lun, 'plots/newdat/'+site+'_monthly_AOD500nm.txt',/get_lun
  printf, lun, 'Month   ','Year    ','AOD 500nm ','N_pts   ' 
  for i=0,n_elements(aod500_mon)-1 do begin
    printf,lun, month[i], month_yr[i], aod500_mon[i],ndat[i], format='(I8, I8, f10.3,I8)'
  endfor
  close,/all


endfor ;site loop

set_plot,'x'
stop
END
