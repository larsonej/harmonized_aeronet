pro aeronet
;aeronet.pro - mauna loa is 712
;restore, 'aod_template.sav'
restore, 'aod_temp.sav'
aot_files=file_search('AOD20/DAILY/','*.lev20')
nf=n_elements(aot_files)
nt=27l*365+6 ; 1993-2018 * 365
nm=27l*12
year=findgen(27)+1993
yrs=[]
days=[]
mons=[]
myrs=[]
for i=0,26 do begin
  if year[i] eq 1996 or year[i] eq 2000 or year[i] eq 2004 or year[i] eq 2008 or year[i] eq 2012 or year[i] eq 2016 then begin
    days=[days, findgen(366)+1]
    yrs=[yrs,fltarr(366)+year[i]]
  endif else begin
    days=[days,findgen(365)+1]
    yrs=[yrs,fltarr(365)+year[i]]
  endelse
myrs=[myrs,fltarr(12)+year[i]]
mons=[mons,findgen(12)+1]
endfor

lat=fltarr(nf)
lon=fltarr(nf)
site=strarr(nf)
elev=fltarr(nf)
aod510=fltarr(nf,nt)*!values.f_nan ; indice 16
aod1020=fltarr(nf,nt)*!values.f_nan ; indice 4
ae=fltarr(nf,nt)*!values.f_nan ;indice 35
maod510=fltarr(nf,nm)*!values.f_nan
maod1020=fltarr(nf,nm)*!values.f_nan
mae=fltarr(nf,nm)*!values.f_nan
norm=fltarr(nf,2)
stdev=fltarr(nf,2)
mnorm=fltarr(nf,2)
mstdev=fltarr(nf,2)
ndat=fltarr(nf)

wav=[1640, 1020, 870, 865, 779, 675, 667, 620, 560, 555, 551, 532, 531, 510, 500, 490, 443, 440, 412, 400, 380, 340];wavlength nm

for si =0,nf-1 do begin ;loop over sites - read data
  data=read_ascii(aot_files[si], template=aod_temp, header=lev_head)
  header=strsplit(lev_head[6],',',/extract)
  site[si]=data.(77)[0]
  lat[si]=data.(78)[0]
  lon[si]=data.(79)[0]
  elev[si]=data.(80)[0]
  nd=n_elements(data.(4))
  ndat[si]=nd
  print, site[si], nd
  if nd gt 99 then begin

;find doy, yr
  dmy=strsplit(data.(0),':',/extract) ; format dd:mm:yyyy
  dmy=dmy.toarray(type='INT')
  yr=reform(dmy[*,2])
  mon=reform(dmy[*,1])
  doy=data.(2)
  aod=fltarr(22,nd)
  a=fltarr(2,nd)

for w=0,21 do begin
  aod[w,*]=data.(3+w)
endfor
  aod[where(aod eq -999)]=!values.f_nan
;calculate AE, AOD510, AOD1020, and place into correct position in time
  for t=0,nd-1 do begin
    fin=where(finite(aod[*,t]) eq 1)
    if n_elements(fin) gt 1 then begin
      fit=linfit(alog(wav[fin]), alog(aod[fin,t]))
      ind=where(yr[t] eq yrs and doy[t] eq days)
      if ind[0] eq -1 then stop
;    print, ind
      aod510[si,ind]=exp(interpol(alog(aod[fin,t]),alog(wav[fin]),alog(510)))
      aod1020[si,ind]=exp(interpol(alog(aod[fin,t]),alog(wav[fin]),alog(1020)))
      ae[si,ind]=-1.*fit[1]
      a[0,t]=aod510[si,ind]
      a[1,t]=aod1020[si,ind]
    endif
  endfor

  for y=0,26 do begin 
  for m=0,11 do begin
    dind=where(yr eq 1993+y and mon eq m+1)
;print, dind
    if n_elements(dind) gt 2 then begin
       maod510[si,12*y+m]=mean(a[0,dind],/nan)
       maod1020[si,12*y+m]=mean(a[1,dind],/nan)
    endif
  endfor
  endfor

;normalize data
  aod510[si,*]=(aod510[si,*]-mean(aod510[si,*],/nan))/stddev(aod510[si,*],/nan)
  aod1020[si,*]=(aod1020[si,*]-mean(aod1020[si,*],/nan))/stddev(aod1020[si,*],/nan)
  maod510[si,*]=(maod510[si,*]-mean(maod510[si,*],/nan))/stddev(maod510[si,*],/nan)
  maod1020[si,*]=(maod1020[si,*]-mean(maod1020[si,*],/nan))/stddev(maod1020[si,*],/nan)
; save averages for climatology
  norm[si,0]=mean(aod510[si,*],/nan)
  norm[si,1]=mean(aod1020[si,*],/nan)
  stdev[si,0]=stddev(aod510[si,*],/nan)
  stdev[si,1]=stddev(aod1020[si,*],/nan)

  mnorm[si,0]=mean(maod510[si,*],/nan)
  mnorm[si,1]=mean(maod1020[si,*],/nan)
  mstdev[si,0]=stddev(maod510[si,*],/nan)
  mstdev[si,1]=stddev(maod1020[si,*],/nan)
;stop
endif ; ndat gt 29
endfor ;site
mae=-1.*alog(maod510/maod1020)/alog(.5)
save, filename='aeronet_all.sav',norm,stdev,ae,aod510,aod1020, lat,lon,site,elev,days,yrs,ndat, mons, myrs, maod1020,maod510,mae,mnorm,mstdev

stop
END
