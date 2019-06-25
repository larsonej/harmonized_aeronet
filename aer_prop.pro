;aer_prop.pro

restore, 'aod_temp.sav'
nt=6562

;dat=read_ascii('AOD20/DAILY/19930101_20190223_Mauna_Loa.lev20',template=aod_temp, header=header)
;mldate=dat.(0)
;mlyr=intarr(nt)
;for i=0,nt-1 do begin
;  dummy=strsplit(mldate[i],':',/extract)
;  mlyr[i]=dummy[2]
;endfor
;mldoy=dat.(2)

volcname=['Calbuco','Nabro','Sarychev','Alu-Dalafilla','Kasatochi']
; calbuco traveled north and east
; kasatochi traveled east and n,s
; nabro traveled east and north over china
; sarychev traveled
; alu-D traveled E
volcyr=[2015, 2011, 2009, 2008,2008]
volcdoy=[112, 164, 166, 308, 220]
volclat=[-41, 13, 48, 14, 52]
volclon=[-73, 42, 153, 40, -175] ;E

list=file_search('AOD20/DAILY/*.lev20',count=ndat)
site=strarr(ndat)
loc=fltarr(2,ndat)
date=fltarr(2,nt,ndat)*!values.f_nan
aod=fltarr(22, nt, ndat)*!values.f_nan
AE=fltarr(6, nt, ndat)*!values.f_nan
Np=fltarr(22, nt, ndat)*!values.f_nan
wat=fltarr(nt,ndat)*!values.f_nan

for i=0,ndat-1 do begin
  print, list[i]
  file=list[i] ;'AOD20/DAILY/19930101_20190223_AAOT.lev20'
  dat=read_ascii(file, template=aod_temp, header=header)
  if i eq 0 then head=strsplit(header[6], ',',/extract)
  nj=n_elements(dat.field03)
  if nj gt 30 then begin
    site[i]=header[1]
    loc[0,i]=dat.field79[0]
    loc[1,i]=dat.field80[0]
    dummy=strsplit(dat.field01,':',/extract)
    for j=0,nj-1 do begin
      date[0,j,i]=float(dummy[j,2])
    endfor
    date[1,0:nj-1,i]=dat.field03
    for k=0,21 do begin
      aod[k,0:nj-1,i]=dat.(3+k)
      np[k,0:nj-1,i]=dat.(39+k)
    endfor
    for k=0,5 do begin
      ae[k,0:nj-1,i]=dat.(34+k)
    endfor
    wat[0:nj-1,i]=dat.(25)
  endif

endfor

save, filename='aeronet_daily.sav',/all
stop
END
