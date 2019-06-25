pro aer_harmonize
;aer_haromonize.pro - mauna loa is 712

restore, 'aeronet_all.sav';,norm,stdev,ae,aod510,aod1020, lat,lon,site,elev,days,yrs,ndat
nf=n_elements(lat)
nt=n_elements(yrs)
nm=n_elements(myrs)
nx=72 ;5deg
ny=45 ;4 deg
glon=findgen(nx)*5+2.5 - 180
glat=findgen(ny)*4+2-90
aeronet=fltarr(nx,ny,nt,3)*!values.f_nan
ndat2=fltarr(nt,3)*!values.f_nan
date=yrs+days/366.

calc=1
if calc eq 1 then begin
for i=0,nt-1 do begin ;loop over time, nt-1

ind=where(finite(aod510[*,i]) eq 1)
if n_elements(ind) gt 5 then begin
ndat2[i,0]=n_elements(ind)
;interpolating in 2D using fancy algorithm,
grid_input, lon[ind],lat[ind],reform(aod510[ind,i]),xyz,newaod,/degrees,/sphere
lons=!radeg*atan(xyz[1,*],xyz[0,*])
lats=!radeg*asin(xyz[2,*])
qhull,lons,lats,tri,/delaunay, sphere=sp
altgrid2=griddata(lons,lats,newaod,/sphere,/grid,/degrees,xout=glon,yout=glat,method='NaturalNeighbor',triangles=tri)
;qhull,lon[ind],lat[ind],tri,/delaunay,sphere=sp
;altgrid2=griddata(lon[ind],lat[ind],reform(aod510[ind,i]),/sphere,/degree,/grid,xout=glon,yout=glat,method='NaturalNeighbor',triangles=tri)
destroy, tri
aeronet[*,*,i,0]=altgrid2 
;stop
endif

ind2=where(finite(aod1020[*,i]) eq 1)
if n_elements(ind2) gt 5 then begin
ndat2[i,1]=n_elements(ind2)
;interpolating in 2D using fancy algorithm,
grid_input, lon[ind2],lat[ind2],reform(aod1020[ind2,i]),xyz,newaod,/degrees,/sphere
lons=!radeg*atan(xyz[1,*],xyz[0,*])
lats=!radeg*asin(xyz[2,*])
qhull, lons,lats,tri,/delaunay, sphere=sp
altgrid2=griddata(lons,lats,newaod,/sphere,/degrees,/grid,xout=glon,yout=glat,method='NaturalNeighbor',triangles=tri)
destroy, tri
aeronet[*,*,i,1]=altgrid2
endif

ind3=where(finite(ae[*,i]) eq 1)
if n_elements(ind3) gt 5 then begin
ndat2[i,2]=n_elements(ind3)
;interpolating in 2D using fancy algorithm,
grid_input, lon[ind3],lat[ind3],reform(ae[ind3,i]),xyz,newaod,/degrees,/sphere
lons=!radeg*atan(xyz[1,*],xyz[0,*])
lats=!radeg*asin(xyz[2,*])
qhull, lons,lats,tri,/delaunay, sphere=sp
altgrid2=griddata(lons,lats,newaod,/sphere,/degrees,/grid,xout=glon,yout=glat,method='NaturalNeighbor',triangles=tri)
destroy, tri
aeronet[*,*,i,2]=altgrid2
endif

endfor



save,filename='aeronet_harm.sav', yrs,days, glon,glat,smaeronet ; aer=lon,lat,time,[aer510,1020,ae]

;for i=0,nt-1 do begin
;contour, aeronet[*,*,i,0],glon,glat,xrange=[180,180],yrange=[-90,90],xstyle=1,ystyle=1,title=strtrim(yrs[i],2)+'_'+strtrim(days[i],2)
;endfor
endif else begin  ;calculations?
restore, 'aeronet_harm.sav'
endelse

smaeronet=aeronet*!values.f_nan
for i=0,nx-1 do begin
for j=0,ny-1 do begin
for k=0,2 do begin
  smaeronet[i,j,*,k]=smooth(aeronet[i,j,*,k],10,/nan)
endfor
endfor
endfor

aerclim=fltarr(nx,ny,2)
for i=0,1 do begin
ind=where(finite(norm[*,i]) eq 1)
grid_input, lon[ind],lat[ind],reform(norm[ind,i]),xyz,newaod,/degrees,/sphere
lons=!radeg*atan(xyz[1,*],xyz[0,*])
lats=!radeg*asin(xyz[2,*])
qhull, lons,lats,tri,/delaunay, sphere=sp
aerclim[*,*,i]=griddata(lons,lats,newaod,/sphere,/degrees,/grid,xout=glon,yout=glat,method='NaturalNeighbor',triangles=tri)
endfor

day=findgen(nt)
a510=reform(smaeronet[*,*,*,0])
a1020=reform(smaeronet[*,*,*,1])
ae2=reform(smaeronet[*,*,*,2])
a510[where(finite(a510) eq 0)]=-999.
a1020[where(finite(a510) eq 0)]=-999.
ae2[where(finite(a510) eq 0)]=-999.
;CREATING NETCDF FILE 
id=ncdf_create('harmonized_aeronet.nc',/clobber)
ncdf_control,id,/fill
;MAKE DIMENSIONS
xid=ncdf_dimdef(id,'lon',nx)
yid=ncdf_dimdef(id,'lat',ny)
tid=ncdf_dimdef(id,'time',/unlimited)
wid=ncdf_dimdef(id,'wav',2)
;DEFINE VARIABLES
did=ncdf_vardef(id,'time',[tid],/float)
lonid=ncdf_vardef(id,'lon',[xid],/float)
latid=ncdf_vardef(id,'lat',[yid],/float)
wavid=ncdf_vardef(id,'wav', [wid],/float)
v1id= ncdf_vardef(id, 'AOD510nm',[xid,yid,tid],/double)
v2id=ncdf_vardef(id, 'AOD1020nm',[xid,yid,tid],/double)
v3id=ncdf_vardef(id, 'AE',[xid,yid,tid],/double)
v4id=ncdf_vardef(id, 'AERCLIM',[xid,yid,wid],/double)

;DEFINE ATTRIBUTES
ncdf_attput, id, did, 'units', 'days since 1993-01-01'
ncdf_attput, id, lonid, 'units','degrees_east'
ncdf_attput, id, lonid, 'axis', 'X'
ncdf_attput, id, latid, 'axis','Y'
ncdf_attput, id, did, 'axis','T'
ncdf_attput, id, latid, 'units','degrees_north'
ncdf_attput, id, v1id, 'units', 'unitless'
;ncdf_attput, id, v1id, 'coordinates', 'lon lat'
ncdf_attput, id, v1id, 'missing_value', -999.
ncdf_attput, id, v2id, 'missing_value', -999.
ncdf_attput, id, v3id, 'missing_value', -999.

ncdf_control,id,/endef
;PUT DATA IN FILE
ncdf_varput, id, did, day
ncdf_varput, id, lonid, glon
ncdf_varput, id, latid, glat
ncdf_varput, id, v1id, a510
ncdf_varput, id, v2id, a1020
ncdf_varput, id, v3id, ae2
ncdf_varput, id, v4id, aerclim
ncdf_close, id

;make global mean
zonaer=mean(smaeronet,dim=1,/nan)
wgts=cos(glat*!pi/180.)
gm=fltarr(nt)
for i=0,nt-1 do begin
gm[i]=total(zonaer[*,i,0]*wgts)/total(wgts)
endfor

gm1=gm
gm2=gm
;deseasonalize
for i=0,365 do begin
ind=where(days eq i)
gm1[ind]=gm1[ind]/mean(gm[ind],/nan)
gm2[ind]=gm2[ind]-mean(gm2[ind],/nan)
endfor

;calculate trend
;fin=where(finite(gm) eq 1)
;x=findgen(nt)
;tr1=linfit(x[fin],gm1[fin])
;tr2=linfit(x[fin],gm2[fin])
;print, tr1, tr2

plot, yrs+days/365.,gm,thick=2, linestyle=1
oplot, yrs+days/365.,gm1,linestyle=2,thick=2
oplot, yrs+days/365.,gm2, thick=2

stop
END
