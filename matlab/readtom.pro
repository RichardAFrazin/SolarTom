pro readtom,dir,file,nrad,nlat,rmin,rmax,N_e,rad,lat,lon
nlon=2*nlat
N_e=fltarr(nrad,nlat,nlon)
openr,1,dir+file
readu,1,N_e
close,1
drad=(rmax-rmin)/nrad
dang=180./nlat
rad = rmin + drad/2. + drad * findgen(nrad)
lat = -90. + dang/2. + dang * findgen(nlat)
lon =   0. + dang/2. + dang * findgen(nlon)
return
end
