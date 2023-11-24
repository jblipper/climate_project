import netCDF4 as nc
import xarray as xr
import numpy as np
from math import sin,cos,atan2,pi
point_sel=(-64.25,-20)
file='sa3-p.nc'
path=f'/Users/justinlipper/project/data/{file}' #define path to input netCDF file
ds=xr.open_dataset(path)
#print(ds)
def get_time_series(point):
     p_sel=ds.sel(latitude=point[1],longitude=point[0])
     if file[4]=='p':
          timeseries=p_sel["tp"].values.tolist()
     if file[4]=='t':
          timeseries=p_sel["t2m"].values.tolist()
     return timeseries
sel_timeseries=get_time_series(point_sel)

def get_R(point):
     ser1=sel_timeseries
     ser2=get_time_series(point)
 #    print(list(np.array(ser1)-np.array(ser2)))
     R=np.corrcoef(ser1,ser2)[0,1]
     return R

def get_D(point):
     R=6371000 # radius of Earth(m)
     x1=point_sel[0]*pi/180
     y1=point_sel[1]*pi/180
     x2=point[0]*pi/180
     y2=point[1]*pi/180
     a=sin((y2-y1)/2)**2+cos(y1)*cos(y2)*sin((x2-x1)/2)**2
     c=2*atan2(a**0.5,(1-a)**0.5)
     d=R*c
     return d

x=get_R(point_sel)
print(x)

c=get_D(point_sel)
print(c)

def get_Rs(d):
     lats=ds.variables['latitude'][:].where(abs(point_sel[1]-ds['latitude'])<=d).dropna('latitude')
     lons=ds.variables['longitude'][:].where(abs(point_sel[0]-ds['longitude'])<=d).dropna('longitude')
     outpath='out.nc'
     outfile=nc.Dataset(outpath,'w')
     outfile.createDimension('latitude',len(lats))
     outfile.createDimension('longitude',len(lons))
     latitude=outfile.createVariable('latitude','f4',('latitude',))
     longitude=outfile.createVariable('longitude','f4',('longitude',))
     R=outfile.createVariable('R','f4',('latitude','longitude'))
     D=outfile.createVariable('D','f4',('latitude','longitude'))
     latitude[:]=lats
     longitude[:]=lons
     outfile.close()
     print(lats)
     print(lons)
     pointest=[point_sel[0],point_sel[1]]
     i=1
     out=xr.open_dataset('out.nc')
     out['R'].loc[{'latitude':pointest[1],'longitude': pointest[0]}]=get_R(pointest)
     print(pointest, get_R(pointest))
     while np.abs(pointest[0]-point_sel[0])<=d-0.25 and np.abs(pointest[1]-point_sel[1])<=d-0.25:
         for j in range(1,i+1):
              pointest[1]=pointest[1]+1/4
              out['R'].loc[{'latitude':pointest[1],'longitude': pointest[0]}]=get_R(pointest)
              print(pointest, get_R(pointest))
         for j in range(1,i+1):
              pointest[0]=pointest[0]-1/4
              out['R'].loc[{'latitude':pointest[1],'longitude': pointest[0]}]=get_R(pointest)
              print(pointest,get_R(pointest))
         for j in range(1,i+2):
              pointest[1]=pointest[1]-1/4
              cr=get_R(pointest)
              out['R'].loc[{'latitude':pointest[1],'longitude': pointest[0]}]=get_R(pointest)
              print(pointest,get_R(pointest))
         for j in range(1,i+2):
              pointest[0]=pointest[0]+1/4  
              out['R'].loc[{'latitude':pointest[1],'longitude': pointest[0]}]=get_R(pointest)
              print(pointest,get_R(pointest)) 
         i=i+2
     for j in range(i-1):
         pointest[1]=pointest[1]+1/4
         out['R'].loc[{'latitude':pointest[1],'longitude': pointest[0]}]=get_R(pointest)
         print(pointest, get_R(pointest))
     out.to_netcdf(f'out_({point_sel[0]},{point_sel[1]}){file}')
     
get_Rs(3)     
#t=file.variables["time"]
#i=file.sel(time=time[0])
#print(time[0])
