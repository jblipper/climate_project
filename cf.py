import netCDF4 as nc
import xarray as xr
import numpy as np
point_sel=(-65,-20)
path='/Users/justinlipper/project/data/sa3-p.nc' #define path to input netCDF file
ds=xr.open_dataset(path)
#print(ds)
def get_time_series(point):
     p_sel=ds.sel(latitude=point[1],longitude=point[0])
     timeseries=p_sel["tp"].values.tolist()
     return timeseries
sel_timeseries=get_time_series(point_sel)

def get_R(point):
     ser1=sel_timeseries
     ser2=get_time_series(point)
 #    print(list(np.array(ser1)-np.array(ser2)))
     R=np.corrcoef(ser1,ser2)[0,1]
     return R

x=get_R(point_sel)
print(x)

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
     latitude[:]=lats
     longitude[:]=lons
     outfile.close()
     print(lats)
     print(lons)
     pointest=[point_sel[0],point_sel[1]]
     i=1
     out=xr.open_dataset('out.nc')
     while np.abs(pointest[0]-point_sel[0])<=d-1 and np.abs(pointest[1]-point_sel[1])<=d-1:
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
     out.to_netcdf('outfile.nc')
get_Rs(3)     
#t=file.variables["time"]
#i=file.sel(time=time[0])
#print(time[0])
