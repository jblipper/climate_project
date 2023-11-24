point=(1,20)
print(point[1]);
import netcdf as nc
path='/Users/justinlipper/lproject/data/sa3-p.nc' #define path to input netCDF file
file=nc.Dataset(path,'r') #open netCDF file in 'read' mode
