import xarray as xr
import matplotlib.pyplot as plt
import sys
z=sys.argv[1]
path=f'/Users/justinlipper/project/data/{z}' #define path to input netCDF file
dataset=xr.open_dataset(path)
data = dataset.variables['R'][:]
x=dataset.variables['longitude'][:]
y=dataset.variables['latitude'][:]
plt.contourf(x,y,data,cmap='viridis') #creates contour plot of the t2m data
plt.colorbar(label='Two Metre Temperature (K)') #plots a color coded legend to match colors to temperature values
plt.xlabel("Longitude (Degrees)") #labels longitude axis
plt.ylabel("Latitude (Degrees)") #labels latitude axis
plt.title("Correlation Coefficient (R)") #adds title to plot
plt.legend(loc='upper left') 
plt.show() #displays the plot
