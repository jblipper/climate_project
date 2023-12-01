import xarray as xr
import matplotlib.pyplot as plt
import sys
import scipy.optimize as so
import numpy as np
z=sys.argv[1]
start=z.find('(')
mid=z.find(',')
end=z.find(')')
point_sel=z[start:end+1]
lon_sel=z[start+1:mid]
lat_sel=z[mid+1:end]
if z[-4]=='t':
     var='Two-Metre Temperature'
     mode='t2m'
if z[-4]=='p':
     var='Total Precipitation'
     mode='tp'
path=f'/Users/justinlipper/project/data/results/{point_sel}/{z}' #define path to input netCDF file
outpath=f'/Users/justinlipper/project/data/results/{point_sel}'
dataset=xr.open_dataset(path)
data = dataset.variables['R'][:]
x=dataset.variables['longitude'][:]
y=dataset.variables['latitude'][:]
plt.contourf(x,y,data,cmap='viridis') #creates contour plot of the t2m data
cbar=plt.colorbar(label='Correlation Coefficient (R)') #plots a color coded legend to match colors to temperature values
cbar.set_label('Correlation Coefficient (R)',fontsize=25)
plt.xlabel("Longitude (Degrees)",size=30) #labels longitude axis
plt.ylabel("Latitude (Degrees)",size=30) #labels latitude axis
plt.plot(float(lon_sel),float(lat_sel),'x',color='red',label=point_sel)
plt.title(f'Correlation with {point_sel} in {var}',size=25) #adds title to plot
plt.legend(loc='upper left',fontsize=25) 
plt.savefig(f'{outpath}/plot1_{z[-4]}')
plt.show() #displays the plot

if mode=='tp':
     plt.contourf(x,y,data,cmap='coolwarm',levels=[-1,1/(np.e),1]) #creates contour plot of the t2m data
if mode=='t2m':
     plt.contourf(x,y,data,cmap='coolwarm',levels=[-1,1/(np.e)**(1/4),1]) #creates contour plot of the t2m data
cbar=plt.colorbar(label='Correlation Coefficient (R)') #plots a color coded legend to match colors to temperature values
cbar.set_label('Correlation Coefficient (R)',fontsize=25)
plt.xlabel("Longitude (Degrees)",size=25) #labels longitude axis
plt.ylabel("Latitude (Degrees)",size=25) #labels latitude axis
plt.plot(float(lon_sel),float(lat_sel),'x',color='black',label=point_sel)
if mode=='t2m':
     plt.title(f'Correlation with {point_sel} in {var} Compared to $1/e^{{1/4}}$',size=25) #adds title to plot
if mode=='tp':
     plt.title(f'Correlation with {point_sel} in {var} Compared to 1/e',size=25) #adds title to plot
plt.legend(loc='upper left',fontsize=25)
plt.savefig(f'{outpath}/plot2_{z[-4]}')
plt.show() #displays the plot

datad = dataset.variables['D'][:]
xd=dataset.variables['longitude'][:]
yd=dataset.variables['latitude'][:]

rs=[]
ds=[]
for i in range(len(x)):
     for j in range(len(y)):
          rs.append(data[i,j])
          ds.append(datad[i,j])
plt.plot(ds,rs,'.')

def model (d,CDD): # defines linear model to make fit for linear trend
     q=1
     r=np.e**(-np.array(d)*q/CDD)
     return r
best_params,cov_matrix=so.curve_fit(model,xdata=ds,ydata=rs,p0=[100]) # calculates best fit parameters for linear fit
trend=model(ds,best_params[0]) # calculates temperature anomaly predicted at each month of the 30 year time period based on the linear trend
print(best_params[0])
with open(f'{outpath}/CDD_{z[-4]}:{best_params[0]}km', 'w') as file:
     file.write(f'R:{best_params[0]}km')

if mode=='t2m':
     vv=np.ones(len(ds))*1/(np.e)**(1/4)
     plt.plot(ds,vv,linestyle=':', label='$R=1/e^{{1/4}}$')
     plt.plot(ds,trend, 'r-', label=f'$e^{{-d/CDD}}$, D_crit=CDD/4={best_params[0]/4:.2f}km') # plots linear trend on monthly temperature anomaly plot
     plt.plot(best_params[0]/4,1/(np.e)**(1/4),'X',color='black')
if mode=='tp':
     vv=np.ones(len(ds))*1/np.e
     plt.plot(ds,vv,linestyle=':', label='$R=1/e$')
     plt.plot(ds,trend, 'r-', label=f'$e^{{-d/CDD}}$, D_crit=CDD={best_params[0]:.2f}km') # plots linear trend on monthly temperature anomaly plot
     plt.plot(best_params[0],1/np.e,'X',color='black')
plt.legend(fontsize=25,loc='upper right') # adds key to plot
plt.xlabel(f'Distance From {point_sel} (km)',size=25)
plt.ylabel(f'Correlation (R) in {var}',size=25)
plt.title(f'Correlation in {var} vs. Distance from {point_sel}',size=25)
plt.savefig(f'{outpath}/plot3_{z[-4]}')
plt.show() # displays plot
resid=np.array(rs)-trend
plt.plot(ds,resid,'.')
plt.xlabel(f'Distance From {point_sel} (km)',size=25)
plt.ylabel(f'Residual of Correlation',size=25)
plt.title(f'Residuals of Correlation in {var} vs. Distance from {point_sel}',size=25)
plt.savefig(f'{outpath}/plot4_{z[-4]}')
plt.show()
