

```python
from netCDF4 import Dataset
import numpy as np
import time as python_time
from numpy import nanmedian
from scipy.ndimage.filters import gaussian_filter
```


```python
nc = Dataset('sg579_JanApr_L2_2m.nc','r')
nc.variables
```




    OrderedDict([(u'divenum', <type 'netCDF4._netCDF4.Variable'>
                  int16 divenum(cast)
                  unlimited dimensions: 
                  current shape = (1038,)
                  filling on, default _FillValue of -32767 used),
                 (u'downcast', <type 'netCDF4._netCDF4.Variable'>
                  int8 downcast(cast)
                  unlimited dimensions: 
                  current shape = (1038,)
                  filling on, default _FillValue of -127 ignored),
                 (u'pressure', <type 'netCDF4._netCDF4.Variable'>
                  float64 pressure(pressure)
                      units: db
                      scale_factor: 1.0
                      add_offset: 0.0
                  unlimited dimensions: 
                  current shape = (500,)
                  filling on, default _FillValue of 9.96920996839e+36 used),
                 (u'DAC_U', <type 'netCDF4._netCDF4.Variable'>
                  float64 DAC_U(cast)
                      _FillValue: -999.0
                      units: m/s
                      scale_factor: 1.0
                      add_offset: 0.0
                  unlimited dimensions: 
                  current shape = (1038,)
                  filling on),
                 (u'DAC_V', <type 'netCDF4._netCDF4.Variable'>
                  float64 DAC_V(cast)
                      _FillValue: -999.0
                      units: m/s
                      scale_factor: 1.0
                      add_offset: 0.0
                  unlimited dimensions: 
                  current shape = (1038,)
                  filling on),
                 (u'latitude', <type 'netCDF4._netCDF4.Variable'>
                  float64 latitude(cast, pressure)
                      _FillValue: -999.0
                      units: degrees north
                      scale_factor: 1.0
                      add_offset: 0.0
                  unlimited dimensions: 
                  current shape = (1038, 500)
                  filling on),
                 (u'longitude', <type 'netCDF4._netCDF4.Variable'>
                  float64 longitude(cast, pressure)
                      _FillValue: -999.0
                      units: degrees east
                      scale_factor: 1.0
                      add_offset: 0.0
                  unlimited dimensions: 
                  current shape = (1038, 500)
                  filling on),
                 (u'latitude_flightmodel', <type 'netCDF4._netCDF4.Variable'>
                  float64 latitude_flightmodel(cast, pressure)
                      _FillValue: -999.0
                      units: degrees north
                      scale_factor: 1.0
                      add_offset: 0.0
                  unlimited dimensions: 
                  current shape = (1038, 500)
                  filling on),
                 (u'longitude_flightmodel', <type 'netCDF4._netCDF4.Variable'>
                  float64 longitude_flightmodel(cast, pressure)
                      _FillValue: -999.0
                      units: degrees east
                      scale_factor: 1.0
                      add_offset: 0.0
                  unlimited dimensions: 
                  current shape = (1038, 500)
                  filling on),
                 (u'time', <type 'netCDF4._netCDF4.Variable'>
                  float64 time(cast, pressure)
                      _FillValue: -999.0
                      units: days since January 1, 2013
                      scale_factor: 1
                      add_offset: 0.0
                  unlimited dimensions: 
                  current shape = (1038, 500)
                  filling on),
                 (u'temperature', <type 'netCDF4._netCDF4.Variable'>
                  float64 temperature(cast, pressure)
                      _FillValue: -999.0
                      units: deg C
                      scale_factor: 1
                      add_offset: 0.0
                  unlimited dimensions: 
                  current shape = (1038, 500)
                  filling on),
                 (u'salinity', <type 'netCDF4._netCDF4.Variable'>
                  float64 salinity(cast, pressure)
                      _FillValue: -999.0
                      units: Practical Salinity Units
                      scale_factor: 1
                      add_offset: 0.0
                  unlimited dimensions: 
                  current shape = (1038, 500)
                  filling on),
                 (u'density', <type 'netCDF4._netCDF4.Variable'>
                  float64 density(cast, pressure)
                      _FillValue: -999.0
                      units: kg/m^3
                      scale_factor: 1
                      add_offset: 1000.0
                  unlimited dimensions: 
                  current shape = (1038, 500)
                  filling on),
                 (u'pot_temperature', <type 'netCDF4._netCDF4.Variable'>
                  float64 pot_temperature(cast, pressure)
                      _FillValue: -999.0
                      units: deg C
                      scale_factor: 1
                      add_offset: 0.0
                  unlimited dimensions: 
                  current shape = (1038, 500)
                  filling on),
                 (u'pot_density', <type 'netCDF4._netCDF4.Variable'>
                  float64 pot_density(cast, pressure)
                      _FillValue: -999.0
                      units: kg/m^3
                      scale_factor: 1
                      add_offset: 1000.0
                  unlimited dimensions: 
                  current shape = (1038, 500)
                  filling on),
                 (u'n_time', <type 'netCDF4._netCDF4.Variable'>
                  int8 n_time(cast, pressure)
                      _FillValue: 25
                  unlimited dimensions: 
                  current shape = (1038, 500)
                  filling on),
                 (u'n_temperature', <type 'netCDF4._netCDF4.Variable'>
                  int8 n_temperature(cast, pressure)
                      _FillValue: 25
                  unlimited dimensions: 
                  current shape = (1038, 500)
                  filling on),
                 (u'n_salinity', <type 'netCDF4._netCDF4.Variable'>
                  int8 n_salinity(cast, pressure)
                      _FillValue: 25
                  unlimited dimensions: 
                  current shape = (1038, 500)
                  filling on)])




```python
time_step = 0.25 # days
pressure_step = 1. # meters
pressure = nc.variables['pressure'][:]
time = np.arange(np.ceil(np.nanmin(nc.variables['time'][:])/time_step)*time_step,np.ceil(np.nanmax(nc.variables['time'][:])/time_step)*time_step,time_step)
```


```python
# set up variables
latitude = np.nan*np.empty(shape=(pressure.shape[0],time.shape[0]))
longitude = np.nan*np.empty(shape=(pressure.shape[0],time.shape[0]))
latitude_fm = np.nan*np.empty(shape=(pressure.shape[0],time.shape[0]))
longitude_fm = np.nan*np.empty(shape=(pressure.shape[0],time.shape[0]))
temperature = np.nan*np.empty(shape=(pressure.shape[0],time.shape[0]))
salinity = np.nan*np.empty(shape=(pressure.shape[0],time.shape[0]))
density = np.nan*np.empty(shape=(pressure.shape[0],time.shape[0]))
pot_temperature = np.nan*np.empty(shape=(pressure.shape[0],time.shape[0]))
pot_density = np.nan*np.empty(shape=(pressure.shape[0],time.shape[0]))
```


```python
nc_time = nc.variables['time'][:].reshape(-1,)
nc_pressure = np.tile(nc.variables['pressure'][:],[1,nc.variables['divenum'][:].size]).reshape(-1,)
nc_latitude = nc.variables['latitude'][:].reshape(-1,)
nc_longitude = nc.variables['longitude'][:].reshape(-1,)
nc_latitude_fm = nc.variables['latitude_flightmodel'][:].reshape(-1,)
nc_longitude_fm = nc.variables['longitude_flightmodel'][:].reshape(-1,)
nc_temperature = nc.variables['temperature'][:].reshape(-1,)
nc_salinity = nc.variables['salinity'][:].reshape(-1,)
nc_density = nc.variables['density'][:].reshape(-1,)
nc_pot_temperature = nc.variables['pot_temperature'][:].reshape(-1,)
nc_pot_density = nc.variables['pot_density'][:].reshape(-1,)
```


```python
for i in range(pressure.shape[0]):
    for j in range(time.shape[0]):
        indices = np.logical_and(np.abs(nc_pressure-pressure[i])<pressure_step/2,np.abs(nc_time-time[j])<time_step/2)
        if indices.size==0:
            continue
        latitude[i,j] = nanmedian(nc_latitude[indices])
        longitude[i,j] = nanmedian(nc_longitude[indices])
        latitude_fm[i,j] = nanmedian(nc_latitude_fm[indices])
        longitude_fm[i,j] = nanmedian(nc_longitude_fm[indices])
        temperature[i,j] = nanmedian(nc_temperature[indices])
        salinity[i,j] = nanmedian(nc_salinity[indices])
        density[i,j] = nanmedian(nc_density[indices])
        pot_temperature[i,j] = nanmedian(nc_pot_temperature[indices])
        pot_density[i,j] = nanmedian(nc_pot_density[indices])       
```

    /usr/local/lib/python2.7/site-packages/ipykernel/__main__.py:3: RuntimeWarning: invalid value encountered in less
      app.launch_new_instance()
    /Users/Zach/anaconda2/lib/python2.7/site-packages/numpy/lib/function_base.py:3858: RuntimeWarning: All-NaN slice encountered
      r = func(a, **kwargs)
    /Users/Zach/anaconda2/lib/python2.7/site-packages/numpy/lib/nanfunctions.py:986: RuntimeWarning: Mean of empty slice
      return np.nanmean(a, axis, out=out, keepdims=keepdims)



```python
for i in range(pressure.size):
    bad_data = np.isnan(latitude[i,:])
    if np.logical_and(np.any(bad_data),np.any(~bad_data)):
        latitude[i,bad_data] = np.interp(time[bad_data],time[~bad_data],latitude[i,~bad_data],left=np.nan,right=np.nan)
        
    bad_data = np.isnan(longitude[i,:])
    if np.logical_and(np.any(bad_data),np.any(~bad_data)):
        longitude[i,bad_data] = np.interp(time[bad_data],time[~bad_data],longitude[i,~bad_data],left=np.nan,right=np.nan)
        
    bad_data = np.isnan(latitude_fm[i,:])
    if np.logical_and(np.any(bad_data),np.any(~bad_data)):
        latitude_fm[i,bad_data] = np.interp(time[bad_data],time[~bad_data],latitude_fm[i,~bad_data],left=np.nan,right=np.nan)
        
    bad_data = np.isnan(longitude_fm[i,:])
    if np.logical_and(np.any(bad_data),np.any(~bad_data)):
        longitude_fm[i,bad_data] = np.interp(time[bad_data],time[~bad_data],longitude_fm[i,~bad_data],left=np.nan,right=np.nan)
```


```python
sigma = time_step
DAC_U = nc.variables['DAC_U'][:]
DAC_V = nc.variables['DAC_V'][:]
DAC_time = np.nanmean(nc.variables['time'][:],axis=1)
DAC_U_interp = np.zeros(shape=time.shape)
DAC_V_interp = np.zeros(shape=time.shape)
for i in range(time.size):
    DAC_U_interp[i] = np.nansum(np.exp(-((DAC_time-time[i])/sigma)**2/2)*DAC_U)/np.nansum(np.exp(-((DAC_time-time[i])/sigma)**2/2))
    DAC_V_interp[i] = np.nansum(np.exp(-((DAC_time-time[i])/sigma)**2/2)*DAC_V)/np.nansum(np.exp(-((DAC_time-time[i])/sigma)**2/2))
```


```python
FILL_VALUE = -999
nc = Dataset('sg579_JanApr_L3_2m.nc','w')
nc.title = 'Level 3 Glider data'
nc.mission = 'OSMOSIS (Ocean Surface Mixing, Ocean Sub-mesoscale Interaction Study)'
nc.institution = 'Caltech'
nc.glider = 'SG579'
nc.author = 'Zachary K Erickson'
nc.contact = 'zerickso@caltech.edu'
nc.history = 'Created '+ python_time.ctime(python_time.time())
nc.Conventions = 'CF-1.6'

nc.createDimension('time',time.size)
nc.createDimension('pressure',pressure.size)
nc_time = nc.createVariable('time',np.dtype('float64').char,('time'))
nc_time.units = 'days since January 1, 2013'
nc_time.scale_factor = 1
nc_time.add_offset = 0.
nc_pressure = nc.createVariable('pressure',np.dtype('float64').char,('pressure'))
nc_pressure.units = 'db'
nc_pressure.scale_factor = 1.
nc_pressure.add_offset = 0.
nc_DAC_U = nc.createVariable('DAC_U',np.dtype('float64').char,('time',),fill_value=-999.)
nc_DAC_U.units = 'm/s'
nc_DAC_U.scale_factor = 1.
nc_DAC_U.add_offset = 0.
nc_DAC_U.sigma_Gaussian_interpolation = time_step
nc_DAC_V = nc.createVariable('DAC_V',np.dtype('float64').char,('time',),fill_value=-999.)
nc_DAC_V.units = 'm/s'
nc_DAC_V.scale_factor = 1.
nc_DAC_V.add_offset = 0.
nc_DAC_V.sigma_Gaussian_interpolation = time_step
nc_lat = nc.createVariable('latitude',np.dtype('float64').char,('pressure','time'),fill_value=FILL_VALUE)
nc_lat.units = 'degrees north'
nc_lat.scale_factor = 1.
nc_lat.add_offset = 0.
nc_lon = nc.createVariable('longitude',np.dtype('float64').char,('pressure','time'),fill_value=FILL_VALUE)
nc_lon.units = 'degrees east'
nc_lon.scale_factor = 1.
nc_lon.add_offset = 0.
nc_lat_fm = nc.createVariable('latitude_flightmodel',np.dtype('float64').char,('pressure','time'),fill_value=FILL_VALUE)
nc_lat_fm.units = 'degrees north'
nc_lat_fm.scale_factor = 1.
nc_lat_fm.add_offset = 0.
nc_lon_fm = nc.createVariable('longitude_flightmodel',np.dtype('float64').char,('pressure','time'),fill_value=FILL_VALUE)
nc_lon_fm.units = 'degrees east'
nc_lon_fm.scale_factor = 1.
nc_lon_fm.add_offset = 0.
nc_temperature = nc.createVariable('temperature',np.dtype('float64').char,('pressure','time'),fill_value=FILL_VALUE)
nc_temperature.units = 'deg C'
nc_temperature.scale_factor = 1
nc_temperature.add_offset = 0.
nc_salinity = nc.createVariable('salinity',np.dtype('float64').char,('pressure','time'),fill_value=FILL_VALUE)
nc_salinity.units = 'Practical Salinity Units'
nc_salinity.scale_factor = 1
nc_salinity.add_offset = 0.
nc_density = nc.createVariable('density',np.dtype('float64').char,('pressure','time'),fill_value=FILL_VALUE)
nc_density.units = 'kg/m^3'
nc_density.scale_factor = 1
nc_density.add_offset = 1000.
nc_pot_temperature = nc.createVariable('pot_temperature',np.dtype('float64').char,('pressure','time'),fill_value=FILL_VALUE)
nc_pot_temperature.units = 'deg C'
nc_pot_temperature.scale_factor = 1
nc_pot_temperature.add_offset = 0.
nc_pot_density = nc.createVariable('pot_density',np.dtype('float64').char,('pressure','time'),fill_value=FILL_VALUE)
nc_pot_density.units = 'kg/m^3'
nc_pot_density.scale_factor = 1
nc_pot_density.add_offset = 1000.

nc_pressure[:] = pressure;
nc_time[:] = time
nc_DAC_U[:] = DAC_U_interp
nc_DAC_V[:] = DAC_V_interp
nc_lat[:] = latitude
nc_lon[:] = longitude
nc_lat_fm[:] = latitude_fm;
nc_lon_fm[:] = longitude_fm;
nc_temperature[:] = temperature
nc_salinity[:] = salinity
nc_density[:] = density
nc_pot_temperature[:] = pot_temperature
nc_pot_density[:] = pot_density
```


```python
nc.close()
```
