
# GDIO - Gridded Data IO

A simple and concise gridded data IO library for reading multiples grib and netcdf files, automatic spatial interpolation of the all data to a single resolution.

The library gdio is based on my own professionals and personal needs as a meteorologist. 
The currents libraries always fail when you need to read handle multiples large 
netcdf/grib files, with different resolutions and time steps.

After version 0.1.2 the output data was converted to object with key-values accessible using attribute notation, and after version 0.1.8 a new multilevel dictionary data structure. 
In the version 0.2.5 the latitude and longitude come in mesh array (ny,nx) format to support irregular or lambert projection.

## Instalation
```
conda config --env --add channels conda-forge
conda install -c rodri90y gdio

if you are using pip install, before install manually the requirements

conda create -n envname --file requirements/base.txt
pip install gdio
or
pip install --index-url https://test.pypi.org/simple/ --upgrade --no-cache-dir --extra-index-url=https://pypi.org/simple/ gdio
```

#### Required dependencies

conda config --add channels conda-forge

+ Python (3.7.9=> or later)
+ netCDF4 (1.5.8 or later)
+ eccodes (2.24.2 or later)
+ python-eccodes (1.4.0 or later)
+ pyproj


#### Optional dependencies
+ scipy (1.4.1 or later)

#### Testing
```
python -m unittest 
```


## Reading files
The gdio support the IO of grib1/2 and netcdf file formats, allowing the time and spatial subdomains cut.

This library unifies categories of information (variable, level, members) in a single 
data structure as a multilevel dictionary/attribute, regardless of the format read (netcdf and grib), the 
output format will be standardized in order to simplify access to the data.

In the dataset first level the following parameters are accessible: ref_time, time_units and time in addition to the weather variables.
ds.ref_time, ds.time
At the variable level we have: level_type, param_id, long_name, parameterUnits, latitude and longitude and at vertical level (isobaricInh, surface, etc) the variable data as value and level are exposed.

Structure data:

    + dataset
        + ref_time
        + time_units
        + time
        + variable (u,v,2t,etc) 
            + level_type
            + param_id
            + long_name
            + parameter_units
            + latitude
            + longitude
            + isobaricInhPa/surface/maxWind/sigma (any level type key)
                + value
                + level
                + members
                + grid_type
                + projparams

Example:
            
    ds.time
    ds.time_units
    ds.v.latitude
    ds.v.isobaricInhPa.value
    ds.v.isobaricInhPa.level
    ds.v.isobaricInhPa.members


### Reading multiple files
This class has high level routines for multiple files and type reading, returning the netcdf/grib data as a list of dictionary type.

```
from gdio.core import gdio

ds = gdio(verbose=False)
ds.mload(['tests/data/era5_20191226-27_lev.grib', 'tests/data/era5_20191227_lev.nc'],  
        merge_files=True, uniformize_grid=True, inplace=True)

>>> ds.dataset[0].keys()
dict_keys(['ref_time', 'time_units', 'time', 'longitude', 'latitude', 't', 'u', 'v', 'r'])

>>> print(ds.dataset[0].u.isobaricInhPa.value.shape)
(1, 6, 7, 241, 281)

>>> ds.dataset[0].time
masked_array(data=[datetime.datetime(2019, 12, 26, 0, 0),
                   datetime.datetime(2019, 12, 26, 12, 0),
                   datetime.datetime(2019, 12, 27, 0, 0),
                   datetime.datetime(2019, 12, 27, 12, 0),
                   datetime.datetime(2019, 12, 27, 0, 0),
                   datetime.datetime(2019, 12, 27, 12, 0)],
             mask=False,
       fill_value='?',
            dtype=object)

```
Loading the data into the spatial subdomain between lat -30, lon 300 and lat 10, lon 320, selecting the time between 
timespteps 12 and 24, and changing the variable names t and u to 2t and 10u.

```
from gdio.core import gdio

ds = gdio(verbose=False)
ds.mload(['tests/data/era5_20191226-27_lev.grib', 'tests/data/era5_20191227_lev.nc'],  
        merge_files=True, uniformize_grid=True, 
        cut_domain=(-30, 300, 10, 320), cut_time=(12, 24), 
        rename_vars={'t': '2t', 'u': '10u'}, inplace=True)

>>> ds.dataset[0].keys()
dict_keys(['ref_time', 'time_units', 'time', 'longitude', 'latitude', 'r', '2t', '10u', 'v'])

>>> print(ds.dataset[0]['10u'].isobaricInhPa.value.shape)
(1, 2, 7, 160, 80)

>>> ds.dataset[0].time
masked_array(data=[datetime.datetime(2019, 12, 26, 12, 0),
                   datetime.datetime(2019, 12, 27, 0, 0),
                   datetime.datetime(2019, 12, 27, 12, 0)],
             mask=False,
       fill_value='?',
            dtype=object)

```

The following parameters can be set to operate on the data during reading.

**uniformize_grid:     boolean**\
interpolate all gridded data to first grid data file resolution

**vars:                list**\
variables names

**merge_files:         boolean**\
merge the variables data of all files into a single data array per variable

**cut_time:            tuple**\
range of time to cut ex.: (0,10)/(0,None)/(None,10)

**cut_domain:          tuple**\
range of latitudes and longitudes to cut: (lat1, lon1, lat2, lon2)
ex.: (-45,-90,20,-30)/(-45,None,20,-30)/(None,-90,None,-20)

**level_type:          list**\
type of level (hybrid, isobaricInhPa, surface)

**filter_by:           dictonary**\
dict with grib parameters at form of pair key:values (list or single values)
eg: filter_by={'perturbationNumber': [0,10],'level': [1000,500,250]} or filter_by={'gridType': 'regular_ll'}|
Obs: this parameter only works on grib files

**rename_vars:         dictonary**\
rename the original variable name (key) to a new name (value). 

Eg. {'tmpmdl': 't', 'tmpprs': 't'}

**sort_before:         bool**\
Sort fields before process validityDate, validityTime, paramId, typeOfLevel, perturbationNumber and level. Warning high
consumption of memory, just use when the grib data structure is not standard


### Selecting a sub sample in mload dataset
Select data by coordinates (date, latitude, longitude, levels and members)

```
sub_set = ds.sel(dates=[datetime(2019,12,26,12,0)], latitude=[-23.54,-22], longitude=[-46.64,-42.2], level=[2,6])

>>> print(sub_set[0].get('u').isobaricInhPa.value.shape)
(1, 1, 4, 6, 18)
```

### Showing the data structure
Prints the data structure tree.
```
>>> ds.describe

    +-- ref_time: 2019-12-26 00:00:00
    +-- time_units: hours
    +-- time: <class 'numpy.ma.core.MaskedArray'> (6,)
    +-- r 
        +-- isobaricInhPa 
            +-- value: <class 'numpy.ndarray'> (1, 6, 7, 160, 80)
            +-- level: [200, 300, 500, 700, 800, 950, 1000]
            +-- members: [0]
        +-- centre: 'ecmwf',
        +-- dataType: 'an',
        +-- param_id: 157
        +-- long_name: Relative humidity
        +-- parameter_units: %
        +-- latitude: <class 'numpy.ndarray'> (160, 80)
        +-- longitude: <class 'numpy.ndarray'> (160, 80)
        +-- level_type: ['isobaricInhPa']
        +-- grid_type: 'regular_ll'
        +-- projparams: { 'a': 6371229.0, 'b': 6371229.0, 'proj': 'regular_ll'}
        
    .
    .
    .
    
    +-- v 
    +-- isobaricInhPa 
        +-- value: <class 'numpy.ndarray'> (1, 6, 7, 160, 80)
        +-- level: [200, 300, 500, 700, 800, 950, 1000]
        +-- members: [0]
    +-- centre: 'ecmwf',
    +-- dataType: 'an',
    +-- param_id: 132
    +-- long_name: V component of wind
    +-- parameter_units: m s**-1
    +-- latitude: <class 'numpy.ndarray'> (160, 80)
    +-- longitude: <class 'numpy.ndarray'> (160, 80)
    +-- level_type: ['isobaricInhPa']
    +-- grid_type: 'regular_ll'
    +-- projparams: { 'a': 6371229.0, 'b': 6371229.0, 'proj': 'regular_ll'}
```


Setting the ensemble grouping grib id key

```
ds.fields_ensemble = 'perturbationNumber'
ds.fields_ensemble_exception = [0]
```


#### Grib
The class netcdf encapsulates all grib functions, as well as the cutting of time and spatial domains , returning the netcdf data as a dictionary type.

Simple reading
```
from gdio.grib import grib
gr = grib(verbose=False)
ds = gr.gb_load('data/era5_20191226-27_lev.grib')

>>> ds.keys()
dict_keys(['ref_time', 'time_units', 'time', 't', 'u', 'v', 'r'])
>>> print(ds.u.isobaricInhPa.value.shape)
(1, 4, 7, 241, 281)
>>> print(ds.u.level_type)
['isobaricInhPa']
>>> print(ds.u.keys())
dict_keys(['centre', 'dataType','isobaricInhPa', 'param_id', 'long_name', 'parameter_units', 'latitude', 'longitude', 'level_type', 'grid_type', projparams])
>>> print(ds.u.isobaricInhPa.level)
[200, 300, 500, 700, 800, 950, 1000]
>>> print(ds.u.parameter_units)
m s**-1
>>> print(ds.u.param_id)
131
```

Reading a subsample in time (time 12-24) and space (bbox -30,-60 and 10,-40)

```
ds = gr.gb_load('data/era5_20191226-27_lev.grib', cut_domain=(-30, -60, 10, -40), cut_time=(12, 24))
```

Setting the ensemble grouping grib id key

```
gr.fields_ensemble = 'perturbationNumber'
gr.fields_ensemble_exception = [0]
```

Filtering by a grib key, dict with grib parameters at form of pair key:
values (list or single values)
eg: filter_by={'perturbationNumber': [0,10],'level': [1000,500,250]}
                            or filter_by={'gridType': 'regular_ll'}

```
ds = gr.gb_load('tests/data/era5_20191226-27_lev.grib', 
                cut_domain=(-30, -60, 10, -40), 
                cut_time=(12, 24), 
                filter_by={'perturbationNumber': 0, 'level':[200,500,950]})
>>> print(ds.u.isobaricInhPa.level)
[200, 500, 950]
```
Rename variables
A dictionary input will rename variables names (key) for a new name (value).
Eg. {'tmpmdl': 't', 'tmpprs': 't'}

```
ds = gr.gb_load('data/era5_20191227_lev.nc', rename_vars={'u':'10u'})
>>> ds.keys()
dict_keys(['ref_time', 'time_units', 'time', 't', '10u', 'v', 'r'])
```
Sorting grib parameter before (extra consumption of memory and possible a little slow). 
Fix grib files unstructured or non-standard.
```
ds = gr.gb_load('data/era5_20191227_lev.nc', sort_before=True)
```
#### Writing a netcdf file

From the loaded dataset
```
nc.nc_write('data/output.nc', ds)
```
From a dictionary
```
from gdio.grib import grib
gr = grib(verbose=False)
ds = gr.gb_load('data/era5_20191226-27_lev.grib')

gr.gb_write('output.grib', self.gbr, least_significant_digit=3, packingType='grid_jpeg')
```

#### Netcdf
The class netcdf encapsulates all netcdf functions of reading and writing, as well as the cutting of time and spatial domains, returning the netcdf data as a dictionary type. The returned dictionary contains for each variable the value, param_id, type_level, level and parameter_units property.

Simple reading
```
from gdio.netcdf import netcdf
nc = netcdf(verbose=False)

ds = nc.nc_load('tests/data/era5_20191227_lev.nc')
>>> ds.keys()
dict_keys(['ref_time', 'time_units', 'time', 'r', 't', 'u', 'v'])
>>> print(ds.u.isobaricInhPa.value.shape)
(1, 2, 7, 161, 241)
>>> print(ds.u.level_type)
['isobaricInhPa']
>>> print(ds.u.keys())
dict_keys(['isobaricInhPa', 'param_id', 'long_name', 'parameter_units', 'latitude', 'longitude', 'level_type'])
>>> print(ds.u.isobaricInhPa.level)
[200, 300, 500, 700, 800, 950, 1000]
>>> print(ds.u.parameter_units)
m s**-1
>>> print(ds.u.param_id)
None
```

Reading a subsample in time (time 12-24) and space (bbox -30,-60 and 10,-40). The returned multilevels dictionary/attributes contains for each variable the value, param_id, type_level, level and parameter_units property.

```
ds = nc.nc_load('data/era5_20191227_lev.nc', cut_domain=(-30, -60, 10, -40), cut_time=(12, 24))
>>> print(ds.u.isobaricInhPa.value.shape)
(1, 1, 7, 80, 40)
```
Rename variables
A dictionary input will rename variables names (key) for a new name (value).
Eg. {'tmpmdl': 't', 'tmpprs': 't'}

```
ds = nc.nc_load('data/era5_20191227_lev.nc', rename_vars={'u':'10u'})
>>> ds.keys()
dict_keys(['ref_time', 'time_units', 'time', 't', '10u', 'v', 'r'])
```

#### Writing a netcdf file

From the loaded dataset
```
nc.nc_write('data/output.nc', ds)
```
From a dictionary
```
from datetime import datetime
import numpy as np
from gdio.netcdf import netcdf

nc = netcdf(verbose=False)

ds = {'ref_time': datetime(2019, 12, 27, 0, 0), 
      'time_units': 'hours', 
      'time': np.array([12]),
      'u': {'isobaricInhPa': {  'value': np.random.random((1, 1, 7, 80, 40)),
                                'level': [200, 300, 500, 700, 800, 950, 1000]
                              },
            'param_id': None, 
            'long_name': 'U component of wind', 
            'level_type': ['isobaricInhPa'],
            'parameter_units': 'm s**-1',
            'longitude': np.array([300. , 300.5, 301. , 301.5, 302. , 302.5, 303. , 303.5,
               304. , 304.5, 305. , 305.5, 306. , 306.5, 307. , 307.5,
               308. , 308.5, 309. , 309.5, 310. , 310.5, 311. , 311.5,
               312. , 312.5, 313. , 313.5, 314. , 314.5, 315. , 315.5,
               316. , 316.5, 317. , 317.5, 318. , 318.5, 319. , 319.5]),
            'latitude': np.array([-30. , -29.5, -29. , -28.5, -28. , -27.5, -27. , -26.5,
               -26. , -25.5, -25. , -24.5, -24. , -23.5, -23. , -22.5,
               -22. , -21.5, -21. , -20.5, -20. , -19.5, -19. , -18.5,
               -18. , -17.5, -17. , -16.5, -16. , -15.5, -15. , -14.5,
               -14. , -13.5, -13. , -12.5, -12. , -11.5, -11. , -10.5,
               -10. ,  -9.5,  -9. ,  -8.5,  -8. ,  -7.5,  -7. ,  -6.5,
                -6. ,  -5.5,  -5. ,  -4.5,  -4. ,  -3.5,  -3. ,  -2.5,
                -2. ,  -1.5,  -1. ,  -0.5,   0. ,   0.5,   1. ,   1.5,
                 2. ,   2.5,   3. ,   3.5,   4. ,   4.5,   5. ,   5.5,
                 6. ,   6.5,   7. ,   7.5,   8. ,   8.5,   9. ,   9.5]),
            }
      }

nc.nc_write('data/output.nc', ds)
```


## Routines
### gdio.mload
Load multiple files (netcdf/grib) returning the data as a list of dictionary type interpolating the data to a same grid

```
mload(files, vars=None, merge_files=True, cut_time=None,
      cut_domain=None, level_type=None, filter_by={},
      uniformize_grid=True, sort_before=False, inplace=False)
```          
**files:               list**

files names
                                    
**uniformize_grid:     boolean**\
interpolate all ncs to first nc grid specification

**vars:                list**\
variables names

**merge_files:         boolean**\
merge files

**cut_time:            tuple**\
                        range of time to cut ex.: (0,10)/(0,None)/(None,10)

**cut_domain:          tuple**\
                        range of latitudes and longitudes to cut: (lat1, lon1, lat2, lon2)
                        ex.: (-45,-90,20,-30)/(-45,None,20,-30)/(None,-90,None,-20)

**level_type:          list**\
                        type of level (hybrid, isobaricInhPa, surface)

**filter_by:           dictonary**\
dict with grib parameters at form of pair key:values (list or single values)
eg: filter_by={'perturbationNumber': [0,10],'level': [1000,500,250]} or filter_by={'gridType': 'regular_ll'}|

**rename_vars:         dictonary**\
rename variables names (key) for a new name (value). Eg. {'tmpmdl': 't', 'tmpprs': 't'}

**sort_before:         bool**\
Sort fields before process validityDate, validityTime, paramId, typeOfLevel, 
perturbationNumber and level. Warning extra consumption of memory and time, 
just use when the grib data structure is not standard (default False)

**return:**                    list of dictionaries

### gdio.sel
Select data by coordinates (date, latitude, longitude, levels and members)

```
sel(data=None, latitude=None, longitude=None, 
    dates=None, level=None, member=None, date_format="%Y-%m-%d %H:%M")
```


**data:       list of dictionaries**\
                             raw dataset
                             
**latitude:     list of floats**\
                             latitudes
                             range of latitudes to select: [lat1, lat2]
                             especific latitudes (1 or >2) [lat1, lat2, lat3, ...]

**longitude:    list of floats**\
                             range of longitudes to select: [lon1, lon2]
                             especific longitudes (1 or >2) [lon1, lon2, lon3, ...]

**dates:        list of datetime/string**\
                             datetime/string date
                             range of dates to select: [date1, date2]
                             especific dates (1 or >2) [date1, date2, date3, ...]

**level:        list of int**\
                             range of levels to select: [level1, level2]
                             especific levels (1 or >2) [level1, level2, level3, ...]

**member:       list of int**\
                             range of levels to select: [member, member]
                             especific levels (1 or >2) [level1, level2, level3, ...]

**return:**     list of dictionaries

### gdio.grib.gb_load
Load grib file
```
def gb_load(ifile, vars=None, level_type=None,
            cut_time=None, cut_domain=None, filter_by={},
            rename_vars={}, sort_before=False)
```

**ifile:       string**\
                            grib 1 or 2 file name

**vars:        list**\
                            variables short name or id parameter number

**cut_time:    tuple**\
                            range of time to cut ex.: (0,10)/(0,None)/(None,10)

**cut_domain:  tuple**\
                            range of latitudes and longitudes to cut: (lat1, lon1, lat2, lon2)
                            ex.: (-45,290,20,330)/(-45,None,20,330)/(None,290,None,320)

**level_type:  list**\
                            type of level (hybrid, isobaricInhPa, surface)

**filter_by:   dictonary**\
                            dict with grib parameters at form of pair key:values (list or single values)
                            eg: filter_by={"perturbationNumber": [0,10],"level": [1000,500,250]}
                            or filter_by={"gridType": "regular_ll"}

**rename_vars: dictonary**\
                            rename variables names (key) for a new name (value).
                            Eg. {"tmpmdl": "t", "tmpprs": "t"}

**sort_before: bool**\
                            Sort fields before process validityDate, validityTime, paramId, typeOfLevel, perturbationNumber and level
                            Warning high consumption of memory, just use when the grib data structure is not standard

**return: dictonary/attributes**\
multiple time data container

### gdio.grib.gb_write
Write grib2 file
```
def gb_write(ofile, data, packingType='grid_simple', least_significant_digit=3, **kwargs))
```
**ifile: string**\
file path

**data: dict**\
dataset

**packingType: string**\
packingType\
- Type of packing:
  - grid_simple
  - spectral_simple
  - grid_simple_matrix
  - grid_jpeg
  - grid_png
  - grid_ieee
  - grid_simple_log_preprocessing
  - grid_second_order

**least_significant_digit: int (default 3)**\
Specify the power of ten of the smallest decimal place in the data that is a
reliable value that dramatically improve the compression by quantizing
(or truncating) the data

### gdio.netcdf.nc_load
Load netcdf files
```
nc_load(ifile, vars=None, cut_time=None, cut_domain=None, level_type=None, rename_vars={}):
```

**ifile:       string**\
                    netcdf file name
                    
**vars:        list**\
                    variables short name
                    
**cut_time:    tuple**\
                    range of time (absolute) to cut ex.: (0,10)/(0,None)/(None,10)
                    
**cut_domain:  tuple**\
                    range of latitudes and longitudes to cut: (lat1, lon1, lat2, lon2)
                    ex.: (-45,290,20,330)/(-45,None,20,330)/(None,290,None,320)
                    
**level_type:  list**\
                    type of level (hybrid, isobaricInhPa, surface)

**rename_vars: dictonary**\
                            rename variables names (key) for a new name (value).
                            Eg. {"tmpmdl": "t", "tmpprs": "t"}
                            
**return: dictonary/attributes**\
multiple time data container

### gdio.netcdf.nc_write
Write netcdf file
```
nc_write(ifile, data, zlib=True, netcdf_format='NETCDF4')
```



**ifile:           string**\
                                file path
                                
**data:            dict**\
                                dataset
                                
**zlib:            bool**\
                                enable compression
                                
**netcdf_format:   string**\
                                netcdf format: NETCDF4, NETCDF4_CLASSIC, NETCDF3_CLASSIC or NETCDF3_64BIT


**complevel:      int**\
 compression level (default 4)

**least_significant_digit: int**\
specify the power of ten of the smallest decimal place in the data that is a
                reliable value that dramatically improve zlib compression by quantizing
                (or truncating) the data (default None)
### gdio.remapbil
```
remapbil(data, lon, lat, lon_new, lat_new, order=1, masked=False)
```

Interpolate data to new domain resolution

**data: array**\
                        3D data (time,lon,lat)

**lon: array**

**lat: array**

**lon_new: array**\
                        new grid logitudes

**lat_new: array**\
                        new grid latitudes

**order:   int**\
                        0- nearest-neighbor, 1 - bilinear, 2 - cubic spline

**masked: boolean**\
                        If True, points outside the range of xin and yin
                        are masked (in a masked array). If masked is set to a number

**return: 3D array**


## Dev utils
Docker compose to support development

### Commands
 - make build
   - Build the container
 - make up
   - Start container
 - make stop
   - Stop container
 - make test
   - Run unit tests in container
 - make bash
   - Access container
 - make ipython
   - Run ipython in container
 - make fix
   - Run autopep to fix code format

## Release History


## Meta
Rodrigo Yamamoto codes@rodrigoyamamoto.com

https://github.com/rodri90y/gdio

## Contributing

* 0.2.6
    * alpha release
    

## License

MIT
