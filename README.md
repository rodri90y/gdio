
# GDIO - Gridded Data IO

A simple and concise gridded data IO library for reading multiples grib and netcdf files, automatic spatial interpolation of the all data to a single resolution.

The library gdio is based on my own professionals and personal needs as a meteorologist. The commons libraries always fail when you need to read handle multiples large netcdf/grib files, with differents resolutions and timesteps.

After version 0.1.2 the output data was converted to object with key-values accessible using attribute notation

## Instalation
```
conda config --env --add channels conda-forge
conda install -c rodri90y gdio

if you are using pip install, before install manually the requirements
```

#### Required dependencies

conda config --add channels conda-forge

+ Python (3.7.* or later)
+ numpy (1.18.4 or later)
+ netCDF4 (1.5.1.2 or later)
+ eccodes (2.12.3 or later)
+ python-eccodes
+ pyproj


#### Optional dependencies
+ scipy (1.4.1 or later)

#### Testing
```
python -m unittest 
```


## Reading files
The gdio support the IO of grib1/2 and netcdf file formats, allowing the time and spatial subdomains cut.

#### Netcdf
The class netcdf encapsulates all netcdf functions of reading and writing, as well as the cutting of time and spatial domains, returning the netcdf data as a dictionary type. The returned dictionary contains for each variable the value, param_id, type_level, level and parameter_units property.

Simple reading
```
from gdio.netcdf import netcdf
nc = netcdf(verbose=False)

ds = nc.nc_load('data/era5_20191227_lev.nc')

>>> print(ds.u.value.shape)
(2, 7, 161, 241)
>>> print(ds.u.type_level)
millibars
>>> print(ds.u.level)
[200, 300, 500, 700, 800, 950, 1000]
>>> print(ds.u.parameter_units)
m s**-1
>>> print(ds.u.param_id)
None
```

Reading a subsample in time (time 12-24) and space (bbox -30,-60 and 10,-40). The returned dictionary contains for each variable the value, param_id, type_level, level and parameter_units property.

```
ds = nc.nc_load('data/era5_20191227_lev.nc', cut_domain=(-30, -60, 10, -40), cut_time=(12, 24))
```

Writing a netcdf file

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

      'u': {'param_id': None, 
            'type_level': 'millibars', 
            'level': [200, 300, 500, 700, 800, 950, 1000], 
            'parameter_units': 'm s**-1',
            'value': np.random.random((1, 7, 80, 40))
            }
      }

nc.nc_write('data/output.nc', ds)
```

#### Grib
The class netcdf encapsulates all grib functions, as well as the cutting of time and spatial domains , returning the netcdf data as a dictionary type.

Simple reading
```
from gdio.grib import grib
gr = grib(verbose=False)
ds = gr.gb_load('data/era5_20191226-27_lev.grib')

>>> print(ds.u.value.shape)
(4, 7, 161, 241)
>>> print(ds.u.type_level)
isobaricInhPa
>>> print(ds.u.level)
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
ds = gr.gb_load('data/era5_20191226-27_lev.grib', 
                cut_domain=(-30, -60, 10, -40), 
                cut_time=(12, 24), 
                filter_by={'perturbationNumber': 0})
```

### Reading multiple files
This class has high level routines for multiple files and type reading, returning the netcdf/grib data as a list of dictionary type.

```
from gdio.core import gdio

ds = gdio(verbose=False)
ds.mload(['tests/data/era5_20191226-27_lev.grib', 'tests/data/era5_20191227_lev.nc'],  
        merge_files=True, uniformize_grid=True, cut_domain=(-30, 300, 10, 320), inplace=True)
```

Setting the ensemble grouping grib id key

```
ds.fields_ensemble = 'perturbationNumber'
ds.fields_ensemble_exception = [0]
```

## Release History


## Meta
Rodrigo Yamamoto codes@rodrigoyamamoto.com

https://github.com/rodri90y/gdio

## Contributing
* 0.1.5
    * alpha release
    

## License

MIT
