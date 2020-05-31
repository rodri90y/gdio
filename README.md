
# GDIO - Gridded Data IO

A simple and concise gridded data IO library to read multiples grib and netcdf files, automatic spatial interpolation of the all data to a single resolution.

The library gdio is based on my own professionals and personal needs as a meteorologist. The commons libraries always fail when you need to read handle multiples large netcdf/grib files, with differents resolutions and timesteps.
## Instalation
```
conda install gdio
```

#### Required dependencies

+ Python (3.6 or later)
+ numpy (1.18.4 or later)
+ pygrib (2.0.4 or later)
+ netCDF4 (1.5.1.2 or later)
+ Texttable (1.6.2 or later)

#### Optional dependencies
+ scipy (1.4.1 or later)

#### Testing
```
python -m unittest 
```


## Reading files
The gdio support the IO of grib1/2 and netcdf file formats, allowing the time and spatial subdomains cut.

#### Netcdf
The class netcdf encapsulate all netcdf functions, as well as the cutting of time and spatial domains, returning the netcdf data as a dictionary type.

Simple reading
```
from gdio.netcdf import netcdf
nc = netcdf(verbose=False)

ds = nc.nc_load('data/era5_20191227_lev.nc')
```

Reading a subsample in time (time 12-24) and space (bbox -30,-60 and 10,-40)

```
ds = nc.nc_load('data/era5_20191227_lev.nc', cut_domain=(-30, -60, 10, -40), cut_time=(12, 24))
```


#### Grib
The class netcdf encapsulate all netcdf functions, as well as the cutting of time and spatial domains , returning the netcdf data as a dictionary type.

Simple reading
```
from gdio.grib import grib
gr = grib(verbose=False)
ds = gr.gb_load('data/era5_20191226-27_lev.grib')
```

Reading a subsample in time (time 12-24) and space (bbox -30,-60 and 10,-40)

```
ds = gr.gb_load('data/era5_20191226-27_lev.grib', cut_domain=(-30, -60, 10, -40), cut_time=(12, 24))
```

Listing the grib content

```
gr.gb_info(ifile)
```

### Reading multiple files
This class has high level routines for multiple files and type reading, returning the netcdf/grib data as a list of dictionary type.

```
from gdio.core import gdio

ds = gdio(verbose=False)
ds.mload(['tests/data/era5_20191226-27_lev.grib', 'tests/data/era5_20191227_lev.nc'],  
        merge_files=True, uniformize_grid=True, cut_domain=(-30, 300, 10, 320), inplace=True)
```
## Release History


## Meta
Rodrigo Yamamoto codes@rodrigoyamamoto.com

https://github.com/rodri90y/gdio

## Contributing
* 0.0.2
    * alpha release
    

## License

MIT
