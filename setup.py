from setuptools import setup, find_packages

setup(
    name='gdio',
    version='0.1.5',
    description='Gridded data io library',
    long_description='Gridded data io library. Grib and netcdf files handler, ',
    long_description_content_type='text/markdown',
    url='https://github.com/rodri90y/gdio',
    download_url="https://github.com/rodri90y/gdio/archive/v0.1.5.tar.gz",
    license='MIT',
    packages=find_packages(),
    author='Rodrigo Yamamoto',
    author_email='codes@rodrigoyamamoto.com',
    keywords=['gdio','grib','netcdf'],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
    install_requires=['numpy','netCDF4','texttable','pygrib'],
)