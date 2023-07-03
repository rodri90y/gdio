from setuptools import setup, find_packages

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setup(
    name='gdio',
    version='0.3.3',
    description='Gridded data io library',
    long_description=long_description,
    long_description_content_type='text/markdown',
    url='https://github.com/rodri90y/gdio',
    download_url="https://github.com/rodri90y/gdio/archive/v0.3.3.tar.gz",
    license='MIT',
    packages=find_packages(),
    author='Rodrigo Yamamoto',
    author_email='codes@rodrigoyamamoto.com',
    keywords=['gdio', 'grib', 'netcdf', 'hdf5'],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.8.5',
    install_requires=['numpy', 'netCDF4', 'h5py', 'eccodes', 'pyproj'],
)
