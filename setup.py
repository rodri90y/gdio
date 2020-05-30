from setuptools import setup

setup(
    name='gdio',
    version='0.0.1',
    description='Gridded data io library',
    long_description='Gridded data io library. Grib and netcdf files handler, ',
    license='MIT',
    packages=['gdio'],
    author='Rodrigo Yamamoto',
    author_email='codes@rodrigoyamamoto.com',
    keywords=['gdio','grib', 'netcdf'],
    url='https://github.com/bgse/hellostackoverflow',
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
    install_requires=['numpy'],
)