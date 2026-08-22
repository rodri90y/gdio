
# Anaconda

### Install
```commandline
conda install anaconda-client
anaconda login

conda install conda-build
conda config --set anaconda_upload no
```
### Build
```
export PYTHON=$(which python)

conda build . (in the project root dir)
conda build name_project
conda build --override-channels -c conda-forge --python 3.11 .
```
em caso de erro:
```
conda clean -a

```
### Install locally
conda install --use-local

anaconda upload /home/user/anaconda3/conda-bld/linux-64/gdio-x.x.x-py37h39e3cac_0.tar.bz2 --force

### Clear build
```
conda build purge
conda build purge-all
```
----
# PIP register
```
PyPi: https://pypi.org/account/register/
Required Tools

sudo python -m pip install --upgrade pip setuptools wheel
sudo python -m pip install tqdm
sudo python -m pip install --user --upgrade twine

or

conda install twine
```
## Setup the project

### Create a setup file setup.py in the package directory
#### Install on Your Local Machine
```commandline
from setuptools import setup, find_packages

	
with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setup(
	name='gdio',
	version='0.0.4',
	description='Gridded data io library',
    	long_description=long_description,
    	long_description_content_type='text/markdown',
	license='MIT',
download_url="https://github.com/rodri90y/gdio/archive/master.zip",
	packages=find_packages('gdio'),
	author='Rodrigo Yamamoto',
	author_email='codes@rodrigoyamamoto.com',
	keywords=['gdio','grib','netcdf'],
	url='https://github.com/rodri90y/gdio',
	classifiers=[
    	"Programming Language :: Python :: 3",
    	"License :: OSI Approved :: MIT License",
    	"Operating System :: OS Independent",
	],
	python_requires='>=3.6',
	install_requires=['numpy','netCDF4','pygrib','texttable'],
)
```
### Install on Your Local Machine
Add a LICENSE

### Compiling the package
```
python setup.py bdist_wheel
```
will create a  structure: build/, dist/ and project.egg.info/

### Install on Your Local Machine
```commandline
python -m pip install dist/gdio-0.1-py3-none-any.whl
```


### Upload on pip

Create on home directory the pypirc file, which stores the PyPi repository information
Windows :  C:\Users\UserName\.pypirc
linux :   ~/.pypirc	
```
[pypi]
repository = https://upload.pypi.org/legacy/
username = user
[pypitest]
repository=https://testpypi.python.org/pypi
username = user
``

To upload your dist/*.whl file on PyPi, use Twine:
``
python -m twine upload dist/*

or

python -m twine upload --repository testpypi dist/*
``
