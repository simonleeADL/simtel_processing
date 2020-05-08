# simtel_processing
Code used to process CTA telescope simulations made with sim-telarray (.simtel.gz files) in to .hdf5 files with Hillas parameterization etc.

## Install
```git clone github.com/simonleeADL/simtel_processing/```

## Dependencies
[ctapipe](github.com/cta-observatory/ctapipe)

[pyfact](github.com/fact-project/pyfact) for HDF5 writing

All other dependencies are part of the ctapipe environment

## Usage
```simtel.gz``` files should be stored in ```gamma```, ```proton```, or ```gamma-diffuse``` folders in the same directory.

To run simtel_processing within this directory:

```python /PATH/TO/simtel_processing/simtel_processing.py -l LATITUDE,LOGITUDE```

or, from within the simtel_processing folder:

```python simtel_processing.py -i /PATH/TO/SIMTEL_FOLDER -l LATITUDE,LONGITUDE```

where ```LATITUDE,LONGITUDE``` are the latitude and longitude of the telescope site.


There are many options available, viewable with ```python simtel_processing.py --help```
