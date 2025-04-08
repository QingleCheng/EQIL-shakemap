# Time-history analysis-based prompt assessment method for regional earthquake-induced landslides

## Introduction

- The proposed method uses the ShakeMap method to construct the response spectrum of the ground motion field for the target region and then uses the continuous wavelet transform to generate the ground motion field. The time-history-based Newmark sliding block method is then adopted to calculate the probability of landslide occurrence. 

### Framework

![Framework](./doc/tr_paper.png)



## Setup

> Since this program needs to use the [ShakeMap](https://github.com/usgs/shakemap) program, it needs to be run on a system supported by [ShakeMap](https://github.com/usgs/shakemap)!

1. Install the ShakeMap software by referring to the ShakeMap Installation Guide available at https://github.com/usgs/shakemap/wiki/Installation.
2.Activate the ShakeMap environment.
3.Clone the project repository and compile the required features. Execute the program from the root directory.

```shell
$ git clone https://gitee.com/guyi2000/land-slide landslide
$ cd ./landslide
$ mkdir -v ./utilities/build
$ pushd ./utilities/build
$ cmake ..
$ make && popd
```


4. Once the aforementioned commands are executed, the programs `GetStationData.exe` and `NewmarkDisp.exe` will be generated in the `./utilities` directory.

## Quick Start

 
1. Modify the information in `config.py`, such as the number of cores used for computation, the path to data files, etc. For example:
```python
NUM_THREADS = 8       # number of cores in multi-threading
NUM_SIDE_STATION = 2  # Number of sampling stations per side
STATION_STEP = 10     # Resolution of landslide calculation, where 1 is the highest
```



2. Activate the  `ShakeMap` environment.
3. Prepare the [required data] (see Data Description).
4. Run `python landslide.py` in the command line.

- Example


```shell
$ vim ./config.py
$ conda activate shakemap
(shakemap) $ python landslide.py
```

### Data Description

- Ground motion data of stations
  - Stored in the form of folders, with the folder name being the station name, which cannot have illegal characters. The folder contains three files: `EW.txt`,`NS.txt`,`UD.txt`, which are in the same format as general seismic data.
- Station location data
  - Stored in the form of a text file, containing key-value pairs of the station name and its latitude and longitude.
- Geological data
  - Stored in the form of a text file containing regional geological data.
  - The format is `ID terrain_type longitude latitude`, with spaces as separators.
- Slope data
  - tored in `GeoTiff`  format containing regional slope data.
- Shear wave velocity at 30m underground
  - Download from [here](https://earthquake.usgs.gov/data/vs30/) or use the data provided by ShakeMap.
  - Stored in grd format as a raster dataset.
  
- Example

```shell
$ mkdir -v {tiff_data,lithologic_data,vs30_data,sample_data}
$ cp ../shakemap_data/vs30/global_vs30.grd ./vs30_data/
$ cp -r /datadir/station ./sample_data/
$ cp /datadir/slope.tif ./tiff_data/
$ cp /datadir/litho.txt ./lithologic_data/
```



### To generate an API documentation

1. Activate the `ShakeMap`  environment.
2. Install [Doxygen](https://www.doxygen.nl/index.html)
3. Install [Sphinx](https://www.sphinx-doc.org/en/master/index.html)
4. Install [Sphinx RTD Theme](https://readthedocs.org/)
5. Run `make html` in the `./doc` directory to generate Python and C++ documentation.

- `Ubuntu` example code

```shell
$ sudo apt install doxygen doxygen-doc -y
$ conda activate shakemap
(shakemap) $ pip install -U Sphinx
(shakemap) $ pip install -U sphinx_rtd_theme
(shakemap) $ cd doc
(shakemap) $ make html
```


### How to contribute the project?

- Just fork the project and raise a pull request
