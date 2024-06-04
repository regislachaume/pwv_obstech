## Brief intro

### Purpose

Precipitable water vapour (PWV) estimation from Global Navigation Satellite Systems (GNSS) measurements using precise point positioning (PPP) plus nearby and/or modelled meteorological data.

### Installation

Within a [virtual environment](https://docs.python.org/3/library/venv.html "venv package"), `pip3 install -e .` to use it locally.  The package is not available from python repositories. 

### Pipeline

Scripts to run the pipeline:
* `ublox_record`: records raw GNSS measurements
* `csrs_pipeline`: scans the disk to submit GNNS measurements and obtain PPP products
* `pwv_pipeline`: scans the disk for PPP products to determine the PWV

### Code structure

* `ublox.py`: records raw GNSS data from a U-Blox receiver
* `csrs.py`: submit to the online PPP service of the Canadian Spatial Reference System and handle results
* `pwv.py`: determines the PWV
* `meteo/`: retrieve and process meteorological data
