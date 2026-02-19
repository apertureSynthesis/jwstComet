This package provides utilities for analyzing comet IFU spectra acquired with the JWST. 
Its primary purpose is to serve as a utility tool for generating extracts with various 
viewing geometries from JWST IFU spectra using tools available in the astropy and photutils 
packages. Its secondary purpose is to automate communication with the NASA Planetary Spectrum 
Generator (PSG) API to facilitate performing retrievals or other analyses of the aforementioned 
JWST comet spectra. 

Note that a local installation of the PSG is highly recommended for the mapping routine, which 
will make hundreds of calls to the API in a short timeframe.

Extraction:
-Beam: extracts a spectrum from an IFU cube at a specified position
-Annulus: extracts an azimuthally averaged spectrum at a specified position

Modeling:
-runPSG: sends an extracted spectrum to the PSG for modeling, including running forward models and retrievals

Mapping:
-Mapping: Extracts spectra at each pixel in the IFU followed by modeling with PSG to create a map of column density or temperature

Requirements:
Anaconda
python >= 3.11
numpy
pandas
matplotlib
astropy
scipy
photutils
jupyter
astroquery

Installation can use conda or pyenv:
Clone the repository and navigate to the jwstComet directory. Included is a .yml file to create a suitable conda environment
or a requirements.txt for pyenv:

For conda:
conda env create --file jwstComet.yml

Activate the new conda environment (named jwstComet), and add the parent directory to your conda path:

conda develop /path-to-jwstComet/

For pyenv, after installing pyenv and pyenv-virtualenv:
pyenv virtualenv 3.11 jwstComet
pyenv activate jwstComet
pip install -r requirements.txt

Usage:
Jupyter notebooks are included to provide examples. Example files contain JWST NIRSpec IFU
observations of comet C/2017 K2 (Woodward et al. 2025, Planetary Science Journal, 6, 139)