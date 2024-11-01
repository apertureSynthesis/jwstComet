This package provides utilities for analyzing comet IFU spectra acquired with the JWST.

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
Jupyter notebooks are included to provide examples