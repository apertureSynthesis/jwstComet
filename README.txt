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

Installation:
Clone the repository and navigate to the jwstComet directory. Included is a .yml file to create a suitable conda environment:

conda env create --file jwstComet.yml

Activate the new conda environment, and add the parent directory to your conda path:

conda develop /path-to-jwstComet/

Usage:
Jupyter notebooks are included to provide examples