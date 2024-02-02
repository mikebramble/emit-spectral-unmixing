# Code written by: Michael Bramble | michael.s.bramble@jpl.nasa.gov
# load unmixing endmembers for EMIT AMD investigation
# 20240124 - initial version

import os
import numpy as np 
import pandas as pd
from scipy.interpolate import CubicSpline


# # # # load usgs wavelength arrays
file = '/Users/bramble/My Drive/_JPL_AMD/endmembers/usgs/splib07a_Wavelengths_BECK_Beckman_0.2-3.0_microns.txt'
df_usgs_beck = pd.read_csv(file,header=None,skiprows=1,names=["wavelength"])

file = '/Users/bramble/My Drive/_JPL_AMD/endmembers/usgs/splib07a_Wavelengths_ASD_0.35-2.5_microns_2151_ch.txt'
df_usgs_asd = pd.read_csv(file,header=None,skiprows=1,names=["wavelength"])


# # # # load spectra

# jarosite
file = '/Users/bramble/My Drive/_JPL_AMD/endmembers/jarosite/c1jb709.txt'
df_c1jb709 = pd.read_csv(file,header=0,skiprows=2,delimiter="\t",usecols=[0, 1],names=["wavelength","reflectance"])

# schwertmannite
file = '/Users/bramble/My Drive/_JPL_AMD/endmembers/usgs/splib07a_Schwertmannite_BZ93-1_BECKb_AREF.txt'
df_schwertmannite_BZ93_1 = pd.read_csv(file,header=None,skiprows=1,names=["reflectance"])

# # # # load EMIT wavelengths
file = '/Users/bramble/My Drive/_JPL_AMD/endmembers/emit_wavelengths.txt'
wavelengths_emit_nm = pd.read_csv(file,header=None)
wavelengths_emit_um = wavelengths_emit[0]/1000

# # # # resample to EMIT wavelenths

# jarosite
spl = CubicSpline(df_c1jb709['wavelength'], df_c1jb709['reflectance'])
y_c1jb709 = spl(wavelengths_emit_um)

# schwertmannite
spl = CubicSpline(df_usgs_beck['wavelength'], df_schwertmannite_BZ93_1['reflectance'])
y_schwertmannite_BZ93_1 = spl(wavelengths_emit_um)


# set bad bands to -0.01 or NaN



# format endmembers into array