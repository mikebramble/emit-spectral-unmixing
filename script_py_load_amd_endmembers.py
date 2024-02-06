# Code written by: Michael Bramble | michael.s.bramble@jpl.nasa.gov
# load unmixing endmembers for EMIT AMD investigation
# 20240124 - initial version
# 20240202 - finalized initial version, exports endmembers in a matrix

import os
import numpy as np 
import pandas as pd
import xarray as xr
from scipy.interpolate import CubicSpline

# # # # 
# # # # load usgs wavelength arrays
# # # # 

file = '/Users/bramble/My Drive/_JPL_AMD/endmembers/usgs/splib07a_Wavelengths_BECK_Beckman_0.2-3.0_microns.txt'
df_usgs_beck = pd.read_csv(file,header=None,skiprows=1,names=["wavelength"])

file = '/Users/bramble/My Drive/_JPL_AMD/endmembers/usgs/splib07a_Wavelengths_ASD_0.35-2.5_microns_2151_ch.txt'
df_usgs_asd = pd.read_csv(file,header=None,skiprows=1,names=["wavelength"])

# # # # 
# # # # load spectra
# # # # 

# pyrite
file = '/Users/bramble/My Drive/_JPL_AMD/endmembers/usgs/splib07a_Pyrite_S29-4_BECKc_AREF.txt'
df_pyrite_S29_4 = pd.read_csv(file,header=None,skiprows=1,names=["reflectance"])

file = '/Users/bramble/My Drive/_JPL_AMD/endmembers/usgs/splib07a_Pyrite_LV95-6A_Weath_on_Tail_BECKb_AREF.txt'
df_pyrite_LV95_6A = pd.read_csv(file,header=None,skiprows=1,names=["reflectance"])

# goethite
file = '/Users/bramble/My Drive/_JPL_AMD/endmembers/usgs/splib07a_Goethite_WS222_Medium_Gr._BECKa_AREF.txt'
df_goethite_WS222_medium_gr = pd.read_csv(file,header=None,skiprows=1,names=["reflectance"])

# jarosite
file = '/Users/bramble/My Drive/_JPL_AMD/endmembers/jarosite/c1jb709.txt'
df_c1jb709 = pd.read_csv(file,header=0,skiprows=2,delimiter="\t",usecols=[0, 1],names=["wavelength","reflectance"])

file = '/Users/bramble/Downloads/usgs_splib07/ASCIIdata/ASCIIdata_splib07a/ChapterM_Minerals/splib07a_Jarosite_GDS99_K_200C_Syn_BECKa_AREF.txt'
df_jarosite_GDS99_K_200C = pd.read_csv(file,header=None,skiprows=1,names=["reflectance"])

# schwertmannite
file = '/Users/bramble/My Drive/_JPL_AMD/endmembers/usgs/splib07a_Schwertmannite_BZ93-1_BECKb_AREF.txt'
df_schwertmannite_BZ93_1 = pd.read_csv(file,header=None,skiprows=1,names=["reflectance"])

# copiapite
file = '/Users/bramble/My Drive/_JPL_AMD/endmembers/usgs/splib07a_Copiapite_GDS21_BECKb_AREF.txt'
df_copiapite_GDS21 = pd.read_csv(file,header=None,skiprows=1,names=["reflectance"])

# hematite
file = '/Users/bramble/My Drive/_JPL_AMD/endmembers/usgs/splib07a_Hematite_GDS27_BECKa_AREF.txt'
df_hematite_GDS27 = pd.read_csv(file,header=None,skiprows=1,names=["reflectance"])

# ferrihydrite
file = '/Users/bramble/My Drive/_JPL_AMD/endmembers/usgs/splib07a_Ferrihydrite_GDS75_Syn_F6_BECKb_AREF.txt'
df_ferrihydrite_GDS75_syn_f6 = pd.read_csv(file,header=None,skiprows=1,names=["reflectance"])

# # # # 
# # # # load EMIT wavelengths
# # # # 

file = '/Users/bramble/My Drive/_JPL_AMD/endmembers/emit_wavelengths.txt'
wavelengths_emit_nm = pd.read_csv(file,header=None)
wavelengths_emit_um = wavelengths_emit_nm[0]/1000

# # # # 
# # # # resample to EMIT wavelenths
# # # # 

# pyrite
spl = CubicSpline(df_usgs_beck['wavelength'], df_pyrite_S29_4['reflectance'].replace(-1.23e+34,0))
y_pyrite_S29_4 = spl(wavelengths_emit_um)

spl = CubicSpline(df_usgs_beck['wavelength'], df_pyrite_LV95_6A['reflectance'].replace(-1.23e+34,0))
y_pyrite_LV95_6A = spl(wavelengths_emit_um)

# goethite
spl = CubicSpline(df_usgs_beck['wavelength'], df_goethite_WS222_medium_gr['reflectance'].replace(-1.23e+34,0))
y_goethite_WS222_medium_gr = spl(wavelengths_emit_um)

# jarosite
spl = CubicSpline(df_c1jb709['wavelength'], df_c1jb709['reflectance'])
y_c1jb709 = spl(wavelengths_emit_um)

spl = CubicSpline(df_usgs_beck['wavelength'], df_jarosite_GDS99_K_200C['reflectance'].replace(-1.23e+34,0))
y_jarosite_GDS99_K_200C = spl(wavelengths_emit_um)

# schwertmannite
spl = CubicSpline(df_usgs_beck['wavelength'], df_schwertmannite_BZ93_1['reflectance'].replace(-1.23e+34,0))
y_schwertmannite_BZ93_1 = spl(wavelengths_emit_um)

# copiapite
spl = CubicSpline(df_usgs_beck['wavelength'], df_copiapite_GDS21['reflectance'].replace(-1.23e+34,0))
y_copiapite_GDS21 = spl(wavelengths_emit_um)

# hematite
spl = CubicSpline(df_usgs_beck['wavelength'], df_hematite_GDS27['reflectance'].replace(-1.23e+34,0))
y_hematite_GDS27 = spl(wavelengths_emit_um)

# ferrihydrite
spl = CubicSpline(df_usgs_beck['wavelength'], df_ferrihydrite_GDS75_syn_f6['reflectance'].replace(-1.23e+34,0))
y_ferrihydrite_GDS75_syn_f6 = spl(wavelengths_emit_um)

# # # # 
# # # # 
# # # # 
# set bad bands to -0.01 or NaN


# # # # 
# format desired endmembers into an array
# # # # 

endmembers = np.column_stack((y_pyrite_LV95_6A, y_pyrite_LV95_6A, y_hematite_GDS27, y_goethite_WS222_medium_gr, y_jarosite_GDS99_K_200C, y_schwertmannite_BZ93_1, y_copiapite_GDS21 ,y_ferrihydrite_GDS75_syn_f6))