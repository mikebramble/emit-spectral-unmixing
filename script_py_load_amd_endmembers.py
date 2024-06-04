# Code written by: Michael Bramble | michael.s.bramble@jpl.nasa.gov
# load unmixing endmembers for EMIT AMD investigation. ASCII spectra files are loaded and then resampled to
# the EMIT wavelenths. 
# 20240124 - initial version
# 20240202 - finalized initial version, exports endmembers in a matrix
# 20240222 - completed first set of endmembers

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

file = '/Users/bramble/My Drive/_JPL_AMD/endmembers/usgs/splib07a_Pyrite_S26-8_BECKc_AREF.txt'
df_pyrite_S26_8 = pd.read_csv(file,header=None,skiprows=1,names=["reflectance"])

file = '/Users/bramble/My Drive/_JPL_AMD/endmembers/usgs/splib07a_Pyrite_S30_BECKc_AREF.txt'
df_pyrite_S30 = pd.read_csv(file,header=None,skiprows=1,names=["reflectance"])

file = '/Users/bramble/My Drive/_JPL_AMD/endmembers/usgs/splib07a_Pyrite_HS35.3_BECKb_AREF.txt'
df_pyrite_HS35_3 = pd.read_csv(file,header=None,skiprows=1,names=["reflectance"])

file = '/Users/bramble/My Drive/_JPL_AMD/endmembers/usgs/splib07a_Pyrite_GDS483.c_30-60um_ASDFRc_AREF.txt'
df_pyrite_GDS483 = pd.read_csv(file,header=None,skiprows=1,names=["reflectance"])

file = '/Users/bramble/My Drive/_JPL_AMD/endmembers/pyrite/c1sc77.txt'
df_c1sc77 = pd.read_csv(file,header=0,skiprows=2,delimiter="\t",usecols=[0, 1],names=["wavelength","reflectance"])

file = '/Users/bramble/My Drive/_JPL_AMD/endmembers/pyrite/c1sh41.txt'
df_c1sh41 = pd.read_csv(file,header=0,skiprows=2,delimiter="\t",usecols=[0, 1],names=["wavelength","reflectance"])

file = '/Users/bramble/My Drive/_JPL_AMD/endmembers/pyrite/cash52.txt'
df_cash52 = pd.read_csv(file,header=0,skiprows=2,delimiter="\t",usecols=[0, 1],names=["wavelength","reflectance"])

file = '/Users/bramble/My Drive/_JPL_AMD/endmembers/pyrite/cbsh52.txt'
df_cbsh52 = pd.read_csv(file,header=0,skiprows=2,delimiter="\t",usecols=[0, 1],names=["wavelength","reflectance"])

file = '/Users/bramble/My Drive/_JPL_AMD/endmembers/pyrite/ccsh52.txt'
df_ccsh52 = pd.read_csv(file,header=0,skiprows=2,delimiter="\t",usecols=[0, 1],names=["wavelength","reflectance"])

# goethite
file = '/Users/bramble/My Drive/_JPL_AMD/endmembers/usgs/splib07a_Goethite_WS222_Medium_Gr._BECKa_AREF.txt'
df_goethite_WS222_medium_gr = pd.read_csv(file,header=None,skiprows=1,names=["reflectance"])

file = '/Users/bramble/My Drive/_JPL_AMD/endmembers/usgs/splib07a_Goethite_WS222_Coarse_Gr._BECKa_AREF.txt'
df_goethite_WS222_coarse_gr = pd.read_csv(file,header=None,skiprows=1,names=["reflectance"])

file = '/Users/bramble/My Drive/_JPL_AMD/endmembers/usgs/splib07a_Goethite_Thin_Film_WS222_BECKa_AREF.txt'
df_goethite_WS222_thin_film = pd.read_csv(file,header=None,skiprows=1,names=["reflectance"])

file = '/Users/bramble/My Drive/_JPL_AMD/endmembers/usgs/splib07a_Goethite_GDS134_ASDFRb_AREF.txt'
df_goethite_GDS134 = pd.read_csv(file,header=None,skiprows=1,names=["reflectance"])

file = '/Users/bramble/My Drive/_JPL_AMD/endmembers/goethite/c1go01.txt'
df_c1go01 = pd.read_csv(file,header=0,skiprows=2,delimiter="\t",usecols=[0, 1],names=["wavelength","reflectance"])

file = '/Users/bramble/My Drive/_JPL_AMD/endmembers/goethite/c1jb54.txt'
df_c1jb54 = pd.read_csv(file,header=0,skiprows=2,delimiter="\t",usecols=[0, 1],names=["wavelength","reflectance"])

file = '/Users/bramble/My Drive/_JPL_AMD/endmembers/goethite/c1jbh58a.txt'
df_c1jbh58a = pd.read_csv(file,header=0,skiprows=2,delimiter="\t",usecols=[0, 1],names=["wavelength","reflectance"])

file = '/Users/bramble/My Drive/_JPL_AMD/endmembers/goethite/caho03.txt'
df_caho03 = pd.read_csv(file,header=0,skiprows=2,delimiter="\t",usecols=[0, 1],names=["wavelength","reflectance"])

# jarosite
file = '/Users/bramble/My Drive/_JPL_AMD/endmembers/jarosite/c1jb709.txt'
df_c1jb709 = pd.read_csv(file,header=0,skiprows=2,delimiter="\t",usecols=[0, 1],names=["wavelength","reflectance"])

file = '/Users/bramble/My Drive/_JPL_AMD/endmembers/usgs/splib07a_Jarosite_GDS99_K_200C_Syn_BECKa_AREF.txt'
df_jarosite_GDS99_K_200C = pd.read_csv(file,header=None,skiprows=1,names=["reflectance"])

file = '/Users/bramble/My Drive/_JPL_AMD/endmembers/usgs/splib07a_Jarosite_GDS98_K_90C_Syn_BECKa_AREF.txt'
df_jarosite_GDS98_K_90C = pd.read_csv(file,header=None,skiprows=1,names=["reflectance"])

file = '/Users/bramble/My Drive/_JPL_AMD/endmembers/usgs/splib07a_Jarosite_GDS636_K_Penalt325um_ASDNGa_AREF.txt'
df_jarosite_GDS636 = pd.read_csv(file,header=None,skiprows=1,names=["reflectance"])

file = '/Users/bramble/My Drive/_JPL_AMD/endmembers/usgs/splib07a_Jarosite_GDS732_K_200CSyn6hr_ASDNGb_AREF.txt'
df_jarosite_GDS732 = pd.read_csv(file,header=None,skiprows=1,names=["reflectance"])

file = '/Users/bramble/My Drive/_JPL_AMD/endmembers/usgs/splib07a_Jarosite_JR2501_(K)_BECKb_AREF.txt'
df_jarosite_JR2501 = pd.read_csv(file,header=None,skiprows=1,names=["reflectance"])

file = '/Users/bramble/My Drive/_JPL_AMD/endmembers/usgs/splib07a_Jarosite_SJ-1_H3O_10-20%_BECKb_AREF.txt'
df_jarosite_SJ_1 = pd.read_csv(file,header=None,skiprows=1,names=["reflectance"])

file = '/Users/bramble/My Drive/_JPL_AMD/endmembers/usgs/splib07a_Jarosite_Thin_Film_GDS243_BECKb_AREF.txt'
df_jarosite_GDS243 = pd.read_csv(file,header=None,skiprows=1,names=["reflectance"])

file = '/Users/bramble/My Drive/_JPL_AMD/endmembers/usgs/splib07a_Jarosite_GDS24_Na_BECKb_AREF.txt'
df_jarosite_GDS24 = pd.read_csv(file,header=None,skiprows=1,names=["reflectance"])

file = '/Users/bramble/My Drive/_JPL_AMD/endmembers/jarosite/c1jb709.txt'
df_c1jb709 = pd.read_csv(file,header=0,skiprows=2,delimiter="\t",usecols=[0, 1],names=["wavelength","reflectance"])

file = '/Users/bramble/My Drive/_JPL_AMD/endmembers/jarosite/c1jb991.txt'
df_c1jb991 = pd.read_csv(file,header=0,skiprows=2,delimiter="\t",usecols=[0, 1],names=["wavelength","reflectance"])

file = '/Users/bramble/My Drive/_JPL_AMD/endmembers/jarosite/c1jba78.txt'
df_c1jba78 = pd.read_csv(file,header=0,skiprows=2,delimiter="\t",usecols=[0, 1],names=["wavelength","reflectance"])

file = '/Users/bramble/My Drive/_JPL_AMD/endmembers/jarosite/c1jbb79.txt'
df_c1jbb79 = pd.read_csv(file,header=0,skiprows=2,delimiter="\t",usecols=[0, 1],names=["wavelength","reflectance"])

file = '/Users/bramble/My Drive/_JPL_AMD/endmembers/jarosite/c1pj02.txt'
df_c1pj02 = pd.read_csv(file,header=0,skiprows=2,delimiter="\t",usecols=[0, 1],names=["wavelength","reflectance"])

file = '/Users/bramble/My Drive/_JPL_AMD/endmembers/jarosite/c1pj03.txt'
df_c1pj03 = pd.read_csv(file,header=0,skiprows=2,delimiter="\t",usecols=[0, 1],names=["wavelength","reflectance"])

# schwertmannite
file = '/Users/bramble/My Drive/_JPL_AMD/endmembers/usgs/splib07a_Schwertmannite_BZ93-1_BECKb_AREF.txt'
df_schwertmannite_BZ93_1 = pd.read_csv(file,header=None,skiprows=1,names=["reflectance"])

file = '/Users/bramble/My Drive/_JPL_AMD/endmembers/schwertmannite/c1jb130a.txt'
df_c1jb130a = pd.read_csv(file,header=0,skiprows=2,delimiter="\t",usecols=[0, 1],names=["wavelength","reflectance"])

file = '/Users/bramble/My Drive/_JPL_AMD/endmembers/schwertmannite/c1jb130b.txt'
df_c1jb130b = pd.read_csv(file,header=0,skiprows=2,delimiter="\t",usecols=[0, 1],names=["wavelength","reflectance"])

file = '/Users/bramble/My Drive/_JPL_AMD/endmembers/schwertmannite/c1jb131.txt'
df_c1jb131 = pd.read_csv(file,header=0,skiprows=2,delimiter="\t",usecols=[0, 1],names=["wavelength","reflectance"])

# copiapite
file = '/Users/bramble/My Drive/_JPL_AMD/endmembers/usgs/splib07a_Copiapite_GDS21_BECKb_AREF.txt'
df_copiapite_GDS21 = pd.read_csv(file,header=None,skiprows=1,names=["reflectance"])

file = '/Users/bramble/My Drive/_JPL_AMD/endmembers/copiapite/c1jb620a.txt'
df_c1jb620a = pd.read_csv(file,header=0,skiprows=2,delimiter="\t",usecols=[0, 1],names=["wavelength","reflectance"])

file = '/Users/bramble/My Drive/_JPL_AMD/endmembers/copiapite/c1jb834.txt'
df_c1jb834 = pd.read_csv(file,header=0,skiprows=2,delimiter="\t",usecols=[0, 1],names=["wavelength","reflectance"])

file = '/Users/bramble/My Drive/_JPL_AMD/endmembers/copiapite/c1jb835.txt'
df_c1jb835 = pd.read_csv(file,header=0,skiprows=2,delimiter="\t",usecols=[0, 1],names=["wavelength","reflectance"])

file = '/Users/bramble/My Drive/_JPL_AMD/endmembers/copiapite/c1jb969.txt'
df_c1jb969 = pd.read_csv(file,header=0,skiprows=2,delimiter="\t",usecols=[0, 1],names=["wavelength","reflectance"])

file = '/Users/bramble/My Drive/_JPL_AMD/endmembers/copiapite/c1jb970.txt'
df_c1jb970 = pd.read_csv(file,header=0,skiprows=2,delimiter="\t",usecols=[0, 1],names=["wavelength","reflectance"])

file = '/Users/bramble/My Drive/_JPL_AMD/endmembers/copiapite/c1jb971.txt'
df_c1jb971 = pd.read_csv(file,header=0,skiprows=2,delimiter="\t",usecols=[0, 1],names=["wavelength","reflectance"])

file = '/Users/bramble/My Drive/_JPL_AMD/endmembers/copiapite/c1jb972.txt'
df_c1jb972 = pd.read_csv(file,header=0,skiprows=2,delimiter="\t",usecols=[0, 1],names=["wavelength","reflectance"])

file = '/Users/bramble/My Drive/_JPL_AMD/endmembers/copiapite/c1jb974.txt'
df_c1jb974 = pd.read_csv(file,header=0,skiprows=2,delimiter="\t",usecols=[0, 1],names=["wavelength","reflectance"])

file = '/Users/bramble/My Drive/_JPL_AMD/endmembers/copiapite/c1jba57.txt'
df_c1jba57 = pd.read_csv(file,header=0,skiprows=2,delimiter="\t",usecols=[0, 1],names=["wavelength","reflectance"])

file = '/Users/bramble/My Drive/_JPL_AMD/endmembers/copiapite/c1jba79.txt'
df_c1jba79 = pd.read_csv(file,header=0,skiprows=2,delimiter="\t",usecols=[0, 1],names=["wavelength","reflectance"])

file = '/Users/bramble/My Drive/_JPL_AMD/endmembers/copiapite/c1jbb73.txt'
df_c1jbb73 = pd.read_csv(file,header=0,skiprows=2,delimiter="\t",usecols=[0, 1],names=["wavelength","reflectance"])

file = '/Users/bramble/My Drive/_JPL_AMD/endmembers/copiapite/c1pc30.txt'
df_c1pc30 = pd.read_csv(file,header=0,skiprows=2,delimiter="\t",usecols=[0, 1],names=["wavelength","reflectance"])

file = '/Users/bramble/My Drive/_JPL_AMD/endmembers/copiapite/c1rm35.txt'
df_c1rm35 = pd.read_csv(file,header=0,skiprows=2,delimiter="\t",usecols=[0, 1],names=["wavelength","reflectance"])

file = '/Users/bramble/My Drive/_JPL_AMD/endmembers/copiapite/cacc13.txt'
df_cacc13 = pd.read_csv(file,header=0,skiprows=2,delimiter="\t",usecols=[0, 1],names=["wavelength","reflectance"])

file = '/Users/bramble/My Drive/_JPL_AMD/endmembers/copiapite/casf39.txt'
df_casf39 = pd.read_csv(file,header=0,skiprows=2,delimiter="\t",usecols=[0, 1],names=["wavelength","reflectance"])

# hematite
file = '/Users/bramble/My Drive/_JPL_AMD/endmembers/usgs/splib07a_Hematite_GDS27_BECKa_AREF.txt'
df_hematite_GDS27 = pd.read_csv(file,header=None,skiprows=1,names=["reflectance"])

file = '/Users/bramble/My Drive/_JPL_AMD/endmembers/usgs/splib07a_Hematite_Thin_Film_GDS27_BECKa_AREF.txt'
df_hematite_GDS27_thin_film = pd.read_csv(file,header=None,skiprows=1,names=["reflectance"])

file = '/Users/bramble/My Drive/_JPL_AMD/endmembers/usgs/splib07a_Hematite_FE2602_BECKb_AREF.txt'
df_hematite_FE2602 = pd.read_csv(file,header=None,skiprows=1,names=["reflectance"])

file = '/Users/bramble/My Drive/_JPL_AMD/endmembers/usgs/splib07a_Hematite_GDS69.b_104-150u_BECKb_AREF.txt'
df_hematite_GDS69b = pd.read_csv(file,header=None,skiprows=1,names=["reflectance"])

file = '/Users/bramble/My Drive/_JPL_AMD/endmembers/usgs/splib07a_Hematite_GDS69.c_60-104um_BECKb_AREF.txt'
df_hematite_GDS69c = pd.read_csv(file,header=None,skiprows=1,names=["reflectance"])

file = '/Users/bramble/My Drive/_JPL_AMD/endmembers/usgs/splib07a_Hematite_GDS69.d_30-45um_BECKb_AREF.txt'
df_hematite_GDS69d = pd.read_csv(file,header=None,skiprows=1,names=["reflectance"])

file = '/Users/bramble/My Drive/_JPL_AMD/endmembers/usgs/splib07a_Hematite_GDS69.e_20-30um_BECKb_AREF.txt'
df_hematite_GDS69e = pd.read_csv(file,header=None,skiprows=1,names=["reflectance"])

file = '/Users/bramble/My Drive/_JPL_AMD/endmembers/usgs/splib07a_Hematite_GDS69.f_10-20um_BECKb_AREF.txt'
df_hematite_GDS69f = pd.read_csv(file,header=None,skiprows=1,names=["reflectance"])

file = '/Users/bramble/My Drive/_JPL_AMD/endmembers/usgs/splib07a_Hematite_GDS69.g_lt10um_BECKb_AREF.txt'
df_hematite_GDS69g = pd.read_csv(file,header=None,skiprows=1,names=["reflectance"])

file = '/Users/bramble/My Drive/_JPL_AMD/endmembers/hematite/c1cy11.txt'
df_c1cy11 = pd.read_csv(file,header=0,skiprows=2,delimiter="\t",usecols=[0, 1],names=["wavelength","reflectance"])

file = '/Users/bramble/My Drive/_JPL_AMD/endmembers/hematite/c1jb256.txt'
df_c1jb256 = pd.read_csv(file,header=0,skiprows=2,delimiter="\t",usecols=[0, 1],names=["wavelength","reflectance"])

file = '/Users/bramble/My Drive/_JPL_AMD/endmembers/hematite/c1jb257.txt'
df_c1jb257 = pd.read_csv(file,header=0,skiprows=2,delimiter="\t",usecols=[0, 1],names=["wavelength","reflectance"])

file = '/Users/bramble/My Drive/_JPL_AMD/endmembers/hematite/c1jb993.txt'
df_c1jb993 = pd.read_csv(file,header=0,skiprows=2,delimiter="\t",usecols=[0, 1],names=["wavelength","reflectance"])

file = '/Users/bramble/My Drive/_JPL_AMD/endmembers/hematite/cafe02.txt'
df_cafe02 = pd.read_csv(file,header=0,skiprows=2,delimiter="\t",usecols=[0, 1],names=["wavelength","reflectance"])

# ferrihydrite
file = '/Users/bramble/My Drive/_JPL_AMD/endmembers/usgs/splib07a_Ferrihydrite_GDS75_Syn_F6_BECKb_AREF.txt'
df_ferrihydrite_GDS75_syn_f6 = pd.read_csv(file,header=None,skiprows=1,names=["reflectance"])

file = '/Users/bramble/My Drive/_JPL_AMD/endmembers/ferrihydrite/c1jb45.txt'
df_c1jb45 = pd.read_csv(file,header=0,skiprows=2,delimiter="\t",usecols=[0, 1],names=["wavelength","reflectance"])

file = '/Users/bramble/My Drive/_JPL_AMD/endmembers/ferrihydrite/c1jb46.txt'
df_c1jb46 = pd.read_csv(file,header=0,skiprows=2,delimiter="\t",usecols=[0, 1],names=["wavelength","reflectance"])

file = '/Users/bramble/My Drive/_JPL_AMD/endmembers/ferrihydrite/c1jb564a.txt'
df_c1jb564a = pd.read_csv(file,header=0,skiprows=2,delimiter="\t",usecols=[0, 1],names=["wavelength","reflectance"])

file = '/Users/bramble/My Drive/_JPL_AMD/endmembers/ferrihydrite/c1jb565a.txt'
df_c1jb565a = pd.read_csv(file,header=0,skiprows=2,delimiter="\t",usecols=[0, 1],names=["wavelength","reflectance"])

file = '/Users/bramble/My Drive/_JPL_AMD/endmembers/ferrihydrite/c1jm67.txt'
df_c1jm67 = pd.read_csv(file,header=0,skiprows=2,delimiter="\t",usecols=[0, 1],names=["wavelength","reflectance"])

file = '/Users/bramble/My Drive/_JPL_AMD/endmembers/ferrihydrite/c2jb255.txt'
df_c2jb255 = pd.read_csv(file,header=0,skiprows=2,delimiter="\t",usecols=[0, 1],names=["wavelength","reflectance"])

# melanterite
file = '/Users/bramble/My Drive/_JPL_AMD/endmembers/melanterite/c1lh41.txt'
df_c1lh41 = pd.read_csv(file,header=0,skiprows=2,delimiter="\t",usecols=[0, 1],names=["wavelength","reflectance"])

file = '/Users/bramble/My Drive/_JPL_AMD/endmembers/melanterite/casf44.txt'
df_casf44 = pd.read_csv(file,header=0,skiprows=2,delimiter="\t",usecols=[0, 1],names=["wavelength","reflectance"])

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

spl = CubicSpline(df_usgs_beck['wavelength'], df_pyrite_S26_8['reflectance'].replace(-1.23e+34,0))
y_pyrite_S26_8 = spl(wavelengths_emit_um)

spl = CubicSpline(df_usgs_beck['wavelength'], df_pyrite_S30['reflectance'].replace(-1.23e+34,0))
y_pyrite_S30 = spl(wavelengths_emit_um)

spl = CubicSpline(df_usgs_beck['wavelength'], df_pyrite_HS35_3['reflectance'].replace(-1.23e+34,0))
y_pyrite_HS35_3 = spl(wavelengths_emit_um)

spl = CubicSpline(df_usgs_asd['wavelength'], df_pyrite_GDS483['reflectance'].replace(-1.23e+34,0))
y_pyrite_GDS483 = spl(wavelengths_emit_um)

spl = CubicSpline(df_c1sc77['wavelength'], df_c1sc77['reflectance'])
y_pyrite_c1sc77 = spl(wavelengths_emit_um)

spl = CubicSpline(df_c1sh41['wavelength'], df_c1sh41['reflectance'])
y_pyrite_c1sh41 = spl(wavelengths_emit_um)

spl = CubicSpline(df_cash52['wavelength'], df_cash52['reflectance'])
y_pyrite_cash52 = spl(wavelengths_emit_um)

spl = CubicSpline(df_cbsh52['wavelength'], df_cbsh52['reflectance'])
y_pyrite_cbsh52 = spl(wavelengths_emit_um)

spl = CubicSpline(df_ccsh52['wavelength'], df_ccsh52['reflectance'])
y_pyrite_ccsh52 = spl(wavelengths_emit_um)

# goethite
spl = CubicSpline(df_usgs_beck['wavelength'], df_goethite_WS222_medium_gr['reflectance'].replace(-1.23e+34,0))
y_goethite_WS222_medium_gr = spl(wavelengths_emit_um)

spl = CubicSpline(df_usgs_beck['wavelength'], df_goethite_WS222_coarse_gr['reflectance'].replace(-1.23e+34,0))
y_goethite_WS222_coarse_gr = spl(wavelengths_emit_um)

spl = CubicSpline(df_usgs_beck['wavelength'], df_goethite_WS222_thin_film['reflectance'].replace(-1.23e+34,0))
y_goethite_WS222_thin_film = spl(wavelengths_emit_um)

spl = CubicSpline(df_usgs_asd['wavelength'], df_goethite_GDS134['reflectance'].replace(-1.23e+34,0))
y_goethite_GDS134 = spl(wavelengths_emit_um)

spl = CubicSpline(df_c1go01['wavelength'], df_c1go01['reflectance'])
y_goethite_c1go01 = spl(wavelengths_emit_um)

spl = CubicSpline(df_c1jb54['wavelength'], df_c1jb54['reflectance'])
y_goethite_c1jb54 = spl(wavelengths_emit_um)

spl = CubicSpline(df_c1jbh58a['wavelength'], df_c1jbh58a['reflectance'])
y_goethite_c1jbh58a = spl(wavelengths_emit_um)

spl = CubicSpline(df_caho03['wavelength'], df_caho03['reflectance'])
y_goethite_caho03 = spl(wavelengths_emit_um)

# jarosite
spl = CubicSpline(df_c1jb709['wavelength'], df_c1jb709['reflectance'])
y_jarosite_c1jb709 = spl(wavelengths_emit_um)

spl = CubicSpline(df_usgs_beck['wavelength'], df_jarosite_GDS99_K_200C['reflectance'].replace(-1.23e+34,0))
y_jarosite_GDS99_K_200C = spl(wavelengths_emit_um)

spl = CubicSpline(df_usgs_beck['wavelength'], df_jarosite_GDS98_K_90C['reflectance'].replace(-1.23e+34,0))
y_jarosite_GDS98_K_90C = spl(wavelengths_emit_um)

spl = CubicSpline(df_usgs_asd['wavelength'], df_jarosite_GDS732['reflectance'].replace(-1.23e+34,0))
y_jarosite_GDS732 = spl(wavelengths_emit_um)

spl = CubicSpline(df_usgs_asd['wavelength'], df_jarosite_GDS636['reflectance'].replace(-1.23e+34,0))
y_jarosite_GDS636 = spl(wavelengths_emit_um)

spl = CubicSpline(df_usgs_beck['wavelength'], df_jarosite_JR2501['reflectance'].replace(-1.23e+34,0))
y_jarosite_JR2501 = spl(wavelengths_emit_um)

spl = CubicSpline(df_usgs_beck['wavelength'], df_jarosite_SJ_1['reflectance'].replace(-1.23e+34,0))
y_jarosite_SJ_1 = spl(wavelengths_emit_um)

spl = CubicSpline(df_usgs_beck['wavelength'], df_jarosite_GDS243['reflectance'].replace(-1.23e+34,0))
y_jarosite_GDS243 = spl(wavelengths_emit_um)

spl = CubicSpline(df_usgs_beck['wavelength'], df_jarosite_GDS24['reflectance'].replace(-1.23e+34,0))
y_jarosite_GDS24 = spl(wavelengths_emit_um)

spl = CubicSpline(df_c1jb991['wavelength'], df_c1jb991['reflectance'])
y_jarosite_c1jb991 = spl(wavelengths_emit_um)

spl = CubicSpline(df_c1jba78['wavelength'], df_c1jba78['reflectance'])
y_jarosite_c1jba78 = spl(wavelengths_emit_um)

spl = CubicSpline(df_c1jbb79['wavelength'], df_c1jbb79['reflectance'])
y_jarosite_c1jbb79 = spl(wavelengths_emit_um)

spl = CubicSpline(df_c1pj02['wavelength'], df_c1pj02['reflectance'])
y_jarosite_c1pj02 = spl(wavelengths_emit_um)

spl = CubicSpline(df_c1pj03['wavelength'], df_c1pj03['reflectance'])
y_jarosite_c1pj03 = spl(wavelengths_emit_um)

# schwertmannite
spl = CubicSpline(df_usgs_beck['wavelength'], df_schwertmannite_BZ93_1['reflectance'].replace(-1.23e+34,0))
y_schwertmannite_BZ93_1 = spl(wavelengths_emit_um)

spl = CubicSpline(df_c1jb130a['wavelength'], df_c1jb130a['reflectance'])
y_schwertmannite_c1jb130a = spl(wavelengths_emit_um)

spl = CubicSpline(df_c1jb130b['wavelength'], df_c1jb130b['reflectance'])
y_schwertmannite_c1jb130b = spl(wavelengths_emit_um)

spl = CubicSpline(df_c1jb131['wavelength'], df_c1jb131['reflectance'])
y_schwertmannite_c1jb131 = spl(wavelengths_emit_um)

# copiapite
spl = CubicSpline(df_usgs_beck['wavelength'], df_copiapite_GDS21['reflectance'].replace(-1.23e+34,0))
y_copiapite_GDS21 = spl(wavelengths_emit_um)

spl = CubicSpline(df_cacc13['wavelength'], df_cacc13['reflectance'])
y_copiapite_cacc13 = spl(wavelengths_emit_um)

spl = CubicSpline(df_c1jb969['wavelength'], df_c1jb969['reflectance'])
y_copiapite_c1jb969 = spl(wavelengths_emit_um)

spl = CubicSpline(df_c1jb970['wavelength'], df_c1jb970['reflectance'])
y_copiapite_c1jb970 = spl(wavelengths_emit_um)

spl = CubicSpline(df_c1jb971['wavelength'], df_c1jb971['reflectance'])
y_copiapite_c1jb971 = spl(wavelengths_emit_um)

spl = CubicSpline(df_c1jb972['wavelength'], df_c1jb972['reflectance'])
y_copiapite_c1jb972 = spl(wavelengths_emit_um)

spl = CubicSpline(df_c1jb974['wavelength'], df_c1jb974['reflectance'])
y_copiapite_c1jb974 = spl(wavelengths_emit_um)

spl = CubicSpline(df_c1pc30['wavelength'], df_c1pc30['reflectance'])
y_copiapite_c1pc30 = spl(wavelengths_emit_um)

spl = CubicSpline(df_casf39['wavelength'], df_casf39['reflectance'])
y_copiapite_casf39 = spl(wavelengths_emit_um)

spl = CubicSpline(df_c1jba57['wavelength'], df_c1jba57['reflectance'])
y_copiapite_c1jba57 = spl(wavelengths_emit_um)

spl = CubicSpline(df_c1jba79['wavelength'], df_c1jba79['reflectance'])
y_copiapite_mg_c1jba79 = spl(wavelengths_emit_um)

spl = CubicSpline(df_c1jbb73['wavelength'], df_c1jbb73['reflectance'])
y_copiapite_mg_c1jbb73 = spl(wavelengths_emit_um)

spl = CubicSpline(df_c1jb620a['wavelength'], df_c1jb620a['reflectance'])
y_copiapite_fe_c1jb620a = spl(wavelengths_emit_um)

spl = CubicSpline(df_c1rm35['wavelength'], df_c1rm35['reflectance'])
y_copiapite_fe_c1rm35 = spl(wavelengths_emit_um)

spl = CubicSpline(df_c1jb834['wavelength'], df_c1jb834['reflectance'])
y_copiapite_sy_c1jb834 = spl(wavelengths_emit_um)

spl = CubicSpline(df_c1jb835['wavelength'], df_c1jb835['reflectance'])
y_copiapite_sy_c1jb835 = spl(wavelengths_emit_um)

# hematite
spl = CubicSpline(df_usgs_beck['wavelength'], df_hematite_GDS27['reflectance'].replace(-1.23e+34,0))
y_hematite_GDS27 = spl(wavelengths_emit_um)

spl = CubicSpline(df_usgs_beck['wavelength'], df_hematite_GDS27_thin_film['reflectance'].replace(-1.23e+34,0))
y_hematite_GDS27_thin_film = spl(wavelengths_emit_um)

spl = CubicSpline(df_usgs_beck['wavelength'], df_hematite_FE2602['reflectance'].replace(-1.23e+34,0))
y_hematite_FE2602 = spl(wavelengths_emit_um)

spl = CubicSpline(df_usgs_beck['wavelength'], df_hematite_GDS69b['reflectance'].replace(-1.23e+34,0))
y_hematite_GDS69b = spl(wavelengths_emit_um)

spl = CubicSpline(df_usgs_beck['wavelength'], df_hematite_GDS69c['reflectance'].replace(-1.23e+34,0))
y_hematite_GDS69c = spl(wavelengths_emit_um)

spl = CubicSpline(df_usgs_beck['wavelength'], df_hematite_GDS69d['reflectance'].replace(-1.23e+34,0))
y_hematite_GDS69d = spl(wavelengths_emit_um)

spl = CubicSpline(df_usgs_beck['wavelength'], df_hematite_GDS69e['reflectance'].replace(-1.23e+34,0))
y_hematite_GDS69e = spl(wavelengths_emit_um)

spl = CubicSpline(df_usgs_beck['wavelength'], df_hematite_GDS69f['reflectance'].replace(-1.23e+34,0))
y_hematite_GDS69f = spl(wavelengths_emit_um)

spl = CubicSpline(df_usgs_beck['wavelength'], df_hematite_GDS69g['reflectance'].replace(-1.23e+34,0))
y_hematite_GDS69g = spl(wavelengths_emit_um)

spl = CubicSpline(df_c1jb256['wavelength'], df_c1jb256['reflectance'])
y_hematite_c1jb256 = spl(wavelengths_emit_um)

spl = CubicSpline(df_c1jb257['wavelength'], df_c1jb257['reflectance'])
y_hematite_c1jb257 = spl(wavelengths_emit_um)

spl = CubicSpline(df_c1jb993['wavelength'], df_c1jb993['reflectance'])
y_hematite_c1jb993 = spl(wavelengths_emit_um)

spl = CubicSpline(df_cafe02['wavelength'], df_cafe02['reflectance'])
y_hematite_cafe02 = spl(wavelengths_emit_um)

spl = CubicSpline(df_c1cy11['wavelength'], df_c1cy11['reflectance'])
y_hematite_c1cy11 = spl(wavelengths_emit_um)

# ferrihydrite
spl = CubicSpline(df_usgs_beck['wavelength'], df_ferrihydrite_GDS75_syn_f6['reflectance'].replace(-1.23e+34,0))
y_ferrihydrite_GDS75_syn_f6 = spl(wavelengths_emit_um)

spl = CubicSpline(df_c1jb45['wavelength'], df_c1jb45['reflectance'])
y_ferrihydrite_c1jb45 = spl(wavelengths_emit_um)

spl = CubicSpline(df_c1jb46['wavelength'], df_c1jb46['reflectance'])
y_ferrihydrite_c1jb46 = spl(wavelengths_emit_um)

spl = CubicSpline(df_c1jb564a['wavelength'], df_c1jb564a['reflectance'])
y_ferrihydrite_c1jb564a = spl(wavelengths_emit_um)

spl = CubicSpline(df_c1jb565a['wavelength'], df_c1jb565a['reflectance'])
y_ferrihydrite_c1jb565a = spl(wavelengths_emit_um)

spl = CubicSpline(df_c1jm67['wavelength'], df_c1jm67['reflectance'])
y_ferrihydrite_c1jm67 = spl(wavelengths_emit_um)

spl = CubicSpline(df_c2jb255['wavelength'], df_c2jb255['reflectance'])
y_ferrihydrite_c2jb255 = spl(wavelengths_emit_um)

# melanterite
spl = CubicSpline(df_c1lh41['wavelength'], df_c1lh41['reflectance'])
y_melanterite_c1lh41 = spl(wavelengths_emit_um)

spl = CubicSpline(df_casf44['wavelength'], df_casf44['reflectance'])
y_melanterite_casf44 = spl(wavelengths_emit_um)

# # # # 
# # # # 
# # # # 
# set bad bands to -0.01 or NaN


# # # # 
# format desired endmembers into an array
# # # # 

endmembers = np.column_stack((y_pyrite_S29_4, y_pyrite_LV95_6A, y_hematite_GDS27, y_goethite_WS222_medium_gr, y_jarosite_GDS99_K_200C, y_schwertmannite_BZ93_1, y_copiapite_GDS21 ,y_ferrihydrite_GDS75_syn_f6))