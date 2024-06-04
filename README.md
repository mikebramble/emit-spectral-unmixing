# EMIT spectral unmixing 
 A personal application of spectral unmixing to the EMIT reflectance data set.

The key script is scripts_emit_spectral_unmixing, which calls on script_py_load_amd_endmembers.

The unmixing script is currently a minimum working model that loads and EMIT image, loads the unmixing endmembers, resamples the data sets, and then performs an unmixing using scipy.optimize.lsq_linear.

The remaining scripts are for testing purposes.