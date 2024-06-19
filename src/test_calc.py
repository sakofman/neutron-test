""" 
Test calculation going through the tasks provided 

TODO: The sequence of steps is written out, need to finish writing the processing module 
"""

import numpy as np
import pandas as pd

import processing_module as pm
import beam_module as b

""" Load constants and isotope data """
# proton_pulse_gap, flight_path, trigger_delay, time_bin, min_E, max_E, crs_skip_pts
constants = pm.CONST_DICT
print(constants)

# create isotope dataframe
isotope_df = pm.create_isotope_dataframe()
print("Isotope data:")
print(isotope_df.to_string())


# Load interpolated cross section data, between min_E and max_E (specified in the constants)
# Store in a dictionary with key as the cross-section name, and value as the interpolated function
interp_crs = pm.create_crs_dict(isotope_df)

# Create energy, time meshgrids with appropriate spacing
E_mesh = pm.create_E_mesh()
T_mesh = pm.create_T_mesh()

# Calculate ideal transmission on a properly spaced Energy mesh
ideal_tr = pm.ideal_transmission_solid(interp_crs, E_mesh)

# Convolve ideal Theoretical transmission with beam pulse profile
beam = b.Beam()
# convolve()

# Call fitting functinos below